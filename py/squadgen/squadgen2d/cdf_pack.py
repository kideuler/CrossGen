"""
cdf_pack.py -- presaved fields: every layout's CDF, DCDF, unit gradients, offsets and SDF at its token and query
points, computed once, so that SQ-VAE and SQ-Diffuse train from arrays instead of evaluating fields per sample.

A CDF NPZ (cdf_data.py) holds a layout's split quad mesh and coordinate pairs, the release's format: fields are
evaluated from it (layout.layout_fields, Appendix B). Training used to do that for every sample, reading and
rebuilding the NPZ's mesh each time (11 of a sample's 16 ms on the M2). A pack evaluates the fields once, at exactly
the points training reads -- the NPZ's N_FPS_COPIES token sets of N and 4N points and its area samples -- in the
Geom-AE's frame, and stores them as plain float32 arrays that the DataLoader memory-maps. A training sample is then a
few slices of them.

    python -m squadgen2d cdf-pack py/squadgen/data/cdf_oracle             # -> py/squadgen/data/cdf_oracle/pack_N512
    python -m squadgen2d cdf-pack py/squadgen/data/cdf --tokens 512 --jobs 8

A pack needs only numpy, not the crossgen module, so it can be built wherever the NPZs are (on a cloud box after
uploading them, say). It is rebuilt only when the manifest's layouts have changed (or with --force); the trainers
refuse a pack that no longer matches the manifest beside it.

Augmentation. The fields of a deformed layout are those of the same points of the deformed mesh, and that is what
make_sample() computes, from the points alone. A Deformation (geometry.py) is a smooth map W with a known Jacobian J,
so at a point p with corner c:
  * cdf and dcdf do not change. They are functions of the barycentric coordinates in the point's triangle
    (t0 = 1 - b0, t1 = b2 in layout_fields), which a deformation of the mesh keeps;
  * a unit gradient g becomes J(p)^-T g, renormalised, as the gradient of a field carried along by W does;
  * an offset c - p becomes W(c) - W(p), and the position becomes W(p);
  * then the frame is re-normalised from the deformed dOmega (geometry.boundary_frame), and sdf is measured to it.
The mesh augmentation applied W to the split mesh's vertices and evaluated the fields on the deformed piecewise-
linear mesh, so the two differ by the warp's curvature over one split-mesh triangle: O(h^2 |W''|) in positions and
offsets, ~1e-4 here, and O(h |W''|) in gradient directions. Without augmentation a sample is bit-identical to one
built from the NPZ (layout_fields at the same points, cast to float32 the same way), and the draws are the ones the
NPZ path made, so evaluation sets are the same layouts at the same points as before packs.

Layout of a pack directory:

    meta.json         format, version, tokens N, copies, tau, the channel tables below, the data directory
    layouts.csv       the manifest's status-ok rows (every column) plus npz_bytes, q0 q1 (rows in query.npy),
                      b0 b1 (rows in boundary.npy), centre_x centre_y scale (x_obj = x / scale + centre)
    tok{N}.npy        (layouts, copies, N, 23)  float32, TOKEN_FIELDS per point
    tok{4N}.npy       (layouts, copies, 4N, 23) float32
    query.npy         (all samples, 18) float32, QUERY_FIELDS per point; layout i's are rows q0:q1
    boundary.npy      (all boundary edges, 2, 2) float32, dOmega of the split mesh as oriented segments

At the default N = 512 and 20 000 samples a layout takes ~2.6 MB (the NPZ ~1.2 MB): ~1.1 GB for the 425 ORACLE
layouts.
"""

import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import os
import shutil
import time

import numpy as np

from . import geometry as geo
from .layout import OFFSET_KEYS, Layout, layout_fields, read_manifest
from .paths import DEFAULT_CDF

FORMAT, VERSION = "squadgen2d-field-pack", 1
TOKEN_FIELDS = (("xyz", 2), ("dcdf", 1), ("gdcdf", 2), ("cdf", 1), ("gcdf", 2), ("sdf", 1)) + tuple(
    (k, 2) for k in OFFSET_KEYS)
QUERY_FIELDS = (("xyz", 2), ("dcdf", 1), ("cdf", 1)) + tuple((k, 2) for k in OFFSET_KEYS)
GRAD_KEYS = ("gdcdf", "gcdf")
PACK_COLUMNS = ("npz_bytes", "q0", "q1", "b0", "b1", "centre_x", "centre_y", "scale")


def _slices(fields):
    out, i = {}, 0
    for k, d in fields:
        out[k] = slice(i, i + d)
        i += d
    return out, i


TOKEN, TOKEN_DIM = _slices(TOKEN_FIELDS)
QUERY, QUERY_DIM = _slices(QUERY_FIELDS)


def default_pack(data, n_tokens):
    return os.path.join(data, f"pack_N{n_tokens}")


def _join(f, fields):
    return np.concatenate([f[k].reshape(len(f["xyz"]), d) for k, d in fields], 1).astype(np.float32)


def _split(a, table):
    return {k: a[:, s] for k, s in table.items()}


def record_from_npz(path, n_tokens, tau):
    """One layout's pack arrays, from its NPZ: {"tok": {n: (copies, n, TOKEN_DIM)} for n = N, 4N, "query":
    (samples, QUERY_DIM), "bnd": (edges, 2, 2), "tau": tau}, and the Layout they came from."""
    sizes = (n_tokens, 4 * n_tokens)
    L = Layout(path, sizes)
    A, B = L.V[L.E[:, 0]], L.V[L.E[:, 1]]
    tok = {}
    for n in sizes:
        tri, bary = L.fps[n]
        out = np.empty((len(tri), n, TOKEN_DIM), np.float32)
        for c in range(len(tri)):
            f = layout_fields(L, L.V, tri[c], bary[c])
            f["sdf"] = geo.truncated_sdf(f["xyz"], A, B, tau, known_inside=True)
            out[c] = _join(f, TOKEN_FIELDS)
        tok[n] = out
    query = _join(layout_fields(L, L.V, *L.samples), QUERY_FIELDS)
    return {"tok": tok, "query": query, "bnd": np.stack([A, B], 1).astype(np.float32), "tau": tau}, L


def _deform(f, W, centre, scale, grads):
    """The fields f (float64, one point set) carried along by the deformation W, in the frame (centre, scale) of the
    deformed domain. See the module docstring."""
    p = f["xyz"]
    Wp = W(p)
    corners = W(np.concatenate([p + f[k] for k in OFFSET_KEYS]))
    for i, k in enumerate(OFFSET_KEYS):
        f[k] = (corners[i * len(p):(i + 1) * len(p)] - Wp) * scale
    if grads:
        J = W.jacobian(p)
        # J^-T g = (1 / det J) [[d, -c], [-b, a]] g for J = [[a, b], [c, d]]; the sign of det J matters, its size not.
        sgn = np.sign(J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0])[:, None]
        for k in GRAD_KEYS:
            g = f[k]
            h = sgn * np.stack([J[:, 1, 1] * g[:, 0] - J[:, 1, 0] * g[:, 1],
                                -J[:, 0, 1] * g[:, 0] + J[:, 0, 0] * g[:, 1]], 1)
            f[k] = h / (np.linalg.norm(h, axis=1, keepdims=True) + 1e-10)
    f["xyz"] = (Wp - centre) * scale


def make_sample(rec, cfg, rng, aug, warp_amp=0.06):
    """One example under the release's batch keys (load_patch_data_reformat, 2D), from a layout's pack arrays.
    Encoder input : xyz_fps_{N} and <k>_fps_{N} for every token field (sdf, the fields, their unit gradients, the
                    seven offsets), the same at 4N, one random token set of each size out of the copies.
    Targets       : xyz_query and <k>_query for every other query field (dcdf, cdf, the offsets), at n_query of the
                    layout's area samples.
    With aug, a Deformation first. One that folds the domain (det J changing sign over the layout's 4N tokens and
    dOmega) is rejected and another drawn, up to five times, as geometry.augment() rejects folded meshes."""
    n, tau = cfg["n_tokens"], cfg["tau"]
    bnd = rec["bnd"].astype(np.float64)
    A, B = bnd[:, 0], bnd[:, 1]
    W = None
    if aug:
        probe = np.concatenate([rec["tok"][4 * n][0][:, TOKEN["xyz"]].astype(np.float64), A])
        for _ in range(5):
            d = geo.Deformation(rng, warp_amp)
            J = d.jacobian(probe)
            det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
            if np.all(det > 0) or np.all(det < 0):
                W = d
                break
        if W is not None:
            A, B = W(A), W(B)
            centre, scale = geo.boundary_frame(A, B)
            A, B = (A - centre) * scale, (B - centre) * scale
    out = {}
    for m in (n, 4 * n):
        t = rec["tok"][m]
        f = _split(t[rng.integers(len(t))].astype(np.float64), TOKEN)
        if W is not None:
            _deform(f, W, centre, scale, grads=True)
        if W is not None or rec["tau"] != tau:
            f["sdf"] = geo.truncated_sdf(f["xyz"], A, B, tau, known_inside=True)
        out.update({f"{k}_fps_{m}": v for k, v in f.items()})
    q = rec["query"]
    idx = rng.choice(len(q), cfg["n_query"], replace=cfg["n_query"] > len(q))
    f = _split(q[idx].astype(np.float64), QUERY)
    if W is not None:
        _deform(f, W, centre, scale, grads=False)
    out.update({("xyz_query" if k == "xyz" else f"{k}_query"): v for k, v in f.items()})
    return {k: v.astype(np.float32) for k, v in out.items()}


# =============================================================================
# READING
# =============================================================================


class FieldPack:
    """A pack directory, read-only. `rows` are its layouts.csv as dicts (the manifest's columns as strings, the pack's
    as numbers), with `npz` absolute and `path` the domain OBJ as the manifest has it (repository-relative).
    record(i) is layout i's arrays as memory-mapped views. The maps are opened on first use in each process, so a
    pack handed to spawned DataLoader workers pickles as its rows alone."""

    def __init__(self, path):
        self.path = os.path.abspath(path)
        meta = os.path.join(self.path, "meta.json")
        if not os.path.exists(meta):
            raise SystemExit(f"no field pack at {self.path}")
        with open(meta) as f:
            self.meta = json.load(f)
        if self.meta.get("format") != FORMAT or self.meta.get("version") != VERSION:
            raise SystemExit(f"{self.path} is not a version-{VERSION} field pack: rebuild it with cdf-pack --force")
        self.n_tokens, self.tau = self.meta["tokens"], self.meta["tau"]
        self.data = os.path.normpath(os.path.join(self.path, self.meta["data"]))
        with open(os.path.join(self.path, "layouts.csv")) as f:
            self.rows = list(csv.DictReader(f))
        for r in self.rows:
            for k in PACK_COLUMNS:
                r[k] = float(r[k]) if k in ("centre_x", "centre_y", "scale") else int(r[k])
            r["npz"] = os.path.join(self.data, r["npz"])
        self._maps = None

    def __len__(self):
        return len(self.rows)

    def __getstate__(self):
        return dict(self.__dict__, _maps=None)

    def _open(self):
        if self._maps is None:
            n = self.n_tokens
            load = lambda name: np.load(os.path.join(self.path, name), mmap_mode="r")
            self._maps = {"tok": {n: load(f"tok{n}.npy"), 4 * n: load(f"tok{4 * n}.npy")},
                          "query": load("query.npy"), "bnd": load("boundary.npy")}
        return self._maps

    def record(self, i):
        r, m = self.rows[i], self._open()
        return {"tok": {k: v[i] for k, v in m["tok"].items()}, "query": m["query"][r["q0"]:r["q1"]],
                "bnd": m["bnd"][r["b0"]:r["b1"]], "tau": self.tau}

    def stale(self):
        """Why the pack no longer matches the manifest of its data directory, or None (also when there is no
        manifest there: a pack copied somewhere on its own)."""
        if not os.path.exists(os.path.join(self.data, "manifest.csv")):
            return None
        want = {(r["npz"], int(r["bytes"])) for r in read_manifest(self.data)}
        have = {(r["npz"], r["npz_bytes"]) for r in self.rows}
        if want == have:
            return None
        return (f"{len(want - have)} layouts of {self.data}/manifest.csv are new or changed since the pack was built, "
                f"and {len(have - want)} of the pack's are gone or changed")


def open_pack(data, pack, n_tokens):
    """The trainers' pack: `pack` if given, else DATA/pack_N{n_tokens}; it must exist, match the Geom-AE's token
    count and be up to date with the manifest."""
    path = pack or default_pack(data, n_tokens)
    if not os.path.exists(os.path.join(path, "meta.json")):
        raise SystemExit(f"no field pack at {path}: build it with  python -m squadgen2d cdf-pack {data} "
                         f"--tokens {n_tokens}")
    p = FieldPack(path)
    if p.n_tokens != n_tokens:
        raise SystemExit(f"{path} has {p.n_tokens} tokens, the Geom-AE {n_tokens}: build a pack with --tokens "
                         f"{n_tokens}")
    why = p.stale()
    if why:
        raise SystemExit(f"{path} is out of date ({why}): rebuild it with  python -m squadgen2d cdf-pack {p.data}")
    return p


# =============================================================================
# BUILDING
# =============================================================================


_files = None


def _open_files(files):
    # Unbuffered: a pool worker may end without flushing Python's buffers, and each write is a whole layout anyway.
    global _files
    _files = {k: (open(p, "r+b", buffering=0), off) for k, (p, off) in files.items()}


def _write(name, row_bytes, row, arr):
    f, off = _files[name]
    f.seek(off + row * row_bytes)
    f.write(np.ascontiguousarray(arr, np.float32).tobytes())


def _pack_job(job):
    """Compute one layout and write its token and query rows straight into the pack's files; the boundary, which
    has no fixed size, goes back to the parent."""
    i, npz, q0, n_tokens, tau = job
    rec, L = record_from_npz(npz, n_tokens, tau)
    for n, t in rec["tok"].items():
        _write(f"tok{n}", t[0].nbytes * len(t), i, t)
    _write("query", QUERY_DIM * 4, q0, rec["query"])
    return i, rec["bnd"], L.centre, L.scale


def manifest_digest(data):
    with open(os.path.join(data, "manifest.csv"), "rb") as f:
        return hashlib.sha1(f.read()).hexdigest()[:16]


def build(data, out, n_tokens, tau, jobs, log=print):
    rows = read_manifest(data)
    sizes = (n_tokens, 4 * n_tokens)
    t0 = time.time()
    counts, copies = [], None
    for r in rows:
        with np.load(r["npz"]) as z:
            missing = [f"fps_point_tri_id_{n}" for n in sizes if f"fps_point_tri_id_{n}" not in z.files]
            if missing:
                raise SystemExit(f"{r['npz']} has no token sets of {' '.join(missing)}: --tokens must be one of the "
                                 "sizes cdf_data.py drew (--fps-min .. --fps-max)")
            c = z[f"fps_point_tri_id_{n_tokens}"].shape[0]
            if copies is not None and c != copies:
                raise SystemExit(f"{r['npz']} has {c} token sets per size, others {copies}")
            copies = c
            counts.append(z["sample_point_tri_id"].shape[0])
    q1 = np.cumsum(counts)
    q0 = q1 - np.asarray(counts)
    tmp = out.rstrip(os.sep) + ".tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    files = {}
    for name, shape in [(f"tok{n}", (len(rows), copies, n, TOKEN_DIM)) for n in sizes] + [
            ("query", (int(q1[-1]), QUERY_DIM))]:
        a = np.lib.format.open_memmap(os.path.join(tmp, f"{name}.npy"), "w+", np.float32, shape)
        files[name] = (a.filename, a.offset)
        del a
    log(f"packing {len(rows)} layouts of {data}: {copies} x ({n_tokens} + {4 * n_tokens}) tokens and "
        f"{int(np.median(counts))} samples each, {jobs} processes ...")
    jobs_ = [(i, r["npz"], int(q0[i]), n_tokens, tau) for i, r in enumerate(rows)]
    bnd, frames = [None] * len(rows), [None] * len(rows)
    pool = mp.get_context("spawn").Pool(jobs, _open_files, (files,)) if jobs > 1 else None
    if pool is None:
        _open_files(files)
    try:
        it = pool.imap_unordered(_pack_job, jobs_, chunksize=4) if pool else map(_pack_job, jobs_)
        for k, (i, b, c, s) in enumerate(it, 1):
            bnd[i], frames[i] = b, (c, s)
            if k % 500 == 0 and k < len(rows):
                el = time.time() - t0
                log(f"  {k}/{len(rows)} in {el:.0f} s, ~{el / k * (len(rows) - k):.0f} s left")
    finally:
        if pool:
            pool.close()
            pool.join()
        else:
            for f, _ in _files.values():
                f.close()
    nb = np.array([len(b) for b in bnd])
    b1 = np.cumsum(nb)
    np.save(os.path.join(tmp, "boundary.npy"), np.concatenate(bnd).astype(np.float32))
    keys = list(rows[0]) + [k for k in PACK_COLUMNS if k not in rows[0]]
    with open(os.path.join(tmp, "layouts.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for i, r in enumerate(rows):
            w.writerow(dict(r, npz=os.path.relpath(r["npz"], data), npz_bytes=int(r["bytes"]), q0=int(q0[i]),
                            q1=int(q1[i]), b0=int(b1[i] - nb[i]), b1=int(b1[i]), centre_x=repr(float(frames[i][0][0])),
                            centre_y=repr(float(frames[i][0][1])), scale=repr(float(frames[i][1]))))
    meta = {"format": FORMAT, "version": VERSION, "tokens": n_tokens, "sizes": list(sizes), "copies": copies,
            "tau": tau, "layouts": len(rows), "samples": int(q1[-1]),
            "data": os.path.relpath(os.path.abspath(data), os.path.abspath(out)),
            "manifest_sha1": manifest_digest(data), "created": time.strftime("%Y-%m-%d %H:%M:%S"),
            "token_fields": [list(x) for x in TOKEN_FIELDS], "query_fields": [list(x) for x in QUERY_FIELDS]}
    with open(os.path.join(tmp, "meta.json"), "w") as f:
        json.dump(meta, f, indent=1)
    shutil.rmtree(out, ignore_errors=True)
    os.replace(tmp, out)
    size = sum(os.path.getsize(os.path.join(out, x)) for x in os.listdir(out))
    log(f"done in {time.time() - t0:.0f} s -> {out} ({size / 2 ** 30:.2f} GiB)")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="squadgen2d cdf-pack", description=__doc__.split("\n")[1].strip())
    ap.add_argument("data", nargs="?", default=DEFAULT_CDF, help="a cdf-data output directory (manifest.csv + npz/)")
    ap.add_argument("--out", help="the pack directory (default: DATA/pack_N<tokens>)")
    ap.add_argument("--tokens", type=int, default=512, help="N, the Geom-AE's num_latents; 4N keys/values go too")
    ap.add_argument("--tau", type=float, default=0.1, help="the Geom-AE's SDF truncation (others recompute sdf)")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--force", action="store_true", help="rebuild even if the pack is up to date")
    args = ap.parse_args(argv)
    out = os.path.abspath(args.out or default_pack(args.data, args.tokens))
    if not args.force and os.path.exists(os.path.join(out, "meta.json")):
        p = FieldPack(out)
        if p.stale() is None and (p.n_tokens, p.tau) == (args.tokens, args.tau):
            print(f"{out} is up to date ({len(p)} layouts); --force rebuilds it")
            return 0
    build(os.path.abspath(args.data), out, args.tokens, args.tau, max(1, args.jobs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
