"""
The Geom-AE's domains: OBJ files as normalised meshes, their families and train/val split, and the token-set cache.

Mambo names are <model>_<face> (S39_12 is planar face 12 of CAD model S39), DeepCAD names <model>_<profile>; a family
is the model, and the split holds out whole families, so validation measures generalisation to shapes from a model
the network has not seen in any face. SQ-VAE and SQ-Diffuse validate on the families their frozen Geom-AE held out
(geom_ae_val_families()), so that no network in the chain was trained on a validation domain.

Token sets. Like the release's NPZs (Sec. 2.4, `fps_bary_coords_{n}`), each domain stores N_FPS_COPIES
farthest-point-sampled token sets per size as (triangle, barycentric) pairs. They are computed once, in the
normalised frame, and cached in the Geom-AE's run directory (fps_cache_N{N}.npz). Because they are barycentric they
follow the mesh through any deformation.
"""

import glob
import os
import re
import time
import zlib

import numpy as np

from . import geometry as geo
from .paths import from_repo, to_repo

N_FPS_COPIES = 5        # token sets per domain and size, as the release's NPZs store


def family(path):
    """Mambo S39_12 -> 'S39'; singlemat geom012 -> 'geom012' (one family per model)."""
    name = os.path.splitext(os.path.basename(path))[0]
    return name.rsplit("_", 1)[0] if "_" in name else name


def val_families(paths, val_frac=0.1, seed=0):
    """The families held out of the domains at `paths`: a seeded draw of round(val_frac x families), at least one."""
    fams = sorted({family(p) for p in paths})
    rng = np.random.default_rng(seed)
    return {str(f) for f in rng.choice(fams, max(1, round(val_frac * len(fams))), replace=False)}


def split_by_family(domains, val_frac, seed):
    val = val_families([d.path for d in domains], val_frac, seed)
    return [d for d in domains if family(d.path) not in val], [d for d in domains if family(d.path) in val]


def obj_paths(dirs):
    """The .obj files under dirs. Searched recursively, so a make_mesh_archive.py directory (singlemat/, mambo/,
    deepcad/ inside it) can be passed whole; a set of normalised paths, so a directory given beside its own parent
    counts each domain once."""
    return sorted({os.path.normpath(p) for d in dirs
                   for p in glob.glob(os.path.join(d, "**", "*.obj"), recursive=True)})


def geom_ae_val_families(gcfg):
    """The families a Geom-AE held out, from its config, and where they were found: recorded in the checkpoint
    (since 2026-10-31), else the split its run logged ("N held-out families: ..." in log.txt beside it, the last run
    of that directory), else drawn again from its --data, if those directories are here. (None, why) when none of
    these is available."""
    if gcfg.get("val_families") is not None:
        return set(gcfg["val_families"]), "recorded in the Geom-AE checkpoint"
    log = os.path.join(from_repo(gcfg["out"]), "log.txt") if gcfg.get("out") else None
    if log and os.path.exists(log):
        fams = None
        with open(log) as f:
            for line in f:
                m = re.search(r"held-out families: ([^)]*)\)", line)
                if m:
                    fams = m.group(1).split()
        if fams:
            return set(fams), f"logged in {to_repo(log)}"
    dirs = [from_repo(d) for d in gcfg.get("data") or []]
    if dirs and all(os.path.isdir(d) for d in dirs):
        paths = obj_paths(dirs)
        if paths:
            return val_families(paths, gcfg["val_frac"], gcfg["seed"]), "drawn again from the Geom-AE's --data"
    return None, "not recorded, not logged, and its --data is not here"


class Domain:
    """One OBJ in the normalised frame. `centre`/`scale` map it back: x_obj = x / scale + centre, the transform Sec. 4
    step 1 says to keep."""

    def __init__(self, path=None, V=None, Fc=None):
        self.path = path
        self.name = os.path.basename(path) if path else "domain"
        if path is not None:
            V, Fc = geo.read_obj(path)
        self.F = Fc.astype(np.int32)       # every DataLoader worker holds a copy of the corpus
        self.E = geo.boundary_edges(self.F)
        self.bverts = np.unique(self.E)
        self.centre, self.scale = geo.normalisation(V, Fc, self.bverts)
        self.V = (V - self.centre) * self.scale
        self.cdf = np.cumsum(np.abs(geo.signed_areas(self.V, self.F)))
        self.fps = None      # {n: ((copies, n) triangle ids, (copies, n, 3) barycentrics)}


def fps_sets(d, sizes, rng):
    """N_FPS_COPIES token sets of each size for one domain: FPS on 8n (at least 8192) area samples. The release draws
    them from 50 000 Poisson-disk samples (Sec. 2.4)."""
    out = {}
    for n in sizes:
        tris, barys = [], []
        for _ in range(N_FPS_COPIES):
            tri, bary = geo.sample_area(d.V, d.F, d.cdf, max(8 * n, 8192), rng)
            sel = geo.farthest_point_sampling(geo.point_from_bary(d.V, d.F, tri, bary), n, rng)
            tris.append(tri[sel])
            barys.append(bary[sel])
        out[n] = (np.stack(tris).astype(np.int32), np.stack(barys).astype(np.float32))
    return out


def _fps_job(job):
    path, sizes = job
    return path, fps_sets(Domain(path), sizes, np.random.default_rng(zlib.crc32(os.path.basename(path).encode())))


def default_fps_workers(dev):
    """Processes for the one-off token-set cache when --fps_workers is not given. It is ~2 s of numpy per domain
    (N_FPS_COPIES greedy FPS runs per size), so 10 k DeepCAD domains are an hour on 6 processes -- the cap that suits
    an 8 GB M2, where each spawned worker also imports torch. A CUDA box has the cores and memory to spare and its
    GPU idles until the cache is done, so it gets every core it may use (the affinity mask: a container's share, not
    the host's count)."""
    cores = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    return cores if dev.startswith("cuda") else min(6, cores)


def load_domains(dirs, n_tokens, cache_dir, workers):
    """Load every OBJ and attach its cached FPS token sets, computing them on first use. A cache entry is keyed by
    file name and face count. A regenerated mesh with a different face count gets new sets, because its triangle ids
    would point elsewhere."""
    paths = obj_paths(dirs)
    if not paths:
        raise SystemExit(f"no .obj files in {dirs}")
    names = [os.path.basename(p) for p in paths]
    if len(set(names)) != len(names):
        raise SystemExit("two domains share a file name; the FPS cache is keyed by it")
    domains = [Domain(p) for p in paths]
    sizes = (n_tokens, 4 * n_tokens)
    cache = os.path.join(cache_dir, f"fps_cache_N{n_tokens}.npz")
    stored = dict(np.load(cache)) if os.path.exists(cache) else {}
    todo = [d.path for d in domains
            if f"{d.name}/nf" not in stored or int(stored[f"{d.name}/nf"]) != len(d.F)]
    if todo:
        print(f"computing FPS token sets for {len(todo)} domains ({workers} workers) ...", flush=True)
        t0 = time.time()
        jobs = [(p, sizes) for p in todo]
        results = []
        # imap, not map, so that a run over thousands of domains says how far it is; the cache is written only once
        # all are done, so stopping part way loses them all.
        pool = None
        if workers > 1:
            import multiprocessing as mp
            pool = mp.get_context("spawn").Pool(workers)
        try:
            for k, r in enumerate(pool.imap(_fps_job, jobs, chunksize=4) if pool else map(_fps_job, jobs), 1):
                results.append(r)
                if k % 500 == 0 and k < len(jobs):
                    el = time.time() - t0
                    print(f"  {k}/{len(jobs)} domains in {el:.0f} s, ~{el / k * (len(jobs) - k) / 60:.0f} min left",
                          flush=True)
        finally:
            if pool:
                pool.terminate()
        nf = {d.path: len(d.F) for d in domains}
        for p, out in results:
            name = os.path.basename(p)
            stored[f"{name}/nf"] = np.int64(nf[p])
            for n, (tri, bary) in out.items():
                stored[f"{name}/tri{n}"], stored[f"{name}/bary{n}"] = tri, bary
        os.makedirs(cache_dir, exist_ok=True)
        np.savez(cache, **stored)
        print(f"  done in {time.time() - t0:.0f} s -> {cache}", flush=True)
    for d in domains:
        d.fps = {n: (stored[f"{d.name}/tri{n}"], stored[f"{d.name}/bary{n}"]) for n in sizes}
    return domains
