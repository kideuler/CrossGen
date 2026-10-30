"""SQ-VAE training data: CrossGen block decompositions as SQuadGen chart distance fields.

One task is one (model, method): run the method through the crossgen module, and if its
blocks are a valid decomposition, build their CDF/DCDF in C++ (BlockDecomposition.cdf(),
src/mesh/ChartDistanceField.hxx) and write the release's NPZ (docs/squadgen_2d.md Sec. 2.4),
which data_tools/test_load_new_format.py loads unchanged. make_cdf_data.sh runs every task
over data/meshes/singlemat and data/meshes/mambo in parallel, then `manifest` gathers them.

    python3 py/squadgen/cdf_data.py task --out DIR data/meshes/mambo/S39_12.obj zipline
    python3 py/squadgen/cdf_data.py manifest --out DIR

Layout of DIR:

    npz/<source>/<model>/<method>.npz     one layout's fields; beside the release's arrays,
                                          0-d metadata: model, source, family, method, split,
                                          coverage, regularity, angle_quality, chord_quality,
                                          alignment_quality
    status/<source>/<model>/<method>.json what the task did, whether or not it wrote an NPZ
    manifest.csv                          one row per task, status files merged

A status is "ok" (NPZ written), "invalid" (the method's blocks do not tile the model),
"folded" (the block mesh under the fields has inverted elements; no NPZ unless
--keep-folded), "error" (an exception, in `error`) or "running" (the process died or hit
--timeout before it finished). Each task runs in its own process so that a method that
crashes takes only its own task with it.

`split` is the Geom-AE's train/val split by family (geom_ae_2d.split_by_family, val_frac
0.1, seed 0, over the same two folders), so an SQ-VAE validated on "val" never trains on
a shape family the frozen Geom-AE was validated on either.
"""
import argparse
import csv
import glob
import json
import os
import signal
import sys
import time
import warnings
import zlib

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CORPUS = [os.path.join(ROOT, "data", "meshes", d) for d in ("singlemat", "mambo")]
# "oracle" is ORACLE, the selector's choice among the other five (Mesh.oracle()): one layout per
# model, the kept method's own run, with which method it kept in `oracle_method`.
METHODS = ("zipline", "umber", "meridian", "torsion", "atlas", "oracle")
QUALITIES = ("regularity", "angle_quality", "chord_quality", "alignment_quality")


def family(path):
    """Mambo S39_12 -> 'S39'; singlemat geom012 -> 'geom012' (as geom_ae_2d.family)."""
    name = os.path.splitext(os.path.basename(path))[0]
    return name.rsplit("_", 1)[0] if "_" in name else name


def val_families(val_frac=0.1, seed=0):
    """geom_ae_2d.split_by_family's held-out families, over the corpus' OBJs."""
    fams = sorted({family(p) for d in CORPUS for p in glob.glob(os.path.join(d, "*.obj"))})
    rng = np.random.default_rng(seed)
    return set(rng.choice(fams, max(1, round(val_frac * len(fams))), replace=False))


def task_paths(out, obj, method):
    source = os.path.basename(os.path.dirname(os.path.abspath(obj)))
    model = os.path.splitext(os.path.basename(obj))[0]
    rel = os.path.join(source, model, method)
    return source, model, os.path.join(out, "npz", rel + ".npz"), os.path.join(out, "status", rel + ".json")


def write_status(path, status):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(status, f, indent=1)
    os.replace(tmp, path)


def run_task(args):
    source, model, npz_path, status_path = task_paths(args.out, args.obj, args.method)
    if os.path.exists(status_path) and not args.force:
        return 0
    status = {"model": model, "source": source, "family": family(args.obj), "method": args.method,
              "path": os.path.relpath(os.path.abspath(args.obj), ROOT), "status": "running"}
    status["split"] = "val" if status["family"] in val_families() else "train"
    write_status(status_path, status)
    # No Python handler: SIGALRM's default action ends the process even inside a C++ call,
    # which leaves the "running" status behind as the record of it.
    if args.timeout > 0:
        signal.alarm(args.timeout)

    import crossgen
    warnings.filterwarnings("ignore", "produced no blocks")   # the status says "invalid"
    try:
        t = time.time()
        m = crossgen.load(args.obj)
        b = getattr(m, args.method)()
        status.update(seconds_method=time.time() - t, num_blocks=b.num_blocks, coverage=b.coverage,
                      valid=b.valid)
        if args.method == "oracle":
            status.update(oracle_method=b.report["method"], oracle_decision=b.report["decision"])
        if not b.valid:
            status["status"] = "invalid"
            write_status(status_path, status)
            return 0
        for q in QUALITIES:
            status[q] = float(getattr(b, q)())

        t = time.time()
        # A seed per (model, method), so that a rerun redraws the same samples.
        seed = zlib.crc32(f"{source}/{model}/{args.method}".encode())
        c = b.cdf(relative_edge_length=args.h_rel, samples=args.samples, fps_min=args.fps_min,
                  fps_max=args.fps_max, fps_copies=args.copies, seed=seed)
        r = c.report
        status.update(seconds_cdf=time.time() - t, charts=r["blocks"], subcharts=r["subcharts"],
                      faces=r["faces"], target_edge_length=r["target"], inverted_faces=r["inverted_faces"],
                      cdf_valid=r["valid"])
        if not r["valid"]:
            raise RuntimeError(f"inconsistent fields: {r}")
        if r["inverted_faces"] > 0 and not args.keep_folded:
            status["status"] = "folded"
            write_status(status_path, status)
            return 0
        os.makedirs(os.path.dirname(npz_path), exist_ok=True)
        meta = {k: status[k] for k in ("model", "source", "family", "method", "split", "coverage", "oracle_method")
                + QUALITIES if k in status}
        c.write(npz_path + ".tmp.npz", **meta)
        os.replace(npz_path + ".tmp.npz", npz_path)
        status.update(status="ok", npz=os.path.relpath(npz_path, args.out), bytes=os.path.getsize(npz_path))
    except Exception as e:  # recorded, not raised: one model's failure is data about it
        status.update(status="error", error=f"{type(e).__name__}: {e}")
    write_status(status_path, status)
    return 0


def run_manifest(args):
    rows = []
    for p in sorted(glob.glob(os.path.join(args.out, "status", "*", "*", "*.json"))):
        with open(p) as f:
            rows.append(json.load(f))
    if not rows:
        sys.exit(f"no status files under {args.out}/status")
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(args.out, "manifest.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} tasks, {len({(r['source'], r['model']) for r in rows})} models -> {args.out}/manifest.csv")
    print(f"{'method':10s} {'ok':>5s} {'invalid':>8s} {'folded':>7s} {'error':>6s} {'died':>5s} {'charts':>7s}")
    for m in sorted({r["method"] for r in rows}, key=lambda m: METHODS.index(m) if m in METHODS else 99):
        rs = [r for r in rows if r["method"] == m]
        n = {s: sum(r["status"] == s for r in rs) for s in ("ok", "invalid", "folded", "error", "running")}
        charts = [r["charts"] for r in rs if r["status"] == "ok"]
        print(f"{m:10s} {n['ok']:5d} {n['invalid']:8d} {n['folded']:7d} {n['error']:6d} {n['running']:5d} "
              f"{np.median(charts) if charts else float('nan'):7.0f}")
    kept = [r["oracle_method"] for r in rows if r["status"] == "ok" and r.get("oracle_method")]
    if kept:
        print("oracle kept: " + ", ".join(f"{m} {kept.count(m)}" for m in METHODS if m in kept))
    ok = [r for r in rows if r["status"] == "ok"]
    print(f"NPZ: {len(ok)} layouts ({sum(r['split'] == 'train' for r in ok)} train, "
          f"{sum(r['split'] == 'val' for r in ok)} val), {sum(r['bytes'] for r in ok) / 2**30:.2f} GiB")
    for r in rows:
        if r["status"] in ("error", "running"):
            what = r.get("error", "crashed or timed out (status left at 'running')")
            print(f"  {'error' if r['status'] == 'error' else 'died':5s} {r['source']}/{r['model']} {r['method']}: {what}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    t = sub.add_parser("task", help="one (model, method)")
    t.add_argument("--out", required=True)
    t.add_argument("--h-rel", type=float, default=0.01,
                   help="block-mesh edge length over the bounding-box diagonal (default 0.01)")
    t.add_argument("--samples", type=int, default=20000)
    t.add_argument("--fps-min", type=int, default=512)
    t.add_argument("--fps-max", type=int, default=4096)
    t.add_argument("--copies", type=int, default=5)
    t.add_argument("--timeout", type=int, default=900, help="seconds; 0 for none")
    t.add_argument("--keep-folded", action="store_true", help="write layouts whose block mesh folds")
    t.add_argument("--force", action="store_true", help="redo a task that already has a status")
    t.add_argument("obj")
    t.add_argument("method", choices=METHODS)
    mf = sub.add_parser("manifest", help="merge the status files into manifest.csv")
    mf.add_argument("--out", required=True)
    args = ap.parse_args()
    return run_task(args) if args.command == "task" else run_manifest(args)


if __name__ == "__main__":
    sys.exit(main())
