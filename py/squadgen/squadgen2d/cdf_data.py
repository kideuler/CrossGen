"""
cdf_data.py -- SQ-VAE training data: CrossGen block decompositions as SQuadGen chart distance fields.

One task is one (model, method): run the method through the crossgen module, and if its blocks are a valid
decomposition, build their CDF/DCDF in C++ (BlockDecomposition.cdf(), src/mesh/ChartDistanceField.hxx) and write the
release's NPZ (docs/squadgen_2d.md Sec. 2.4), which data_tools/test_load_new_format.py loads unchanged. `run` runs
every task over directories of OBJs in parallel, each in its own process, then `manifest` gathers them:

    python -m squadgen2d cdf-data run                                     # singlemat + mambo, five methods -> data/cdf
    python -m squadgen2d cdf-data run --methods oracle --out py/squadgen/data/cdf_oracle
    python -m squadgen2d cdf-data run data/objs/2d_data --methods oracle --out py/squadgen/data/cdf_2d --jobs 8
    python -m squadgen2d cdf-data task --out DIR data/meshes/mambo/S39_12.obj zipline     # one task
    python -m squadgen2d cdf-data manifest --out DIR

This is the one command that needs the crossgen module (build it with `cmake --build build --target crossgen` and
run with PYTHONPATH=build/python, or configure once with -DCROSSGEN_PYTHON_DEV_PTH=ON; `run` also finds it in
build/python by itself). Everything after it -- cdf-pack and the trainers -- needs only numpy and torch, so the
directory it writes can be copied to a machine without CrossGen.

Layout of DIR:

    npz/<source>/<model>/<method>.npz     one layout's fields; beside the release's arrays, 0-d metadata: model,
                                          source, family, method, split, coverage, regularity, angle_quality,
                                          chord_quality, alignment_quality
    status/<source>/<model>/<method>.json what the task did, whether or not it wrote an NPZ
    manifest.csv                          one row per task, status files merged

A status is "ok" (NPZ written), "invalid" (the method's blocks do not tile the model), "folded" (the block mesh under
the fields has inverted elements; no NPZ unless --keep-folded), "error" (an exception, in `error`) or "running" (the
process died or hit --timeout before it finished). Each task runs in its own process so that a method that crashes
takes only its own task with it. `run` resumes: a task with a status file is skipped, so stopping and rerunning does
only what is left (tasks stopped by Ctrl-C are forgotten, not recorded as died). Measured on a 44-model sample
(2026-10-28): 87 % of tasks give a layout, ~1.2 MB each at the defaults, ~1.6 s each.

`split` is the train/val split a Geom-AE trained on the same OBJ directories would draw (domains.val_families,
val_frac 0.1, seed 0). The trainers use the split of the Geom-AE they actually freeze, recorded in its checkpoint;
this column is their fallback for a Geom-AE that recorded none.

One layout per domain (--fallback)
----------------------------------
`run --fallback zipline atlas` (the default chain when --fallback is given bare) gives each domain one layout instead
of one per method: ZIPLINE runs first, and its layout is kept when its blocks are a decomposition with 100 % coverage
and their CDF builds unfolded; otherwise ATLAS runs and is held to the same test; a domain neither passes is skipped,
printed as it happens with each method's reason, and listed in DIR/skipped.txt.

    python -m squadgen2d cdf-data run data/objs/2d_data --fallback zipline atlas    # -> data/cdf_zipline_atlas

The test is stricter than `valid`, which accepts 99.9 % coverage: a tiling's area ratio comes out at 1 - 2e-13 from
rounding, while a ZIPLINE layout with a face it could not make a block (a T-junction on a side) falls short by 1.3e-5
or more on the corpus, so anything short of 1 - 1e-9 is a hole the fields would have. A layout whose CDF folds (or
errors) falls through to the next method as well: the chain's purpose is a usable layout per domain, and a folded
one gives no NPZ. Each method is still its own task and process, with its own status file and --timeout, so a method
that crashes or hangs passes the domain on to the next instead of losing it, and a rerun resumes where it stopped.
The manifest then holds one row per method tried and at most one ok row per domain, the kept method in its `method`
column. Measured on 18 corpus domains (2026-10-30): the NPZs are byte-identical to the same methods' in data/cdf.

The meshing and smoothing the layout is judged on are the CDF's own: b.cdf() meshes the blocks with BlockQuadMesh
(ZIPLINE's Stage 6 mesher) and Winslow-smooths the blocks the transfinite blend folds, and "folded" is that mesh
failing. The method's own b.mesh(h) and TMOP (q.smooth()) are not run, since the fields are not read off them.

DIR/fallback.json records the chain. A directory holds one kind of run, every method's or one chain's: a fallback run
refuses a directory with another chain's statuses or a run of every method's, and a run of every method refuses a
fallback directory, since either would break "at most one layout per domain" or mix the two.
"""

import argparse
import concurrent.futures as cf
import csv
import glob
import json
import os
import signal
import subprocess
import sys
import threading
import time
import warnings
import zlib

import numpy as np

from .domains import family, obj_paths, val_families
from .layout import METHODS
from .paths import CORPUS, DATA, DEFAULT_CDF, PROJECT, REPO

QUALITIES = ("regularity", "angle_quality", "chord_quality", "alignment_quality")
DEFAULT_METHODS = ("zipline", "umber", "meridian", "torsion", "atlas")
FALLBACK = ("zipline", "atlas")
# "100 % coverage" to rounding: a tiling's area ratio comes out at 1 - 2e-13 (ATLAS, MERIDIAN on the corpus), and the
# smallest shortfall of a ZIPLINE layout with uncovered faces there is 1.3e-5 (M7_15).
FULL_COVERAGE = 1.0 - 1e-9


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
              "path": os.path.relpath(os.path.abspath(args.obj), REPO), "status": "running"}
    if args.split:
        status["split"] = args.split
    else:
        status["split"] = "val" if status["family"] in val_families(obj_paths(CORPUS)) else "train"
    write_status(status_path, status)
    # No Python handler: SIGALRM's default action ends the process even inside a C++ call, which leaves the "running"
    # status behind as the record of it.
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
        if not b.valid or (args.full_coverage and b.coverage < FULL_COVERAGE):
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
    chain = fallback_chain(args.out)
    if chain:
        fallback_summary(args.out, chain, rows)


def read_status(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def why_not(method, status):
    """Why one method of a fallback chain gave no layout, as the skip line prints it."""
    if status is None:
        return f"{method} wrote no status"
    s = status["status"]
    if s == "invalid":
        c = status.get("coverage", 0.0)
        cover = "100 %" if c >= FULL_COVERAGE else f"{np.floor(1e4 * c) / 100:.2f} %"   # 99.999 % is not "100.00"
        return f"{method} {'valid but' if status.get('valid') else 'invalid,'} coverage {cover}"
    if s == "folded":
        return f"{method} folded ({status.get('inverted_faces')} inverted faces in the CDF mesh)"
    if s == "error":
        return f"{method} error: {status.get('error')}"
    return f"{method} died (crashed or timed out)" if s == "running" else f"{method} {s}"


def fallback_chain(out):
    """The chain of the fallback run that made `out` (DIR/fallback.json), or None for a run of every method."""
    s = read_status(os.path.join(out, "fallback.json"))
    return s["methods"] if s else None


def fallback_summary(out, chain, rows):
    """Per domain: the method kept (the first of the chain whose task is ok), or why none was. Writes the skipped
    domains to DIR/skipped.txt; a domain whose chain was stopped before its end is unfinished, not skipped."""
    by = {}
    for r in rows:
        by.setdefault((r["source"], r["model"]), {})[r["method"]] = r
    kept, skipped, unfinished = {m: 0 for m in chain}, [], 0
    for (source, model), rs in sorted(by.items()):
        ok = [m for m in chain if rs.get(m, {}).get("status") == "ok"]
        if ok:
            kept[ok[0]] += 1
        elif all(m in rs for m in chain):
            path = next(iter(rs.values()))["path"]
            skipped.append(f"{source}/{model}\t{path}\t" + "; ".join(why_not(m, rs[m]) for m in chain))
        else:
            unfinished += 1
    with open(os.path.join(out, "skipped.txt"), "w") as f:
        f.write("".join(s + "\n" for s in skipped))
    print(f"per domain ({' then '.join(chain)}): " + ", ".join(f"{m} {n}" for m, n in kept.items())
          + f", skipped {len(skipped)} (listed in {os.path.join(out, 'skipped.txt')})"
          + (f", {unfinished} unfinished (run again to finish them)" if unfinished else ""))


def check_kind(out, chain):
    """A directory holds one kind of run: a fallback chain's (at most one layout per domain) or every method's. Refuse
    to mix them, or two chains, in one directory; record the chain of a new fallback directory."""
    have = fallback_chain(out)
    started = next(glob.iglob(os.path.join(out, "status", "*", "*", "*.json")), None) is not None
    if chain is None:
        if have:
            raise SystemExit(f"{out} holds a fallback run ({' then '.join(have)}): a run of every method needs its "
                             "own --out")
        return
    if started and have != list(chain):
        what = f"a fallback run of {' then '.join(have)}" if have else "a run of every method"
        raise SystemExit(f"{out} already holds {what}: a fallback run of {' then '.join(chain)} needs its own --out")
    os.makedirs(out, exist_ok=True)
    write_status(os.path.join(out, "fallback.json"), {"methods": list(chain), "full_coverage": FULL_COVERAGE})


def forget(out, status_path):
    """Remove a task's status and its NPZ: a later method of a chain left over from a --force rerun in which an
    earlier one now gives the layout."""
    s = read_status(status_path)
    if s and s.get("npz") and os.path.exists(os.path.join(out, s["npz"])):
        os.remove(os.path.join(out, s["npz"]))
    if os.path.exists(status_path):
        os.remove(status_path)


def crossgen_env():
    """The environment tasks run in: this one with the package importable, and build/python added when that is
    where the crossgen module is. None when no crossgen module with BlockDecomposition.cdf can be found."""
    env = dict(os.environ)
    paths = [PROJECT] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p]
    probe = [sys.executable, "-c", "import crossgen; crossgen.BlockDecomposition.cdf"]
    for extra in ([], [os.path.join(REPO, "build", "python")]):
        env["PYTHONPATH"] = os.pathsep.join(paths + extra)
        if subprocess.run(probe, env=env, capture_output=True).returncode == 0:
            return env
    return None


def run_all(args):
    """Every (OBJ, method) task under args.dirs, args.jobs at a time, then the manifest. Each task is its own process
    in its own session: a crash or a SIGALRM timeout takes only that task, and Ctrl-C here stops the tasks running
    and forgets their statuses, so that the next run does them again rather than counting them as died.

    With --fallback, a job is a domain rather than a task: its chain's tasks one after another, until one is ok."""
    env = crossgen_env()
    if env is None:
        raise SystemExit(f"{sys.executable} cannot import a crossgen module with BlockDecomposition.cdf. Build it:  "
                         "cmake --build build --target crossgen -j8  and run with PYTHONPATH=build/python (or "
                         "configure once with -DCROSSGEN_PYTHON_DEV_PTH=ON)")
    chain = None if args.fallback is None else (args.fallback or list(FALLBACK))
    methods = chain or args.methods or list(DEFAULT_METHODS)
    if args.out is None:
        args.out = os.path.join(DATA, "cdf_" + "_".join(chain)) if chain else DEFAULT_CDF
    objs = obj_paths(args.dirs)
    # Drawn over every domain of the directories, as a Geom-AE trained on them would, before --models picks some.
    val = val_families(objs)
    if args.models:
        objs = [p for p in objs if os.path.splitext(os.path.basename(p))[0] in args.models]
    if not objs:
        raise SystemExit(f"no OBJ files to run (under {' '.join(args.dirs)}" + (", --models)" if args.models else ")"))
    check_kind(args.out, chain)
    common = ["--out", args.out, "--h-rel", str(args.h_rel), "--samples", str(args.samples), "--fps-min",
              str(args.fps_min), "--fps-max", str(args.fps_max), "--copies", str(args.copies), "--timeout",
              str(args.timeout)] + (["--keep-folded"] if args.keep_folded else []) + (["--force"] if args.force else [])
    if chain:
        common.append("--full-coverage")

    # A job: (label, [(method, status path, command), ...]), run in order until a status is ok. Done already when its
    # steps' statuses exist up to an ok one or to the end (a task with a status is not rerun without --force).
    jobs = []
    for obj in objs:
        split = "val" if family(obj) in val else "train"
        source, model = task_paths(args.out, obj, methods[0])[:2]
        steps = [(m, task_paths(args.out, obj, m)[3], [sys.executable, "-m", "squadgen2d.cdf_data", "task"] + common
                  + ["--split", split, obj, m]) for m in methods]
        jobs += [(f"{source}/{model}", steps)] if chain else [(f"{source}/{model} {s[0]}", [s]) for s in steps]

    def settled(steps):
        for _, status, _ in steps:
            s = read_status(status)
            if s is None:
                return False
            if s["status"] == "ok":
                return True
        return True

    todo = [j for j in jobs if args.force or not settled(j[1])]
    os.makedirs(args.out, exist_ok=True)
    if chain:
        print(f"{len(objs)} models, each {' then '.join(chain)} until one gives a layout (valid, 100 % coverage, "
              f"unfolded), else skipped; {len(jobs) - len(todo)} done already; {len(todo)} to run, {args.jobs} at a "
              f"time -> {args.out}", flush=True)
    else:
        print(f"{len(objs)} models x {{{' '.join(methods)}}} = {len(jobs)} tasks, {len(jobs) - len(todo)} done "
              f"already; {len(todo)} to run, {args.jobs} at a time -> {args.out}", flush=True)

    # Ctrl-C sets `stopping` under the lock that launches go through, so that no thread starts a chain's next task
    # once the running ones have been collected for termination.
    running, lock, stopping = {}, threading.Lock(), threading.Event()

    def launch(status, cmd):
        with lock:
            if stopping.is_set():
                return None
            p = subprocess.Popen(cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                 start_new_session=True)
            running[status] = p
        try:
            return p.wait()
        finally:
            running.pop(status, None)

    def one(steps):
        """The job's tasks in order until one is ok: [(method, status), ...] for those reached, and how many died."""
        seen, died = [], 0
        for k, (method, status, cmd) in enumerate(steps):
            if args.force or not os.path.exists(status):
                rc = launch(status, cmd)
                if rc is None or stopping.is_set():
                    return seen, died
                died += rc != 0
            s = read_status(status)
            seen.append((method, s))
            if s is not None and s["status"] == "ok":
                if args.force:
                    for _, later, _ in steps[k + 1:]:
                        forget(args.out, later)
                break
        return seen, died

    t0, done, died, skipped = time.time(), 0, 0, 0
    pool = cf.ThreadPoolExecutor(args.jobs)
    try:
        futs = {pool.submit(one, steps): label for label, steps in todo}
        for f in cf.as_completed(futs):
            seen, d = f.result()
            done += 1
            died += d
            if chain and not any(s is not None and s["status"] == "ok" for _, s in seen):
                skipped += 1
                print(f"  skipped {futs[f]}: " + "; ".join(why_not(m, s) for m, s in seen), flush=True)
            if done % max(1, len(todo) // 20) == 0 and done < len(todo):
                el = time.time() - t0
                print(f"  {done}/{len(todo)} in {el / 60:.1f} min, ~{el / done * (len(todo) - done) / 60:.0f} min "
                      f"left ({died} tasks died" + (f", {skipped} skipped)" if chain else ")"), flush=True)
    except KeyboardInterrupt:
        with lock:
            stopping.set()
            stopped = dict(running)
        pool.shutdown(wait=False, cancel_futures=True)
        for p in stopped.values():
            p.terminate()
        for status, p in stopped.items():
            p.wait()
            s = read_status(status)
            if s is not None and s.get("status") == "running":
                os.remove(status)
        raise SystemExit(f"\nstopped after {done} {'domains' if chain else 'tasks'}; {len(stopped)} running tasks "
                         "forgotten. Run the same command again to continue.")
    pool.shutdown()
    print(f"done in {(time.time() - t0) / 60:.1f} min", flush=True)
    run_manifest(args)
    return 0


def add_cdf_options(p, out=DEFAULT_CDF):
    p.add_argument("--out", default=out, help="the data directory (default: py/squadgen/data/cdf; with --fallback, "
                   "py/squadgen/data/cdf_<chain>, e.g. cdf_zipline_atlas)")
    p.add_argument("--h-rel", type=float, default=0.01,
                   help="block-mesh edge length over the bounding-box diagonal (default 0.01)")
    p.add_argument("--samples", type=int, default=20000, help="area samples per layout")
    p.add_argument("--fps-min", type=int, default=512, help="token-set sizes from fps-min to fps-max, doubling")
    p.add_argument("--fps-max", type=int, default=4096)
    p.add_argument("--copies", type=int, default=5, help="token sets of each size")
    p.add_argument("--timeout", type=int, default=900, help="seconds per task; 0 for none")
    p.add_argument("--keep-folded", action="store_true", help="write layouts whose block mesh folds")
    p.add_argument("--force", action="store_true", help="redo tasks that already have a status (default: resume)")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="squadgen2d cdf-data", description=__doc__.split("\n")[1].strip())
    sub = ap.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run", help="every (model, method) under some OBJ directories, in parallel, then the manifest")
    r.add_argument("dirs", nargs="*", default=CORPUS,
                   help="directories of .obj domains, searched recursively (default: data/meshes/singlemat and mambo)")
    r.add_argument("--methods", nargs="+", choices=METHODS,
                   help="default: the five; oracle is ORACLE's choice among them, one layout per model")
    r.add_argument("--fallback", nargs="*", choices=METHODS, metavar="METHOD",
                   help="one layout per domain instead: these methods in order (bare: zipline atlas), keeping the "
                        "first whose blocks are valid with 100%% coverage and whose CDF is unfolded; a domain none "
                        "passes is skipped and printed")
    r.add_argument("--models", nargs="+", help="only these models (OBJ names without .obj)")
    r.add_argument("--jobs", type=int, default=os.cpu_count() or 4,
                   help="tasks (with --fallback: domains) at a time (default: logical CPUs)")
    add_cdf_options(r, out=None)
    t = sub.add_parser("task", help="one (model, method)")
    add_cdf_options(t)
    t.add_argument("--split", choices=("train", "val"),
                   help="the domain's split (default: drawn over data/meshes/singlemat + mambo)")
    t.add_argument("--full-coverage", action="store_true",
                   help="count valid blocks with less than 100%% coverage as invalid (--fallback's test)")
    t.add_argument("obj")
    t.add_argument("method", choices=METHODS)
    mf = sub.add_parser("manifest", help="merge the status files into manifest.csv")
    mf.add_argument("--out", default=DEFAULT_CDF)
    args = ap.parse_args(argv)
    if args.command == "run":
        if args.fallback is not None and args.methods:
            ap.error("--methods (every method's layout) and --fallback (one layout per domain) are two kinds of run")
        if args.fallback and len(set(args.fallback)) < len(args.fallback):
            ap.error(f"--fallback names a method twice: {' '.join(args.fallback)}")
        return run_all(args)
    return run_task(args) if args.command == "task" else run_manifest(args)


if __name__ == "__main__":
    raise SystemExit(main())
