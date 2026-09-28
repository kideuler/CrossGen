#!/usr/bin/env python3
"""Build the algorithm-selection dataset: boundary features per mesh (and, with
--eigs, a Shape-DNA spectrum), and for each of the five methods whether it
produced a block decomposition and how good that decomposition is.

    PYTHONPATH=build/python python3 py/build_dataset.py [--out dataset] [--jobs 4]
        [--eigs 0] [--copies 4] [--sources singlemat mambo] [--limit N]

Writes, in --out:
  dataset.csv   one row per mesh: mesh_id, group, source, area, lam_1..lam_K
                (with --eigs K), feat_<name> for each boundary feature,
                probe_<PROBE>_<column> (below), then for each method
                <method>_valid, _regularity, _angle_quality, _chord_quality,
                _coverage, _num_blocks and _seconds.
  spectra.csv   mesh_id, group, lam_1..lam_K         } the two files memo_data.py
  results.csv   mesh_id, algorithm, metric, value    } writes, for train_selector.py
  cache.jsonl   every task's raw result, one line per mesh; a rerun skips the
                tasks already in it, so an interrupted run resumes, a run that
                asks for less (no spectrum, fewer eigenvalues) only rewrites the
                CSVs, and one that asks for more runs only what is missing.
  logs/, tmp/   each mesh's C++ output and its workers' result files while it
                runs; both are deleted when the script finishes.

Every number is written exactly as the C++ side returns it, averaged over the
copies below; nothing is rescaled here.

The spectrum is off unless --eigs asks for it. On the 2026-09-28 data (442
meshes) it made the selector worse -- regret with a fallback, 5-fold CV by
shape, 0.219 with it and 0.186 without -- and it was a third of the build's
CPU time: a median 11 s a mesh, more than any one method's median. The area
column comes from the boundary features, so it is there either way.

The spectrum is ShapeDNA with normalization="weyl_ratio": lambda_k |Omega| /
(4 pi k), the area-normalised eigenvalue over Weyl's leading term. Left
area-normalised only, lam_100 is ~100x lam_1 on every shape and the columns
mostly encode k; divided by Weyl's line they sit near 1 and differ between
shapes by what the leading term does not know -- the boundary term (~ perimeter
/ sqrt(area k), the Dirichlet correction) and below. SPECTRUM_NORMALIZATION
takes any other ShapeDNA normalisation ("area", "none", ...) instead.

The boundary features are Mesh.boundary_features() (mesh::BoundaryFeatures):
the corners by their ideal block count, the holes, T3's singularity bound, the
fewest quarter turns any quad layout of the model can have, and the like -- what
a layout's topology is decided by, and what a spectrum resolves least. On the
2026-09-26 data they predicted the right method better than the spectrum did.

PROBE's first copy is also written on its own, as probe_<PROBE>_valid,
_regularity, _angle_quality, _chord_quality, _coverage and _num_blocks. The
selector runs that method first (ZIPLINE: ~0.1 s a mesh) and reads its result
as inputs, so in training those have to be one run, as they are in use, and
not the mean over copies. The metrics are empty where the run was not valid;
num_blocks is the run's own count, valid or not (the <method>_num_blocks columns
average the valid copies only); a run that raised has coverage 0 and no count.
On the 2026-09-28 data this input halved the selector's regret: a cross-field
method's own layout says whether the others will split a face into many
blocks, which nothing in the boundary features does.

The metrics are BlockDecomposition's regularity(), angleQuality() and
chordQuality(): in [0, 1], higher better, defined only for a decomposition that
covers the model. Whether it does is <method>_valid, a separate column, and a
method that raises, crashes the process or makes no progress for METHOD_TIMEOUT
is not valid either. Where no run of a method was valid its metrics are empty
(NaN), not 0: a zero there used to repeat the failure in every quality column
and was most of what separated the methods.

Each method runs on --copies copies of the mesh, each an exact symmetry of the
square applied to the vertices (a quarter turn, a mirror, ...; nothing
rounds), and ZIPLINE and ATLAS also get a different seed on each. On the same
shape the methods do not always agree with themselves: exact rotations and
mirrors changed ATLAS's result on 9 of 16 test faces and UMBER's on 6, and the
metrics note asks for three seeds and a median for that reason. A mesh's
<method>_valid is then the fraction of copies that were valid and its metrics
the mean over those, a label a shape can be expected to predict rather than
one draw of it. Vertices no triangle uses are dropped before anything runs.

Each mesh runs in its own child process (this script with --worker), so a
segfault or a hang in one method costs that method and not the run, and the
parent retries the tasks after it in a fresh child. Children are single-
threaded, and --jobs of them run at once.
"""
import argparse
import concurrent.futures as cf
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

K_EIGS = 0                         # the length of the spectrum; 0 = no spectrum (see above)
SPECTRUM_NORMALIZATION = "weyl_ratio"   # a ShapeDNA normalisation; see crossgen.Mesh.shape_dna
SHAPE_DNA_OPTIONS = dict(degree=3, boundary="dirichlet")
SOURCES = {"singlemat": REPO / "data/meshes/singlemat", "mambo": REPO / "data/meshes/mambo"}
METHODS = ["zipline", "umber", "meridian", "torsion", "atlas"]
METRICS = ["regularity", "angle_quality", "chord_quality"]   # BlockDecomposition methods = CSV names
# Mesh.boundary_features() keys written as feat_<key>: the ones that do not
# change with the model's scale (area and perimeter do).
FEATURES = ["regions", "holes", "euler", "corners", "corners_one_block", "corners_two_blocks",
            "corners_three_blocks", "corners_four_blocks", "acute_corners", "ambiguous_corners",
            "corner_defect", "singularity_bound", "minimum_defect", "isoperimetric_ratio",
            "curved_fraction", "shortest_run", "interface_length"]
# The method the selector runs first and reads as inputs (None: no probe
# columns); its copy 0 is written on its own as probe_<PROBE>_<column>.
PROBE = "zipline"
PROBE_COLUMNS = ["valid", *METRICS, "coverage", "num_blocks"]
COPIES = 4                        # runs of each method per mesh, one per symmetry below
# Symmetries of the square, as the rows of a 2x2 matrix: the first COPIES of
# them are applied. The first four cover both axis orders in both handednesses.
SYMMETRIES = [((1, 0), (0, 1)),     # as given
              ((0, -1), (1, 0)),    # a quarter turn
              ((-1, 0), (0, 1)),    # mirrored in x
              ((0, 1), (1, 0)),     # mirrored across the diagonal
              ((-1, 0), (0, -1)),   # a half turn
              ((0, 1), (-1, 0)),    # three quarter turns
              ((1, 0), (0, -1)),    # mirrored in y
              ((0, -1), (-1, 0))]   # mirrored across the other diagonal
SEED_OPTIONS = {"zipline": "field_seed", "atlas": "anneal_seed"}   # copy k adds k to the default
CACHE_VERSION = 3                  # records of any other version are ignored (and rerun)
METHOD_TIMEOUT = 300.0            # seconds without progress before a child is killed
SPECTRUM_TIMEOUT = 60.0


# ── the child: one mesh, a list of tasks ────────────────────────────────────
def import_crossgen():
    try:
        import crossgen
    except ImportError:
        sys.path.insert(0, str(REPO / "build/python"))
        import crossgen
    return crossgen


def worker(mesh_path, tasks, k, result_path):
    """Run `tasks` ("spectrum", "features" and "<method>#<copy>") on one mesh,
    appending one JSON line per task to `result_path` -- a "start" line before
    it, so the parent knows which task a crash or hang belongs to, and the
    result after."""
    import numpy as np
    crossgen = import_crossgen()
    m = crossgen.load(mesh_path)
    # Only the vertices some triangle uses: mesh::Mesh keeps every vertex of
    # the file, and ones in no triangle change what UMBER and ATLAS compute.
    used = np.unique(m.triangles)
    if len(used) < m.num_vertices:
        remap = np.full(m.num_vertices, -1, dtype=np.int64)
        remap[used] = np.arange(len(used))
        m = crossgen.Mesh(m.vertices[used], remap[m.triangles], m.materials)
    copies = {}

    def copy_of(c):
        if c not in copies:
            R = np.array(SYMMETRIES[c], dtype=float)
            copies[c] = m if c == 0 else crossgen.Mesh(m.vertices @ R.T, m.triangles, m.materials)
        return copies[c]

    with open(result_path, "a") as out:
        def emit(**rec):
            out.write(json.dumps(rec) + "\n")
            out.flush()

        for task in tasks:
            emit(task=task, status="start")
            t0 = time.time()
            try:
                if task == "spectrum":
                    lam, rep = m.shape_dna(normalization=SPECTRUM_NORMALIZATION, count=k, threads=1,
                                           report=True, **SHAPE_DNA_OPTIONS)
                    values = {"lam": [float(x) for x in lam], "area": float(rep["area"])}
                elif task == "features":
                    values = {key: float(v) for key, v in m.boundary_features().items()}
                else:
                    method, c = task.split("#")
                    c = int(c)
                    kwargs = {}
                    if method in SEED_OPTIONS:
                        key = SEED_OPTIONS[method]
                        kwargs[key] = crossgen.options(method)[key] + c
                    b = getattr(copy_of(c), method)(**kwargs)
                    valid = bool(b.valid)
                    values = {"valid": valid, "coverage": float(b.coverage), "num_blocks": int(b.num_blocks)}
                    for name in METRICS:
                        v = float(getattr(b, name)())
                        values[name] = v if valid and math.isfinite(v) else None
                    # The decomposition holds its whole pipeline; drop it now,
                    # not when the next method's result rebinds the name.
                    del b
                emit(task=task, status="ok", seconds=time.time() - t0, values=values)
            except Exception as e:  # a method that raises is a method that failed
                emit(task=task, status="error", seconds=time.time() - t0, error=f"{type(e).__name__}: {e}")


# ── the parent ──────────────────────────────────────────────────────────────
def mesh_list(sources, limit):
    meshes = []
    for src in sources:
        d = SOURCES[src]
        if not d.is_dir():
            print(f"warning: {d} does not exist, skipping '{src}'"
                  + (" (make it with ./make_meshes.sh --mambo)" if src == "mambo" else ""), file=sys.stderr)
            continue
        for p in sorted(d.glob("*.obj")):
            # Faces of one mambo part are one group (<part>_<face>.obj), so a
            # grouped split keeps a part's faces on one side; a singlemat model
            # is its own group.
            group = f"{src}/{p.stem.rsplit('_', 1)[0]}" if src == "mambo" else f"{src}/{p.stem}"
            meshes.append({"mesh_id": f"{src}/{p.stem}", "group": group, "source": src, "path": str(p)})
    return meshes[:limit] if limit else meshes


def task_list(k, copies):
    return (["spectrum"] if k > 0 else []) + ["features"] + [f"{m}#{c}" for m in METHODS for c in range(copies)]


def missing_tasks(prior, k, copies):
    """The tasks of task_list(k, copies) that `prior`, a cache record of the
    mesh (or None), has no result for. A spectrum shorter than k counts as
    missing; one longer is cut to k when the CSVs are written. A task that
    failed has a result: a failure is final, as in the run that recorded it."""
    have = prior["results"] if prior else {}
    return [t for t in task_list(k, copies)
            if t not in have or (t == "spectrum" and prior["k"] < k)]


def run_mesh(mesh, k, copies, out_dir, progress, prior=None):
    """The tasks of one mesh that `prior` (an earlier, incomplete cache record
    of it, or None) lacks, in as many children as crashes and hangs need.
    Returns the mesh's cache record, prior's other results included. Nothing of
    the mesh is loaded here: the child loads it, and the OS takes all of it
    back when the child exits."""
    tasks = missing_tasks(prior, k, copies)
    results = {t: r for t, r in (prior["results"] if prior else {}).items() if t not in tasks}
    tmp = out_dir / "tmp"
    log = out_dir / "logs" / (mesh["mesh_id"].replace("/", "__") + ".log")
    env = dict(os.environ, OMP_NUM_THREADS="1")
    attempt = 0
    while tasks:
        attempt += 1
        rpath = tmp / f"{mesh['mesh_id'].replace('/', '__')}.{attempt}.jsonl"
        rpath.unlink(missing_ok=True)
        with open(log, "a") as logf:
            proc = subprocess.Popen(
                [sys.executable, __file__, "--worker", mesh["path"], "--eigs", str(k),
                 "--result", str(rpath), "--tasks", *tasks],
                stdout=logf, stderr=subprocess.STDOUT, env=env)
            # Poll for progress: a task that goes quiet for its timeout is hung.
            last, seen, running, counted = time.time(), 0, None, 0
            while proc.poll() is None:
                time.sleep(1.0)
                text = rpath.read_text() if rpath.exists() else ""
                lines = text[:text.rfind("\n") + 1].splitlines()   # complete lines only
                if len(lines) > seen:
                    seen, last = len(lines), time.time()
                    recs = [json.loads(l) for l in lines]
                    running = recs[-1]
                    finished = sum(r["status"] != "start" for r in recs)
                    progress.advance(finished - counted)
                    counted = finished
                limit = SPECTRUM_TIMEOUT if running and running["task"] == "spectrum" else METHOD_TIMEOUT
                if time.time() - last > limit:
                    proc.kill()
                    proc.wait()
                    break
        recs = [json.loads(l) for l in rpath.read_text().splitlines()] if rpath.exists() else []
        rpath.unlink(missing_ok=True)
        finished = [r for r in recs if r["status"] != "start"]
        progress.advance(len(finished) - counted)   # any that landed after the last poll
        for r in finished:
            results[r["task"]] = r
        # The task that started and never finished took the child down with it:
        # record it as failed and go on with the rest in a new child.
        started = [r["task"] for r in recs if r["status"] == "start"]
        stuck = next((t for t in started if t not in results), None)
        if stuck is None and proc.returncode != 0 and not started:
            stuck = tasks[0]            # died before it could start anything (load failed)
        if stuck is not None:
            why = "timeout" if proc.returncode in (-9, None) else f"crashed (exit {proc.returncode})"
            results[stuck] = {"task": stuck, "status": "error", "error": why}
            progress.advance(1)
        tasks = [t for t in tasks if t not in results]
    return {**mesh, "version": CACHE_VERSION, "k": k, "copies": copies, "results": results}


def mean(xs):
    return sum(xs) / len(xs) if xs else math.nan


def row_of(rec, k):
    """One dataset.csv row from a cache record: each method's copies averaged,
    a copy that errored counting as not valid, with coverage 0; PROBE's copy 0
    on its own."""
    row = {"mesh_id": rec["mesh_id"], "group": rec["group"], "source": rec["source"]}
    spec = rec["results"].get("spectrum", {})
    sv = spec["values"] if spec.get("status") == "ok" else {}
    feats = rec["results"].get("features", {})
    fv = feats.get("values", {}) if feats.get("status") == "ok" else {}
    # The two agree to 1e-13; the features are there whether or not the spectrum is.
    row["area"] = fv.get("area", sv.get("area", math.nan))
    lam = (list(sv.get("lam", [])) + [math.nan] * k)[:k]
    row.update({f"lam_{i + 1}": x for i, x in enumerate(lam)})
    row.update({f"feat_{name}": fv.get(name, math.nan) for name in FEATURES})
    if PROBE:
        r = rec["results"].get(f"{PROBE}#0", {})
        v = r["values"] if r.get("status") == "ok" else None
        valid = bool(v and v["valid"])
        row[f"probe_{PROBE}_valid"] = float(valid)
        for name in METRICS:
            row[f"probe_{PROBE}_{name}"] = v[name] if valid and v[name] is not None else math.nan
        row[f"probe_{PROBE}_coverage"] = v["coverage"] if v else 0.0
        row[f"probe_{PROBE}_num_blocks"] = float(v["num_blocks"]) if v else math.nan
    copies = rec["copies"]
    for method in METHODS:
        runs = [rec["results"].get(f"{method}#{c}", {}) for c in range(copies)]
        ok = [r["values"] for r in runs if r.get("status") == "ok"]
        valid = [v for v in ok if v["valid"]]
        row[f"{method}_valid"] = len(valid) / copies
        for name in METRICS:
            row[f"{method}_{name}"] = mean([v[name] for v in valid if v[name] is not None])
        row[f"{method}_coverage"] = mean([v["coverage"] for v in ok] + [0.0] * (copies - len(ok)))
        row[f"{method}_num_blocks"] = mean([v["num_blocks"] for v in valid])
        row[f"{method}_seconds"] = mean([r["seconds"] for r in runs if "seconds" in r])
    return row


def cache_records(cache_path, copies):
    """The cache's records made with `copies` copies, one at a time, so that
    nothing holds all of it. Whether one has every task a run needs is
    missing_tasks()'s question."""
    if not cache_path.exists():
        return
    with open(cache_path) as f:
        for line in f:
            rec = json.loads(line)
            if rec.get("version") == CACHE_VERSION and rec.get("copies") == copies:
                yield rec


METHOD_COLUMNS = ["valid", *METRICS, "coverage", "num_blocks", "seconds"]


def write_csvs(cache_path, wanted, k, copies, out_dir):
    """The three CSVs in one pass over the cache, a row at a time, from the
    first record of each mesh that has every task. Rows come in the order the
    meshes finished; sort on mesh_id downstream if it matters."""
    lam_cols = [f"lam_{i + 1}" for i in range(k)]
    probe_cols = [f"probe_{PROBE}_{c}" for c in PROBE_COLUMNS] if PROBE else []
    header = (["mesh_id", "group", "source", "area"] + lam_cols + [f"feat_{n}" for n in FEATURES]
              + probe_cols + [f"{m}_{c}" for m in METHODS for c in METHOD_COLUMNS])

    def fmt(x):
        if isinstance(x, float):
            return "" if math.isnan(x) else f"{x:.10g}"
        return x

    seen, n_rows, n_spec = set(), 0, 0
    with open(out_dir / "dataset.csv", "w", newline="") as fd, \
         open(out_dir / "spectra.csv", "w", newline="") as fs, \
         open(out_dir / "results.csv", "w", newline="") as fr:
        wd, ws, wr = csv.writer(fd), csv.writer(fs), csv.writer(fr)
        wd.writerow(header)
        ws.writerow(["mesh_id", "group"] + lam_cols)
        wr.writerow(["mesh_id", "algorithm", "metric", "value"])
        for rec in cache_records(cache_path, copies):
            if rec["mesh_id"] not in wanted or rec["mesh_id"] in seen or missing_tasks(rec, k, copies):
                continue
            seen.add(rec["mesh_id"])
            row = row_of(rec, k)
            wd.writerow([fmt(row[c]) for c in header])
            n_rows += 1
            # memo_data.py's pair: a mesh with no spectrum has nothing to learn from.
            if any(math.isnan(row[c]) for c in lam_cols):
                continue
            n_spec += 1
            ws.writerow([fmt(row[c]) for c in ["mesh_id", "group"] + lam_cols])
            for m in METHODS:
                for n in ["valid", *METRICS]:
                    wr.writerow([row["mesh_id"], m, n, fmt(row[f"{m}_{n}"])])
    print(f"{n_rows} meshes" + (f" ({n_rows - n_spec} without a spectrum)" if k else "")
          + f", {n_spec * len(METHODS) * (1 + len(METRICS))} results -> {out_dir}/dataset.csv, spectra.csv, results.csv")


class Progress:
    """One bar over the whole dataset, in tasks (the spectrum if asked for, the
    features and every copy of every method per mesh; cached tasks count as
    done), with a percentage, the mesh count and an ETA from this session's
    rate. On a terminal it redraws in place; into a file it writes a line at
    each whole percent instead."""
    WIDTH = 40

    def __init__(self, meshes_total, meshes_done, tasks_total, tasks_done):
        self.total, self.done = tasks_total, tasks_done
        self.meshes_total, self.meshes_done = meshes_total, meshes_done
        self.start_done, self.t0 = self.done, time.time()
        self.lock = threading.Lock()
        self.tty = sys.stderr.isatty()
        self.last_pct = -1

    def advance(self, n=1):
        if n > 0:
            with self.lock:
                self.done += n

    def mesh_finished(self):
        with self.lock:
            self.meshes_done += 1

    @staticmethod
    def hms(sec):
        sec = int(sec)
        return f"{sec // 3600}:{sec // 60 % 60:02d}:{sec % 60:02d}"

    def line(self):
        with self.lock:
            done, total, md = self.done, self.total, self.meshes_done
        frac = done / total if total else 1.0
        filled = int(self.WIDTH * frac)
        elapsed = time.time() - self.t0
        rate = (done - self.start_done) / elapsed if elapsed > 0 else 0.0
        eta = self.hms((total - done) / rate) if rate > 0 else "--:--:--"
        bar = "█" * filled + "░" * (self.WIDTH - filled)
        return frac, (f"[{bar}] {100 * frac:5.1f}%  {md}/{self.meshes_total} meshes  "
                      f"elapsed {self.hms(elapsed)}  ETA {eta}")

    def draw(self):
        frac, text = self.line()
        if self.tty:
            sys.stderr.write("\r" + text + "\x1b[K")
            sys.stderr.flush()
        elif int(100 * frac) != self.last_pct:
            self.last_pct = int(100 * frac)
            print(text, file=sys.stderr, flush=True)

    def note(self, msg):
        """A message above the bar, which is redrawn under it."""
        if self.tty:
            sys.stderr.write("\r\x1b[K")
        print(msg, file=sys.stderr, flush=True)
        self.draw()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default="dataset")
    ap.add_argument("--eigs", type=int, default=K_EIGS,
                    help=f"eigenvalues of the Shape-DNA spectrum per mesh (default {K_EIGS}; 0: no spectrum)")
    ap.add_argument("--copies", type=int, default=COPIES,
                    help=f"runs of each method per mesh, on symmetric copies (1..{len(SYMMETRIES)})")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--sources", nargs="+", default=list(SOURCES), choices=list(SOURCES))
    ap.add_argument("--limit", type=int, default=0, help="only the first N meshes (for a trial run)")
    ap.add_argument("--worker", help=argparse.SUPPRESS)
    ap.add_argument("--result", help=argparse.SUPPRESS)
    ap.add_argument("--tasks", nargs="+", help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.worker:
        worker(args.worker, args.tasks, args.eigs, args.result)
        return
    if not 1 <= args.copies <= len(SYMMETRIES):
        sys.exit(f"--copies must be in 1..{len(SYMMETRIES)}")
    if args.eigs < 0:
        sys.exit("--eigs must be 0 (no spectrum) or more")

    out_dir = Path(args.out)
    for d in (out_dir, out_dir / "logs", out_dir / "tmp"):
        d.mkdir(parents=True, exist_ok=True)
    cache_path = out_dir / "cache.jsonl"
    meshes = mesh_list(args.sources, args.limit)
    wanted = {m["mesh_id"] for m in meshes}
    # A mesh is done when a record of it has every task. Of the others, keep the
    # record with the fewest tasks missing (a mesh cached before --eigs asked
    # for a spectrum lacks only that), so that only the missing ones run. The
    # records of done meshes are not held.
    done, partial = set(), {}
    for rec in cache_records(cache_path, args.copies):
        mid = rec["mesh_id"]
        if mid not in wanted or mid in done:
            continue
        missing = len(missing_tasks(rec, args.eigs, args.copies))
        if not missing:
            done.add(mid)
            partial.pop(mid, None)
        elif mid not in partial or missing < len(missing_tasks(partial[mid], args.eigs, args.copies)):
            partial[mid] = rec
    todo = [m for m in meshes if m["mesh_id"] not in done]
    per_mesh = len(task_list(args.eigs, args.copies))
    left = sum(len(missing_tasks(partial.get(m["mesh_id"]), args.eigs, args.copies)) for m in todo)
    print(f"{len(meshes)} meshes, {len(meshes) - len(todo)} cached, {len(todo)} to run"
          + (f" ({len(partial)} of them in part)" if partial else "") + f" on {args.jobs} job(s), "
          f"{args.copies} cop{'y' if args.copies == 1 else 'ies'} of each method, "
          + (f"{args.eigs} eigenvalues" if args.eigs else "no spectrum"))
    progress = Progress(len(meshes), len(meshes) - len(todo), per_mesh * len(meshes), per_mesh * len(meshes) - left)
    progress.draw()
    # A finished future is dropped as soon as its record is on disk, so the
    # parent holds only the meshes in flight, never the ones already written.
    with open(cache_path, "a") as cache, cf.ThreadPoolExecutor(args.jobs) as pool:
        pending = {pool.submit(run_mesh, m, args.eigs, args.copies, out_dir, progress, partial.pop(m["mesh_id"], None))
                   for m in todo}
        while pending:
            finished, pending = cf.wait(pending, timeout=1.0)
            for fut in finished:
                rec = fut.result()
                cache.write(json.dumps(rec) + "\n")
                cache.flush()
                progress.mesh_finished()
                bad = [f"{t} ({r.get('error', '?')})" for t, r in rec["results"].items() if r["status"] != "ok"]
                if bad:
                    progress.note(f"{rec['mesh_id']}: failed {'; '.join(bad)}")
                del rec
            progress.draw()
    if progress.tty:
        sys.stderr.write("\n")

    write_csvs(cache_path, wanted, args.eigs, args.copies, out_dir)
    # Scratch only: every result and its error text are in cache.jsonl.
    for d in ("logs", "tmp"):
        shutil.rmtree(out_dir / d, ignore_errors=True)


if __name__ == "__main__":
    main()
