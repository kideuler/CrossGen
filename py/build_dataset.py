#!/usr/bin/env python3
"""Build the algorithm-selection dataset: a Shape-DNA spectrum per mesh, and
three block-decomposition scores for each of the five methods.

    PYTHONPATH=build/python python3 py/build_dataset.py [--out dataset] [--jobs 4]
        [--eigs 100] [--sources singlemat mambo] [--limit N]

Writes, in --out:
  dataset.csv   one row per mesh: mesh_id, group, source, area, lam_1..lam_K,
                then <method>_<metric> for 5 methods x 3 metrics (15 columns).
  spectra.csv   mesh_id, group, lam_1..lam_K         } the two files memo_data.py
  results.csv   mesh_id, algorithm, metric, value    } writes, for train_selector.py
  cache.jsonl   every task's raw result, one line per mesh; a rerun skips the
                meshes already in it, so an interrupted run resumes.
  logs/, tmp/   each mesh's C++ output and its workers' result files while it
                runs; both are deleted when the script finishes.

Every number is written exactly as the C++ side returns it; nothing is
rescaled here.

The spectrum is ShapeDNA with normalization="weyl_ratio": lambda_k |Omega| /
(4 pi k), the area-normalised eigenvalue over Weyl's leading term. Left
area-normalised only, lam_100 is ~100x lam_1 on every shape and the columns
mostly encode k; divided by Weyl's line they sit near 1 and differ between
shapes by what the leading term does not know -- the boundary term (~ perimeter
/ sqrt(area k), the Dirichlet correction) and below. SPECTRUM_NORMALIZATION
takes any other ShapeDNA normalisation ("area", "none", ...) instead.

The metrics are BlockDecomposition's regularity(), angleQuality() and
chordQuality(): in [0, 1], higher better, and 0 for a decomposition that does
not cover the model. A method that raises, crashes the process or makes no
progress for METHOD_TIMEOUT scores 0 on all three too -- no layout is the worst
layout.

Each mesh runs in its own child process (this script with --worker), so a
segfault or a hang in one method costs that method and not the run, and the
parent retries the methods after it in a fresh child. Children are single-
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

K_EIGS = 100                       # the length of the spectrum
SPECTRUM_NORMALIZATION = "weyl_ratio"   # a ShapeDNA normalisation; see crossgen.Mesh.shape_dna
SHAPE_DNA_OPTIONS = dict(degree=3, boundary="dirichlet")
SOURCES = {"singlemat": REPO / "data/meshes/singlemat", "mambo": REPO / "data/meshes/mambo"}
METHODS = ["zipline", "umber", "meridian", "torsion", "atlas"]
METRICS = ["regularity", "angle_quality", "chord_quality"]   # BlockDecomposition methods = CSV names
CACHE_VERSION = 2                  # records of any other version are ignored (and rerun)
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
    """Run `tasks` ("spectrum" and method names) on one mesh, appending one JSON
    line per task to `result_path` -- a "start" line before it, so the parent
    knows which task a crash or hang belongs to, and the result after."""
    crossgen = import_crossgen()
    m = crossgen.load(mesh_path)
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
                else:
                    b = getattr(m, task)()
                    values = {name: float(getattr(b, name)()) for name in METRICS}
                    values["coverage"] = float(b.coverage)
                    values["num_blocks"] = int(b.num_blocks)
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


def run_mesh(mesh, k, out_dir, progress):
    """All tasks for one mesh, in as many children as crashes and hangs need.
    Returns the mesh's cache record. Nothing of the mesh is loaded here: the
    child loads it, and the OS takes all of it back when the child exits."""
    tasks = ["spectrum"] + METHODS
    results = {}
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
    return {**mesh, "version": CACHE_VERSION, "k": k, "results": results}


def row_of(rec, k):
    """One dataset.csv row from a cache record."""
    row = {"mesh_id": rec["mesh_id"], "group": rec["group"], "source": rec["source"]}
    spec = rec["results"].get("spectrum", {})
    if spec.get("status") == "ok":
        v = spec["values"]
        row["area"] = v["area"]
        lam = v["lam"]
    else:
        row["area"], lam = math.nan, []
    lam = (lam + [math.nan] * k)[:k]
    row.update({f"lam_{i + 1}": x for i, x in enumerate(lam)})
    for method in METHODS:
        r = rec["results"].get(method, {})
        for name in METRICS:
            row[f"{method}_{name}"] = r["values"][name] if r.get("status") == "ok" else 0.0
    return row


def cache_records(cache_path, k):
    """The cache, one record at a time, so that nothing holds all of it."""
    if not cache_path.exists():
        return
    with open(cache_path) as f:
        for line in f:
            rec = json.loads(line)
            if rec.get("version") == CACHE_VERSION and rec["k"] == k:
                yield rec


def write_csvs(cache_path, wanted, k, out_dir):
    """The three CSVs in one pass over the cache, a row at a time. Rows come in
    the order the meshes finished; sort on mesh_id downstream if it matters."""
    lam_cols = [f"lam_{i + 1}" for i in range(k)]
    header = ["mesh_id", "group", "source", "area"] + lam_cols + [f"{m}_{n}" for m in METHODS for n in METRICS]

    def fmt(x):
        return f"{x:.10g}" if isinstance(x, float) else x

    seen, n_rows, n_spec = set(), 0, 0
    with open(out_dir / "dataset.csv", "w", newline="") as fd, \
         open(out_dir / "spectra.csv", "w", newline="") as fs, \
         open(out_dir / "results.csv", "w", newline="") as fr:
        wd, ws, wr = csv.writer(fd), csv.writer(fs), csv.writer(fr)
        wd.writerow(header)
        ws.writerow(["mesh_id", "group"] + lam_cols)
        wr.writerow(["mesh_id", "algorithm", "metric", "value"])
        for rec in cache_records(cache_path, k):
            if rec["mesh_id"] not in wanted or rec["mesh_id"] in seen:
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
                for n in METRICS:
                    wr.writerow([row["mesh_id"], m, n, fmt(row[f"{m}_{n}"])])
    print(f"{n_rows} meshes ({n_rows - n_spec} without a spectrum), {n_spec * len(METHODS) * len(METRICS)} "
          f"results -> {out_dir}/dataset.csv, spectra.csv, results.csv")


class Progress:
    """One bar over the whole dataset, in tasks (a spectrum and five methods per
    mesh; cached meshes count as done), with a percentage, the mesh count and an
    ETA from this session's rate. On a terminal it redraws in place; into a
    file it writes a line at each whole percent instead."""
    WIDTH = 40

    def __init__(self, meshes_total, meshes_done):
        per = 1 + len(METHODS)
        self.total, self.done = per * meshes_total, per * meshes_done
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
    ap.add_argument("--eigs", type=int, default=K_EIGS)
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

    out_dir = Path(args.out)
    for d in (out_dir, out_dir / "logs", out_dir / "tmp"):
        d.mkdir(parents=True, exist_ok=True)
    cache_path = out_dir / "cache.jsonl"
    done = {rec["mesh_id"] for rec in cache_records(cache_path, args.eigs)}

    meshes = mesh_list(args.sources, args.limit)
    wanted = {m["mesh_id"] for m in meshes}
    todo = [m for m in meshes if m["mesh_id"] not in done]
    print(f"{len(meshes)} meshes, {len(meshes) - len(todo)} cached, {len(todo)} to run on {args.jobs} job(s)")
    progress = Progress(len(meshes), len(meshes) - len(todo))
    progress.draw()
    # A finished future is dropped as soon as its record is on disk, so the
    # parent holds only the meshes in flight, never the ones already written.
    with open(cache_path, "a") as cache, cf.ThreadPoolExecutor(args.jobs) as pool:
        pending = {pool.submit(run_mesh, m, args.eigs, out_dir, progress) for m in todo}
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

    write_csvs(cache_path, wanted, args.eigs, out_dir)
    # Scratch only: every result and its error text are in cache.jsonl.
    for d in ("logs", "tmp"):
        shutil.rmtree(out_dir / d, ignore_errors=True)


if __name__ == "__main__":
    main()
