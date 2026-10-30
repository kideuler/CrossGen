"""Mesh singlemat, MAMBO and DeepCAD at one NP into a directory of OBJs, then tar it up.

    python3 py/make_mesh_archive.py 100 /Volumes/ext/meshes_np100     # -> meshes_np100/ + meshes_np100.tar.gz
    python3 py/make_mesh_archive.py 64 ~/meshes_np64 --zstd --jobs 8
    python3 py/make_mesh_archive.py 100 /tmp/try --limit 200          # a trial: 200 DeepCAD domains only

The output directory gets one subdirectory per corpus, laid out like data/meshes:

    <out>/singlemat/*.obj          Mesh2Dgmsh on data/geometry/singlemat/*.geo
    <out>/mambo/*.obj, shapes.txt  Mesh3Dto2Dfacesgmsh on every .step/.stp under ~/Documents/mambo
    <out>/deepcad/*.obj            Mesh2Dgmsh on data/geometry/deepcad/*.geo
    <out>/deepcad/manifest.csv     deepcad_to_geo.py's manifest, the rows of the domains meshed here
    <out>/meshing_log.txt          NP, commit, counts, what was skipped and what failed

and the archive is <out>.tar.gz (or .tar.zst with --zstd) beside it, holding the directory under its own name.

### How this differs from make_meshes.sh

The meshers and their arguments are make_meshes.sh's (h = 2/NP; MAMBO's faces through one --registry, so each
shape is written once), but

* nothing is cleaned or rebuilt. make_meshes.sh sources clean.sh and build.sh, which delete data/meshes/**.obj
  and `rm -rf build`; here build/Mesh2Dgmsh and build/Mesh3Dto2Dfacesgmsh must already exist, and data/meshes is
  never touched, so an archive at another NP leaves the working corpus as it is.
* the .geo files are meshed in parallel. DeepCAD's ~98k domains take ~0.35 s each, 9-10 h one after another.
  MAMBO's STEP files stay serial, in make_meshes.sh's order (`find | sort`, which is Python's sorted() on these
  names): the registry is appended to as faces are written and decides which copy of a repeated shape is kept,
  so the chain runs as one task beside the .geo pool and its file names match data/meshes/mambo at the same NP.
* every call has a time limit (--timeout), so one stuck gmsh costs one domain, not the run.

### Which DeepCAD domains are meshed

Mesh3Dto2Dfacesgmsh drops a MAMBO face that is not nice at the NP it is given -- a corner-to-corner stretch of
boundary under 3 mesh edges, or thinner than 4 mesh edges read as 2A/P -- because a layout cannot resolve a
feature the mesh does not. data/geometry/deepcad was filtered by the same rule, but once, at deepcad_to_geo.py's
--np, which was 100 (GEO_NP). Its manifest records each domain's shortest stretch and 2A/P in mesh edges at that
NP (run_h, width_h), and both scale with the mesh size, so at NP = N a domain is meshed only if
run_h N/GEO_NP >= 3 and width_h N/GEO_NP >= 4: the MAMBO faces and the DeepCAD domains then obey one rule at
every N. Nothing is dropped at N >= 100; at N = 50 it is 17 % of them. --all-deepcad meshes every file anyway.
Above 100 the files are still only those nice at 100; a profile that becomes nice only at a finer mesh needs
`deepcad_to_geo.py --np N` first.

### Size

Mesh2Dgmsh writes 17 fixed digits per coordinate, ~0.5 MB per domain at NP = 100, growing as NP^2 (the triangle
count): DeepCAD alone is ~48 GB at NP = 100, and gzip takes it to about a third of that. The free space where the
output goes is checked against that estimate before anything runs (--force skips the check). zstd made the same
size as gzip -6 ten times faster on these files, but gzip is on every machine that might unpack it.

### The archive

tar runs with COPYFILE_DISABLE=1. Without it macOS's tar stores a file's extended attributes as a `._<name>`
entry beside it, which tar on Linux unpacks as a real file -- `._B0_0.obj`, which any `*.obj` glob then reads as
a mesh. It is written as <archive>.part and renamed only once tar and the compressor have both succeeded.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
GEOMETRY = REPO / "data" / "geometry"
MESH2D = REPO / "build" / "Mesh2Dgmsh"
MESH3D = REPO / "build" / "Mesh3Dto2Dfacesgmsh"
DEFAULT_MAMBO = Path.home() / "Documents" / "mambo"

GEO_NP = 100                  # deepcad_to_geo.py's --np when data/geometry/deepcad was written
MIN_RUN, MIN_WIDTH = 3.0, 4.0  # Mesh3Dto2Dfacesgmsh's --min-run / --min-width defaults
MB_AT_100 = 0.5               # an OBJ at NP = 100, measured on 400 DeepCAD, the MAMBO and the singlemat meshes
GZIP_RATIO = 2.8              # tar.gz size / OBJ size on the same 400


def run(cmd, timeout):
    """One mesher call: (ok, the last lines it printed, for the log when it failed)."""
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, "timed out after %d s" % timeout
    tail = p.stdout.decode(errors="replace").strip().splitlines()[-5:]
    return p.returncode == 0, ("exit %d: " % p.returncode if p.returncode else "") + " | ".join(tail)


def mesh_geo(geo, out_dir, n, timeout, stop):
    """Mesh2Dgmsh on one .geo. A zero exit is not enough: the .obj has to be there and not empty."""
    if stop.is_set():
        return geo, None, "cancelled"
    ok, msg = run([str(MESH2D), str(geo), str(out_dir), str(n)], timeout)
    obj = out_dir / (geo.stem + ".obj")
    if ok and obj.is_file() and obj.stat().st_size > 0:
        return geo, True, ""
    obj.unlink(missing_ok=True)   # half-written by a call that timed out or crashed
    return geo, False, msg or "no .obj written"


def mesh_mambo(steps, out_dir, n, timeout, stop):
    """Every STEP file in order, one registry (see the docstring). Returns the failures as (file, message)."""
    registry = out_dir / "shapes.txt"
    failed = []
    for step in steps:
        if stop.is_set():
            failed.append((step, "cancelled"))
            continue
        # The output name is a printf format, so a literal '%' in the stem is doubled, as make_meshes.sh does.
        fmt = out_dir / (step.stem.replace("%", "%%") + "_%d.obj")
        ok, msg = run([str(MESH3D), str(step), str(fmt), "0", str(n), "--registry", str(registry)], timeout)
        if not ok:
            failed.append((step, msg))
    return failed


def deepcad_domains(n, all_deepcad):
    """The DeepCAD .geo files to mesh at NP = n, their manifest rows, and how many the N rule leaves out."""
    geo_dir = GEOMETRY / "deepcad"
    manifest = geo_dir / "manifest.csv"
    geos = sorted(geo_dir.glob("*.geo"))
    if not manifest.is_file():
        if not all_deepcad and n < GEO_NP:
            sys.exit("%s is missing, so the domains not nice at NP %d cannot be told apart; "
                     "rerun py/deepcad_to_geo.py or pass --all-deepcad" % (manifest, n))
        return geos, {}, 0
    with open(manifest, newline="") as fp:   # csv, not split(","): a few sketch names hold commas
        reader = csv.DictReader(fp)
        rows = {r["name"]: r for r in reader}
        cols = reader.fieldnames
    scale = n / GEO_NP
    keep = [g for g in geos if all_deepcad or g.stem not in rows
            or (float(rows[g.stem]["run_h"]) * scale >= MIN_RUN
                and float(rows[g.stem]["width_h"]) * scale >= MIN_WIDTH)]
    return keep, (rows, cols), len(geos) - len(keep)


def free_bytes(path):
    """Free space on the disk that will hold path, which need not exist yet."""
    p = path
    while not p.exists():
        p = p.parent
    return shutil.disk_usage(p).free


def git_state():
    """The commit the meshers were built from, as far as git can say (read-only calls)."""
    try:
        head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, timeout=30).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, timeout=30).stdout.strip()
        return head + (" + local changes" if dirty else "") if head else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def make_archive(out, archive, zstd):
    """tar the directory under its own name through gzip or zstd into archive (see the docstring)."""
    part = archive.with_name(archive.name + ".part")
    env = dict(os.environ, COPYFILE_DISABLE="1")
    compress = ["zstd", "-q", "-T0", "-3"] if zstd else ["gzip", "-6"]
    with open(part, "wb") as fp:
        tar = subprocess.Popen(["tar", "-cf", "-", "-C", str(out.parent), out.name], stdout=subprocess.PIPE, env=env)
        comp = subprocess.Popen(compress, stdin=tar.stdout, stdout=fp)
        tar.stdout.close()   # so that tar sees a broken pipe if the compressor dies
        comp_rc, tar_rc = comp.wait(), tar.wait()
    if tar_rc or comp_rc:
        part.unlink(missing_ok=True)
        sys.exit("archiving failed (tar exit %d, %s exit %d); the meshes are still in %s"
                 % (tar_rc, compress[0], comp_rc, out))
    part.rename(archive)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("np", type=int, help="make_meshes.sh's NP: the mesh size is h = 2/NP")
    ap.add_argument("out", type=Path, help="new or empty directory for the meshes; the archive goes beside it")
    ap.add_argument("--jobs", type=int, default=os.cpu_count(), help="meshers running at once")
    ap.add_argument("--mambo-dir", type=Path, default=DEFAULT_MAMBO, help="searched recursively for .step/.stp")
    ap.add_argument("--all-deepcad", action="store_true", help="mesh every DeepCAD .geo, nice at NP or not")
    ap.add_argument("--limit", type=int, default=0, help="only the first K DeepCAD domains (by name), for a trial")
    ap.add_argument("--timeout", type=int, default=900, help="seconds one mesher call may take")
    ap.add_argument("--zstd", action="store_true", help=".tar.zst instead of .tar.gz (much faster, same size)")
    ap.add_argument("--force", action="store_true", help="skip the free-space check")
    args = ap.parse_args()

    n = args.np
    if n <= 0:
        sys.exit("NP must be a positive integer")
    out = args.out.expanduser().resolve()
    archive = out.with_name(out.name + (".tar.zst" if args.zstd else ".tar.gz"))
    for exe in (MESH2D, MESH3D):
        if not os.access(exe, os.X_OK):
            sys.exit("%s is not built: cmake --build build --target Mesh2Dgmsh Mesh3Dto2Dfacesgmsh -j8" % exe)
    if args.zstd and shutil.which("zstd") is None:
        sys.exit("--zstd needs the zstd program (brew install zstd)")
    # The MAMBO registry has to start empty, and anything already in the directory would go into the archive.
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        sys.exit("%s exists and is not an empty directory; pick a new one" % out)
    if archive.exists():
        sys.exit("%s exists; move it away first" % archive)

    singlemat = sorted((GEOMETRY / "singlemat").glob("*.geo"))
    deepcad, manifest, dropped = deepcad_domains(n, args.all_deepcad)
    if args.limit:
        deepcad = deepcad[:args.limit]
    steps = []
    if args.mambo_dir.is_dir():
        steps = sorted(p for p in args.mambo_dir.rglob("*") if p.is_file() and p.suffix.lower() in (".step", ".stp"))
    if not steps:
        sys.exit("no .step/.stp under %s (--mambo-dir)" % args.mambo_dir)

    # The MAMBO count is data/meshes/mambo's at NP 100, near enough for an estimate at any NP.
    need = MB_AT_100 * 1e6 * (n / 100.0) ** 2 * (len(singlemat) + len(deepcad) + 408) * (1 + 1 / GZIP_RATIO)
    free = free_bytes(out)
    print("NP %d -> %s: %d singlemat, %d DeepCAD (%d not nice at NP %d left out), %d MAMBO STEP files"
          % (n, out, len(singlemat), len(deepcad), dropped, n, len(steps)))
    print("about %.1f GB of meshes + archive, %.1f GB free" % (need / 1e9, free / 1e9), flush=True)
    if need > free and not args.force:
        sys.exit("not enough space for the estimate; pick a directory on a larger disk, lower NP, "
                 "use --limit, or --force")

    dirs = {k: out / k for k in ("singlemat", "mambo", "deepcad")}
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    # A full run takes hours, and a Mac that idles to sleep suspends every mesher with it: a trial at NP 50 sat
    # 18 min on one STEP file that takes 3 s, the 18 min the lid was shut. caffeinate -w holds the system awake
    # until this process exits, however it exits. A closed lid on battery still sleeps.
    if sys.platform == "darwin" and shutil.which("caffeinate"):
        subprocess.Popen(["caffeinate", "-i", "-w", str(os.getpid())])

    stop = threading.Event()
    start = time.time()
    geo_jobs = [(g, dirs["singlemat"]) for g in singlemat] + [(g, dirs["deepcad"]) for g in deepcad]
    geo_failed, done, last = [], 0, start
    pool = ThreadPoolExecutor(max_workers=max(1, args.jobs))
    try:
        # MAMBO first, so its serial chain overlaps the pool rather than trailing it.
        mambo = pool.submit(mesh_mambo, steps, dirs["mambo"], n, args.timeout, stop)
        futures = [pool.submit(mesh_geo, g, d, n, args.timeout, stop) for g, d in geo_jobs]
        for fut in as_completed(futures):
            geo, ok, msg = fut.result()
            done += 1
            if not ok:
                geo_failed.append((geo, msg))
            now = time.time()
            if now - last >= 60 or done == len(futures):
                last = now
                eta = (now - start) / done * (len(futures) - done)
                print("  .geo %d/%d, %d failed, %.0f s, ~%.0f min left%s"
                      % (done, len(futures), len(geo_failed), now - start, eta / 60,
                         "" if mambo.done() else "; MAMBO still running"), flush=True)
        if not mambo.done():
            print("  waiting for MAMBO's serial chain", flush=True)
        mambo_failed = mambo.result()
    except KeyboardInterrupt:
        stop.set()
        pool.shutdown(wait=False, cancel_futures=True)
        sys.exit("interrupted; %s is incomplete and was not archived" % out)
    pool.shutdown()
    elapsed = time.time() - start

    if manifest:
        rows, cols = manifest
        with open(dirs["deepcad"] / "manifest.csv", "w", newline="") as fp:
            w = csv.DictWriter(fp, fieldnames=cols)
            w.writeheader()
            for g in deepcad:
                if g.stem in rows and (dirs["deepcad"] / (g.stem + ".obj")).is_file():
                    w.writerow(rows[g.stem])

    counts = {k: len(list(d.glob("*.obj"))) for k, d in dirs.items()}
    size = sum(f.stat().st_size for d in dirs.values() for f in d.iterdir())
    lines = ["NP %d (h = 2/%d), %s, CrossGen %s" % (n, n, time.strftime("%Y-%m-%d %H:%M"), git_state()),
             "command: %s" % " ".join(sys.argv),
             "meshed in %.0f s with %d jobs, %.2f GB of OBJ" % (elapsed, args.jobs, size / 1e9),
             "singlemat: %d of %d .geo" % (counts["singlemat"], len(singlemat)),
             "mambo: %d faces from %d STEP files in %s (nice at NP %d, each shape once: shapes.txt)"
             % (counts["mambo"], len(steps), args.mambo_dir, n),
             "deepcad: %d of %d .geo%s%s" % (counts["deepcad"], len(deepcad),
                                            "; %d more not nice at NP %d left out" % (dropped, n) if dropped else "",
                                            "; --limit %d" % args.limit if args.limit else ""),
             "failed: %d .geo, %d STEP" % (len(geo_failed), len(mambo_failed))]
    lines += ["  %s: %s" % (p.relative_to(REPO) if p.is_relative_to(REPO) else p, m)
              for p, m in geo_failed + mambo_failed]
    (out / "meshing_log.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:7]), flush=True)
    if geo_failed or mambo_failed:
        print("  (the failures are listed in %s)" % (out / "meshing_log.txt"))

    print("archiving -> %s" % archive, flush=True)
    t0 = time.time()
    make_archive(out, archive, args.zstd)
    print("%s: %.2f GB in %.0f s" % (archive, archive.stat().st_size / 1e9, time.time() - t0))


if __name__ == "__main__":
    main()
