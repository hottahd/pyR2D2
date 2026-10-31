"""r2d2plus-pack・r2d2plus-archive・r2d2plus-restore。

    r2d2plus-pack <ラン> [--steps 81,100-200] [--kinds qf,qs,qt] [--delete-original] [--out DIR] [--dry-run]
    r2d2plus-archive <ラン> <出力先> [--parts restart,time,rest] [--dry-run]
    r2d2plus-restore <ラン> --slot o [--from gfarm:///…/run/<ラン>] [--dry-run]

詳しくは pyR2D2.pack.pack、pyR2D2.pack.archive、pyR2D2.pack.restore と doc/source/data_compression.rst。
"""

import argparse
import subprocess
import sys

from .archive import Archiver
from .pack import KINDS, Packer, parse_steps
from .restore import GFARM_ROOT, Restorer


def build_parser():
    ap = argparse.ArgumentParser(
        prog="r2d2plus-pack",
        description=(
            "Pack 3D (qf), slice (qs) and tau (qt) outputs of a run into one .zarr.zip per step, "
            "checked bit for bit against the binary files. Results go to data/param/pack_report.toml."
        ),
    )
    ap.add_argument("run", help="run directory or its data/")
    ap.add_argument("--steps", help="steps such as 81,100-200 (default: all)")
    ap.add_argument("--kinds", default=",".join(KINDS), help="comma separated subset of qf,qs,qt")
    ap.add_argument(
        "--delete-original",
        action="store_true",
        help="delete the binary files of steps checked in this call (the last step is kept for a "
        "wall_limit run). In a running run, only steps followed by a later time file "
        "(data/time/mhd or tau) are packed and deleted; restart and data/time are never touched",
    )
    ap.add_argument(
        "--out",
        help="write the zips and the report under OUT/data/ instead of the run (nothing is written "
        "in the run; cannot be used with --delete-original)",
    )
    ap.add_argument("--workers", type=int, default=8, help="threads reading 3D binary files")
    ap.add_argument("--dry-run", action="store_true", help="only show what would be done")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    kinds = [k for k in args.kinds.split(",") if k]
    unknown = [k for k in kinds if k not in KINDS]
    if unknown:
        raise SystemExit(f"unknown kinds {unknown}; choose from {KINDS}")
    packer = Packer(args.run, workers=args.workers, log=lambda *m: print(*m, flush=True), out=args.out)
    ok = packer.run(
        kinds=kinds,
        steps=parse_steps(args.steps),
        delete_original=args.delete_original,
        dry_run=args.dry_run,
    )
    return 0 if ok else 1


def build_archive_parser():
    ap = argparse.ArgumentParser(
        prog="r2d2plus-archive",
        description=(
            "Put the rest of a run into the form moved to HPCI: data/restart/<slot>.tar, data/time.zip "
            "and rest.tar under OUT, each checked against the run. The run is not changed. "
            "Pack 3D, slice and tau with r2d2plus-pack first; their zips are sent from the run."
        ),
    )
    ap.add_argument("run", help="run directory")
    ap.add_argument("out", help="output directory (outside the run)")
    ap.add_argument("--parts", default="restart,time,rest", help="comma separated subset of restart,time,rest")
    ap.add_argument("--dry-run", action="store_true", help="only check and show what would be done")
    return ap


def archive_main(argv=None):
    args = build_archive_parser().parse_args(argv)
    parts = [p for p in args.parts.split(",") if p]
    unknown = [p for p in parts if p not in ("restart", "time", "rest")]
    if unknown:
        raise SystemExit(f"unknown parts {unknown}")
    archiver = Archiver(args.run, args.out, log=lambda *m: print(*m, flush=True))
    try:
        archiver.run_all(parts=parts, dry_run=args.dry_run)
    except RuntimeError as exc:
        print(f"r2d2plus-archive: {exc}", file=sys.stderr)
        return 1
    return 0


def build_restore_parser():
    ap = argparse.ArgumentParser(
        prog="r2d2plus-restore",
        description=(
            "Bring a restart slot of a run moved to HPCI back: fetch data/restart/<slot>.tar and "
            "archive_report.toml (made by r2d2plus-archive), check the tar's SHA-256, extract it to "
            "data/restart/<slot>/ and check names, types and sizes against the report. "
            "Run it on csgw2 (gfpcopy and /vol0600 are both there)."
        ),
    )
    ap.add_argument("run", help="run directory on Fugaku")
    ap.add_argument("--slot", required=True, help="restart slot such as o, e or 00000080")
    ap.add_argument(
        "--from", dest="source",
        help="directory written by r2d2plus-archive, gfarm: URL or local path "
        f"(default: {GFARM_ROOT}/<project>/run/<run name>, as Sync's work/<project>/run/<caseid>)",
    )
    ap.add_argument(
        "--project",
        help="project for the default --from (default: taken from the run path …/<project>/run/<run name>)",
    )
    ap.add_argument("--dry-run", action="store_true", help="only read the report and show what would be done")
    return ap


def restore_main(argv=None):
    args = build_restore_parser().parse_args(argv)
    try:
        restorer = Restorer(
            args.run, source=args.source, project=args.project, log=lambda *m: print(*m, flush=True)
        )
        restorer.run_slot(args.slot, dry_run=args.dry_run)
    except (RuntimeError, ValueError, FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"r2d2plus-restore: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
