"""HPCI (Gfarm) に移したランの restart を、富岳のランへ戻す (``r2d2plus-restore``)。

``r2d2plus-archive`` が作った ``data/restart/<slot>.tar`` と ``archive_report.toml`` を
取ってきて (Gfarm なら ``gfpcopy``。csgw2 で動かす)、``data/restart/<slot>/`` に展開する。

照合:

1. tar の SHA-256 が ``archive_report.toml`` と同じ
2. 展開したものの名前・種類・大きさが ``archive_report.toml`` の一覧 (``members``) と同じ

作業は ``data/restart/.restore-<slot>/`` で行い、照合が通ってから ``data/restart/<slot>`` にする
(R2D2plus が途中のものを読むことはない)。``data/restart/<slot>`` がすでにあれば何もしない。
結果は ``data/param/restore_report.toml`` に足していく。
"""

import datetime
import os
import shutil
import socket
import subprocess
import tarfile
import time
from pathlib import Path

try:
    import tomllib
except ImportError:  # pragma: no cover
    tomllib = None

import pyR2D2

from .archive import _entries, _walk, sha256
from .pack import _toml_value

# r2d2plus-archive の出力を置く Gfarm の場所 (R2D2plus-project docs/run_data.md)。
# Sync の work/<project>/run/<caseid> と同じく、<project>/run/<caseid> を後ろに付ける
GFARM_ROOT = "gfarm:///home/hp200137/hpci000196/fugaku"


def default_source(run, project=None):
    """既定の取ってくる場所 ``GFARM_ROOT/<project>/run/<caseid>``。

    ``project`` が無ければ、ランのパス ``…/<project>/run/<caseid>`` から取る
    (富岳の ``/vol0600/…/u10094/R2D2plus/run/d00012`` なら ``R2D2plus``)。
    """
    run = Path(run)
    if project is None:
        if run.parent.name != "run":
            raise ValueError(
                f"cannot tell the project from {run} (not …/<project>/run/<caseid>); give --project or --from"
            )
        project = run.parent.parent.name
    return f"{GFARM_ROOT}/{project}/run/{run.name}"


class Restorer:
    """
    restart の slot を 1 つ戻す

    Parameters
    ----------
    run : str
        戻す先のランのディレクトリ (``data/`` の親)
    source : str
        ``r2d2plus-archive`` の出力先 (``archive_report.toml`` のあるディレクトリ)。
        ``gfarm:`` で始まれば ``gfpcopy`` で取ってくる。既定は :py:func:`default_source`
    project : str
        既定の場所に使うプロジェクトの名前 (既定はランのパスから取る)
    """

    def __init__(self, run, source=None, log=print, gfpcopy="gfpcopy", project=None):
        run = Path(run).resolve()
        if run.name == "data" and (run / "param").is_dir():
            run = run.parent
        self.run = run
        self.restart = run / "data" / "restart"
        self.source = source or default_source(run, project)
        self.gfpcopy = gfpcopy
        self._log = log
        self._t0 = time.time()

    def log(self, *msg):
        self._log(f"[{time.time() - self._t0:8.1f}s]", *msg)

    def _is_gfarm(self):
        return self.source.startswith("gfarm:")

    def _fetch(self, rel, dest_dir):
        """``source`` の ``rel`` を ``dest_dir`` へ持ってきて、そのパスを返す。"""
        dest_dir.mkdir(parents=True, exist_ok=True)
        name = rel.rsplit("/", 1)[-1]
        if self._is_gfarm():
            src = self.source.rstrip("/") + "/" + rel
            subprocess.run([self.gfpcopy, "-f", src, str(dest_dir)], check=True)
        else:
            src = Path(self.source) / rel
            if not src.is_file():
                raise FileNotFoundError(src)
            shutil.copyfile(src, dest_dir / name)
        path = dest_dir / name
        if not path.is_file():
            raise FileNotFoundError(f"{rel} did not arrive in {dest_dir}")
        return path

    def report(self, work):
        path = self._fetch("archive_report.toml", work)
        with open(path, "rb") as f:
            return tomllib.load(f)

    @staticmethod
    def slots(report):
        return {p["slot"]: p for p in report.get("part", []) if p.get("part") == "restart"}

    def run_slot(self, slot, dry_run=False):
        dest = self.restart / slot
        if dest.exists():
            raise RuntimeError(f"{dest} already exists; not restored")
        work = self.restart / f".restore-{slot}"
        if work.exists():
            shutil.rmtree(work)
        try:
            report = self.report(work)
            entry = self.slots(report).get(slot)
            if entry is None:
                raise RuntimeError(
                    f"no restart slot {slot!r} in {self.source}/archive_report.toml "
                    f"(slots: {sorted(self.slots(report))})"
                )
            if dry_run:
                self.log(
                    f"would restore {slot} ({entry['entries']} entries, {entry['bytes'] / 1e9:.3f} GB) "
                    f"from {self.source} to {dest}"
                )
                return None
            t0 = time.time()
            tar_path = self._fetch(entry["tar"], work)
            self.log(f"fetched {entry['tar']} ({tar_path.stat().st_size / 1e9:.3f} GB)")
            if tar_path.stat().st_size != entry["bytes"]:
                raise RuntimeError(f"{tar_path}: {tar_path.stat().st_size} bytes, the report says {entry['bytes']}")
            if "sha256" in entry and sha256(tar_path) != entry["sha256"]:
                raise RuntimeError(f"{tar_path}: SHA-256 differs from archive_report.toml")
            out = work / "x"
            out.mkdir()
            with tarfile.open(tar_path) as tf:
                tf.extractall(out, filter="data")
            got = _entries(out, [slot] + [f"{slot}/{r}" for r in _walk(out / slot)])
            if "members" in entry:
                expected = {(m[0], m[1], int(m[2])) for m in entry["members"]}
                check = "tar SHA-256 and extracted (name, type, size) == archive_report.toml"
            else:  # members を残す前の archive (2026-10-31 以前)
                expected = None
                check = "extracted entry count == archive_report.toml (no member list in the report)"
                if len(got) != entry["entries"]:
                    raise RuntimeError(f"{len(got)} entries extracted, the report says {entry['entries']}")
            if expected is not None and got != expected:
                missing = sorted(expected - got)[:5]
                extra = sorted(got - expected)[:5]
                raise RuntimeError(f"extracted {slot} differs from the report (missing {missing}, extra {extra})")
            os.rename(out / slot, dest)
            self._record(
                slot=slot, source=f"{self.source}/{entry['tar']}", entries=len(got),
                bytes=entry["bytes"], check=check, seconds=round(time.time() - t0, 1),
            )
            self.log(f"restored {dest} ({len(got)} entries)")
            return dest
        finally:
            if work.exists():
                shutil.rmtree(work)

    def _record(self, **entry):
        path = self.run / "data" / "param" / "restore_report.toml"
        new = not path.exists()
        entry["date"] = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
        entry["host"] = socket.gethostname()
        entry["pyR2D2"] = getattr(pyR2D2, "__version__", "")
        with open(path, "a") as f:
            if new:
                f.write("# r2d2plus-restore (pyR2D2.pack.restore) の記録。戻した restart の slot ごとに [[restore]] を 1 つ足す。\n")
            f.write("\n[[restore]]\n")
            for key, value in entry.items():
                f.write(f"{key} = {_toml_value(value)}\n")
