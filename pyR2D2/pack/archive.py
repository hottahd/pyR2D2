"""HPCI (Gfarm) に移すために、ランの残りをまとめる (``r2d2plus-archive``)。

3D・slice・tau は先に ``r2d2plus-pack`` で 1 ステップ 1 個の ``.zarr.zip`` にしておく。
このコマンドは元のランには触らず、出力先に次のものを作る:

* ``data/restart/<slot>.tar``: restart の slot ごと (中のファイルの名前・大きさ・種類が元と一致することを確かめる)
* ``data/time.zip``: ``data/time`` (全メンバーの CRC と、名前・大きさが元と一致することを確かめる)
* ``rest.tar``: 残り全部 (``data/param/source`` と ``bin`` も入れる。ランの記録なので)
* ``archive_report.toml``: 上の照合の結果と、``rest.tar`` に入れずにランから直接送る zip の一覧

rest.tar に入れないのは、restart・time と、3D・slice・tau の元のファイルと zip
(``remap/qq``、``slice/qq*``・``slice/zarr``、``tau/qq*``・``tau/zarr``)。zip はステップごとに
取り出せるよう、そのままの形で (ランから) 送る。元のファイルのうち zip の無いステップが残っていれば、
それが rest.tar から漏れないよう、止める (先に r2d2plus-pack を流す)。
"""

import datetime
import os
import socket
import subprocess
import tarfile
import time
import zipfile
from pathlib import Path

import pyR2D2

from .pack import KINDS, Packer, _toml_value, run_status


def _is_bulk(rel):
    """rest.tar に入れない ``data/`` の下のパス (``rel`` は ``data/`` からの相対、POSIX)。"""
    parts = rel.split("/")
    if parts[0] in ("restart", "time"):
        return True
    if parts[:2] == ["remap", "qq"]:
        return True
    if parts[0] in ("slice", "tau") and len(parts) >= 2:
        return parts[1] == "zarr" or parts[1].startswith("qq")
    return False


def _entries(root, rels):
    """``root`` からの相対パス ``rels`` の (名前, 種類, 大きさ) の集合。tar の中と比べる形。"""
    out = set()
    for rel in rels:
        p = root / rel
        if p.is_symlink():
            out.add((rel, "l", 0))
        elif p.is_dir():
            out.add((rel, "d", 0))
        else:
            out.add((rel, "f", p.stat().st_size))
    return out


def _tar_entries(path):
    out = set()
    with tarfile.open(path) as tf:
        for m in tf:
            name = m.name[2:] if m.name.startswith("./") else m.name
            name = name.rstrip("/")
            kind = "l" if m.issym() else ("d" if m.isdir() else ("f" if m.isfile() else "?"))
            out.add((name, kind, m.size if kind == "f" else 0))
    return out


def _walk(root, skip=None):
    """``root`` の下のファイル・リンク・ディレクトリの相対パス (``skip(rel)`` が真のものは下ごと除く)。"""
    rels = []
    for dirpath, dirnames, filenames in os.walk(root):
        base = Path(dirpath).relative_to(root)
        keep = []
        for d in sorted(dirnames):
            rel = (base / d).as_posix()
            if skip is not None and skip(rel):
                continue
            rels.append(rel)
            if not (Path(dirpath) / d).is_symlink():
                keep.append(d)
        dirnames[:] = keep
        for f in sorted(filenames):
            rel = (base / f).as_posix()
            if skip is None or not skip(rel):
                rels.append(rel)
    return rels


def _make_tar(tar_path, root, rels):
    """``root`` の下の ``rels`` だけを (下を辿らずに) tar にし、中身の名前・種類・大きさを元と比べる。"""
    listing = tar_path.with_name(tar_path.name + ".list")
    listing.write_bytes(b"".join(r.encode() + b"\0" for r in rels))
    work = tar_path.with_name(tar_path.name + ".packing")
    try:
        subprocess.run(
            ["tar", "-cf", str(work), "-C", str(root), "--no-recursion", "--null", "-T", str(listing)],
            check=True,
        )
    finally:
        listing.unlink()
    expected = _entries(root, rels)
    got = _tar_entries(work)
    if got != expected:
        missing = sorted(expected - got)[:5]
        extra = sorted(got - expected)[:5]
        raise RuntimeError(f"{tar_path}: tar differs from the source (missing {missing}, extra {extra})")
    os.replace(work, tar_path)
    return len(expected)


class Archiver:
    """
    1 つのランを、HPCI に移す形に出力先へまとめる

    Parameters
    ----------
    run : str
        ランのディレクトリ (``data/`` の親)
    out : str
        出力先 (ランの外)
    """

    def __init__(self, run, out, log=print):
        run = Path(run).resolve()
        if run.name == "data" and (run / "param").is_dir():
            run = run.parent
        self.run = run
        self.datadir = run / "data"
        self.out = Path(out).resolve()
        if self.out == run or run in self.out.parents or self.out in run.parents:
            raise ValueError("out must be outside the run")
        self._log = log
        self._t0 = time.time()
        self.report = []

    def log(self, *msg):
        self._log(f"[{time.time() - self._t0:8.1f}s]", *msg)

    def _record(self, **entry):
        entry["date"] = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
        self.report.append(entry)

    def unpacked(self):
        """元のファイルが残っていて zip の無いステップ (rest.tar から漏れるもの)。"""
        packer = Packer(self.datadir, log=lambda *m: None)
        out = []
        for unit, _ in packer.units(KINDS):
            if any(p.exists() for p in unit.raw) and not unit.zip.exists():
                out.append(unit.label())
        return out

    def zips(self):
        """ランから直接送る zip (``data/`` からの相対パス)。"""
        out = []
        for sub in ["remap/qq/zarr", "slice/zarr", "tau/zarr"]:
            d = self.datadir / sub
            if d.is_dir():
                out += sorted(p.relative_to(self.datadir).as_posix() for p in d.rglob("*.zarr.zip"))
        return out

    def restart(self):
        src = self.datadir / "restart"
        if not src.is_dir():
            return
        dest = self.out / "data" / "restart"
        dest.mkdir(parents=True, exist_ok=True)
        for slot in sorted(src.iterdir()):
            t0 = time.time()
            rels = [slot.name] + [f"{slot.name}/{r}" for r in _walk(slot)]
            tar_path = dest / f"{slot.name}.tar"
            n = _make_tar(tar_path, src, rels)
            self._record(
                part="restart", slot=slot.name, tar=tar_path.relative_to(self.out).as_posix(),
                entries=n, bytes=tar_path.stat().st_size, seconds=round(time.time() - t0, 1),
                check="tar members (name, type, size) == source",
            )
            self.log(f"restart {slot.name}: {n} entries, {tar_path.stat().st_size / 1e9:.3f} GB")

    def time_zip(self):
        src = self.datadir / "time"
        if not src.is_dir():
            return
        t0 = time.time()
        dest = self.out / "data" / "time.zip"
        dest.parent.mkdir(parents=True, exist_ok=True)
        work = dest.with_name("time.packing.zip")
        files = {r: (src / r).stat().st_size for r in _walk(src) if (src / r).is_file()}
        with zipfile.ZipFile(work, "w", zipfile.ZIP_STORED, allowZip64=True) as zf:
            for rel in files:
                zf.write(src / rel, rel)
        with zipfile.ZipFile(work) as zf:
            bad = zf.testzip()
            members = {i.filename: i.file_size for i in zf.infolist()}
        if bad is not None or members != files:
            raise RuntimeError(f"{work}: CRC error {bad} or members differ from data/time")
        os.replace(work, dest)
        self._record(
            part="time", zip=dest.relative_to(self.out).as_posix(), files=len(files),
            bytes=dest.stat().st_size, seconds=round(time.time() - t0, 1),
            check="testzip (CRC of every member); member names and sizes == data/time",
        )
        self.log(f"time.zip: {len(files)} files")

    def rest(self):
        t0 = time.time()

        def skip(rel):
            return rel.startswith("data/") and _is_bulk(rel[len("data/"):])

        rels = _walk(self.run, skip)
        tar_path = self.out / "rest.tar"
        self.out.mkdir(parents=True, exist_ok=True)
        n = _make_tar(tar_path, self.run, rels)
        self._record(
            part="rest", tar="rest.tar", entries=n, bytes=tar_path.stat().st_size,
            seconds=round(time.time() - t0, 1),
            check="tar members (name, type, size) == source",
            excluded="data/restart, data/time, data/remap/qq, data/slice/qq*, data/slice/zarr, data/tau/qq*, data/tau/zarr",
        )
        self.log(f"rest.tar: {n} entries, {tar_path.stat().st_size / 1e9:.3f} GB")

    def write_report(self, zips):
        path = self.out / "archive_report.toml"
        with open(path, "w") as f:
            f.write("# r2d2plus-archive (pyR2D2.pack.archive) の記録\n")
            f.write(f"run = {_toml_value(str(self.run))}\n")
            f.write(f"run_status = {_toml_value(run_status(self.datadir) or '')}\n")
            f.write(f"host = {_toml_value(socket.gethostname())}\n")
            f.write(f"pyR2D2 = {_toml_value(getattr(pyR2D2, '__version__', ''))}\n")
            f.write("# rest.tar に入れず、ランから直接送る zip (data/ からの相対パス)\n")
            f.write(f"zips = {_toml_value(zips)}\n")
            for entry in self.report:
                f.write("\n[[part]]\n")
                for key, value in entry.items():
                    f.write(f"{key} = {_toml_value(value)}\n")

    def run_all(self, parts=("restart", "time", "rest"), dry_run=False):
        status = run_status(self.datadir)
        if status == "running":
            raise RuntimeError(f"{self.run} is running (param/run_summary.toml); not archived")
        left = self.unpacked()
        if left:
            raise RuntimeError(
                f"binary outputs without a zip remain ({len(left)}: {left[:5]} ...); "
                "run r2d2plus-pack first, or they would be left out of rest.tar"
            )
        zips = self.zips()
        if dry_run:
            self.log(f"run status {status!r}; would write {parts} to {self.out}; {len(zips)} zips stay in the run")
            return
        self.out.mkdir(parents=True, exist_ok=True)
        if "restart" in parts:
            self.restart()
        if "time" in parts:
            self.time_zip()
        if "rest" in parts:
            self.rest()
        self.write_report(zips)
        self.log(f"done; {len(zips)} zips are to be sent from the run itself (see archive_report.toml)")
