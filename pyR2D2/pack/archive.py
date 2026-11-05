"""HPCI (Gfarm) に移すために、ランの残りをまとめる (``r2d2plus-archive``)。

3D・slice・tau は先に ``r2d2plus-pack`` で 1 ステップ 1 個の ``.zarr.zip`` にしておく。
このコマンドは元のランには触らず、出力先に次のものを作る:

* ``data/restart/<slot>.tar``: restart の slot ごと (中のファイルの名前・大きさ・種類が元と一致することを確かめる)
* ``data/time.zip``: ``data/time`` (``time/mhd``・``time/vl``・``time/tau``・``time/slice`` など下を全部。
  中の名前は ``time/<族>/t.dac.N``。全メンバーの CRC と、名前・大きさが元と一致することを確かめる)
* ``rest.tar``: 残り全部のうち、ソースとバイナリを除いたもの (データはグループで読めてよい)
* ``rest_private.tar``: ソースとバイナリ (``data/param/source``・``bin``、``conversion/fortran_source``
  など。:func:`_is_private`)。ランの記録なので捨てないが、**権限 600 で作る**
  (堀田、2026-11-05「課題の中でもソースは読ませません。データは読んでもいいです」)。
  Gfarm へ送った後も ``gfls -l`` で 600 を確かめ、違えば ``gfchmod 600`` する
* ``archive_report.toml``: 上の照合の結果 (restart は slot ごとの中身の一覧と tar の SHA-256。
  ``r2d2plus-restore`` が戻すときに使う) と、``rest.tar`` に入れずにランから直接送る zip の一覧

rest.tar に入れないのは、restart・time と、3D・slice・tau の元のファイルと zip
(``remap/qq``、``slice/qq*``・``slice/zarr``、``tau/qq*``・``tau/zarr``)。zip はステップごとに
取り出せるよう、そのままの形で (ランから) 送る。元のファイルのうち zip の無いステップが残っていれば、
それが rest.tar から漏れないよう、止める (先に r2d2plus-pack を流す)。
"""

import datetime
import hashlib
import os
import socket
import subprocess
import tarfile
import time
import zipfile
from pathlib import Path

import pyR2D2
from pyR2D2 import _term as term

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


# 中身を読まなくてもソースと分かる名前 (大文字小文字を区別する: .F と .f は別の拡張子)
_PRIVATE_SUFFIXES = frozenset(
    ".c .cc .cp .cpp .cxx .c++ .C .h .hh .hpp .hxx .h++ .H .inc .ipp .tpp .inl .cu .cuh "
    ".f .F .for .FOR .ftn .fpp .f77 .F77 .f90 .F90 .f95 .F95 .f03 .F03 .f08 .F08 .mod .o .a .so".split()
)
_PRIVATE_NAMES = frozenset({"Makefile", "makefile", "GNUmakefile", "CMakeLists.txt", "CMakePresets.json"})


def _is_private(rel, path=None):
    """ソース・バイナリか (``rest_private.tar`` に分ける)。

    ``rel`` はランからの相対パス (POSIX)。下のどれかに当たれば真。迷う名前は秘密の側に倒す:
    分けすぎてもランの記録は ``rest_private.tar`` に残るが、漏れると読めてしまうため。

    * ディレクトリの名前が ``bin``、または ``source`` を含む (``fortran_source`` など)
    * ファイルの名前に ``source`` を含む (``*.r2d2plus_source.tar.zst`` など)
    * ソースの拡張子、``Makefile`` など
    * ``path`` を渡せば、ELF (実行ファイル) かも中身の先頭で確かめる
    """
    parts = rel.split("/")
    if any(p == "bin" or "source" in p.lower() for p in parts):
        return True
    name = parts[-1]
    if name in _PRIVATE_NAMES or os.path.splitext(name)[1] in _PRIVATE_SUFFIXES:
        return True
    if path is not None and not path.is_symlink() and path.is_file():
        with open(path, "rb") as f:
            return f.read(4) == b"\x7fELF"
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


def sha256(path, bufsize=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(bufsize):
            h.update(chunk)
    return h.hexdigest()


def _make_tar(tar_path, root, rels, private=False):
    """``root`` の下の ``rels`` だけを (下を辿らずに) tar にし、中身の名前・種類・大きさを元と比べる。

    ``private`` なら tar を権限 600 で作る (作っている途中も他人に読めない)。
    """
    listing = tar_path.with_name(tar_path.name + ".list")
    listing.write_bytes(b"".join(r.encode() + b"\0" for r in rels))
    work = tar_path.with_name(tar_path.name + ".packing")
    try:
        subprocess.run(
            ["tar", "-cf", str(work), "-C", str(root), "--no-recursion", "--null", "-T", str(listing)],
            check=True,
            umask=0o077 if private else -1,
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
    return sorted(expected)


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
        self._log(term.dim(f"[{time.time() - self._t0:8.1f}s]"), *msg)

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
            members = _make_tar(tar_path, src, rels)
            self._record(
                part="restart", slot=slot.name, tar=tar_path.relative_to(self.out).as_posix(),
                entries=len(members), bytes=tar_path.stat().st_size, sha256=sha256(tar_path),
                seconds=round(time.time() - t0, 1),
                check="tar members (name, type, size) == source",
                # r2d2plus-restore が展開した後に照合する (名前, 種類 f/d/l, 大きさ)
                members=[list(m) for m in members],
            )
            self.log(f"restart {slot.name}: {len(members)} entries, {tar_path.stat().st_size / 1e9:.3f} GB")

    def time_zip(self):
        src = self.datadir / "time"
        if not src.is_dir():
            return
        t0 = time.time()
        dest = self.out / "data" / "time.zip"
        dest.parent.mkdir(parents=True, exist_ok=True)
        work = dest.with_name("time.packing.zip")
        files = {r: (src / r).stat().st_size for r in _walk(src) if (src / r).is_file()}
        # Data.zip_time・Data.time_read と同じく time/<族>/t.dac.N の名前で入れる (2026-10-31 までは <族>/… だった)
        with zipfile.ZipFile(work, "w", zipfile.ZIP_STORED, allowZip64=True) as zf:
            for rel in files:
                zf.write(src / rel, f"time/{rel}")
        with zipfile.ZipFile(work) as zf:
            bad = zf.testzip()
            members = {i.filename.removeprefix("time/"): i.file_size for i in zf.infolist()}
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
        private = [r for r in rels if _is_private(r, self.run / r)]
        hidden = set(private)
        rels = [r for r in rels if r not in hidden]
        tar_path = self.out / "rest.tar"
        self.out.mkdir(parents=True, exist_ok=True)
        n = len(_make_tar(tar_path, self.run, rels))
        # 前の版の rest.tar を残さず置き換えたので、ここでは秘密の側だけ作り直す
        self._rest_private(private)
        self._record(
            part="rest", tar="rest.tar", entries=n, bytes=tar_path.stat().st_size,
            seconds=round(time.time() - t0, 1),
            check="tar members (name, type, size) == source",
            excluded="data/restart, data/time, data/remap/qq, data/slice/qq*, data/slice/zarr, data/tau/qq*, data/tau/zarr",
        )
        self.log(f"rest.tar: {n} entries, {tar_path.stat().st_size / 1e9:.3f} GB")

    def _rest_private(self, rels):
        """ソースとバイナリの tar (600)。無ければ作らず、前の回のものは消す。"""
        tar_path = self.out / "rest_private.tar"
        if not rels:
            tar_path.unlink(missing_ok=True)
            return
        t0 = time.time()
        n = len(_make_tar(tar_path, self.run, rels, private=True))
        mode = tar_path.stat().st_mode & 0o777
        if mode != 0o600:
            raise RuntimeError(f"{tar_path}: mode {mode:o}, must be 600")
        self._record(
            part="rest_private", tar="rest_private.tar", entries=n, bytes=tar_path.stat().st_size,
            mode="600", seconds=round(time.time() - t0, 1),
            check="tar members (name, type, size) == source; mode 600",
        )
        self.log(
            f"rest_private.tar (source and binaries, mode 600): {n} entries; "
            f"{term.yellow('keep it 600 on Gfarm too')} (gfls -l; gfchmod 600)"
        )

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
        self.log(f"{term.green('done')}; {len(zips)} zips are to be sent from the run itself (see archive_report.toml)")
