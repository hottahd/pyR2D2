"""リモートの zarr (``.zarr.zip`` を含む) から一部だけを ssh 越しに持ってくる。

リモートでは ``_remote_extract.py`` を ``python -`` で動かし、要る範囲だけを読んで
標準出力に流す (リモートには何も書かない)。手元では、元と同じ形・同じチャンクの
zarr ディレクトリを作り、持ってきた範囲だけを書く。書いていない所は NaN になるので、
:py:class:`pyR2D2.Data` の ``qx``, ``qr`` などで元のデータと同じ添字のまま読める。
持ってきた範囲は属性 ``pyR2D2_partial`` に残す。
"""

import json
import shlex
import struct
import subprocess
from pathlib import Path

import numpy as np

from ..data_io import zarr_util
from ._remote_extract import MAGIC

_SCRIPT = Path(__file__).with_name("_remote_extract.py")

PARTIAL_ATTR = "pyR2D2_partial"


def _covers(entry, keys, region):
    """手元に記録済みの ``entry`` が、要求 (``keys``, ``region``) を含むか。"""
    if entry["keys"] is not None and (keys is None or not set(keys) <= set(entry["keys"])):
        return False
    if entry["region"] is None:
        return True
    if region is None or len(region) != len(entry["region"]):
        return False
    for have, want in zip(entry["region"], region):
        if have is None:
            continue
        if want is None or not (have[0] <= want[0] and want[1] <= have[1]):
            return False
    return True


def local_status(path, keys, region):
    """手元の ``path`` の状態: ``"full"`` (部分でない zarr がある)、``"done"`` (要求を含む)、``None``。"""
    path = Path(path)
    if Path(str(path) + ".zip").exists():
        return "full"
    if not path.is_dir():
        return None
    root = zarr_util.open_zarr_group(path)
    attrs = root.attrs.asdict()
    if PARTIAL_ATTR not in attrs:
        return "full"
    if any(_covers(e, keys, region) for e in attrs[PARTIAL_ATTR]):
        return "done"
    return None


def _read_exact(stream, size):
    data = stream.read(size)
    if len(data) != size:
        raise EOFError("stream ended early")
    return data


def _read_record(stream):
    (size,) = struct.unpack("<Q", _read_exact(stream, 8))
    return json.loads(_read_exact(stream, size))


def _skip_to_magic(stream, limit=1 << 20):
    """シェルの初期化ファイルが標準出力に何か書いても、MAGIC の後ろから読む。"""
    buf = b""
    while not buf.endswith(MAGIC):
        c = stream.read(1)
        if not c or len(buf) > limit:
            raise EOFError("no response from the remote script")
        buf += c


def _open_local(path, header, source):
    """手元の zarr を開き、足りない配列を元と同じ形・チャンク・NaN 埋めで作る。"""
    zarr = zarr_util._require_zarr()
    path.parent.mkdir(parents=True, exist_ok=True)
    root = zarr.open_group(path, mode="a")
    attrs = root.attrs.asdict()
    for key, value in header["attrs"].items():
        if key not in attrs:
            root.attrs[key] = value
    if PARTIAL_ATTR not in attrs:
        root.attrs[PARTIAL_ATTR] = []
        root.attrs["pyR2D2_source"] = source
    names = set(root.array_keys())
    for name, meta in header["arrays"].items():
        if name in names:
            if list(root[name].shape) != meta["shape"]:
                raise ValueError(
                    f"{path}:{name} has shape {root[name].shape}, remote has {meta['shape']}"
                )
            continue
        root.create_array(
            name,
            shape=meta["shape"],
            chunks=meta["chunks"],
            dtype=meta["dtype"],
            fill_value=np.nan,
            compressors=zarr_util._codec_f32(),
        )
    return root


def fetch(server, files, ssh="ssh", remote_python="python3", block=1, verbose=True):
    """
    リモートの zarr の一部を、1 回の ssh で手元の zarr ディレクトリに書く。

    Parameters
    ----------
    server : str
        ssh の接続先
    files : list of dict
        ``{"remote": リモートのパス, "local": 手元のパス, "keys": 変数名の list か None,
        "region": 各軸の [始め, 終わり) か None の list、または None}``
    ssh : str
        ssh のコマンド (``"ssh -p 22"`` のような引数付きでもよい)
    remote_python : str
        リモートの Python (numpy と zarr 3 が入っているもの)
    block : int
        軸 0 方向に一度に送るチャンクの数。リモートのメモリはこれに比例する

    Returns
    -------
    list of Path
        書いた手元のパス (リモートにないものは除く)
    """
    if not files:
        return []
    request = {
        "files": [{"path": str(f["remote"]), "keys": f["keys"], "region": f["region"]} for f in files],
        "block": block,
    }
    by_remote = {str(f["remote"]): f for f in files}
    script = _SCRIPT.read_text() + f"\nmain({json.dumps(request)!r})\n"
    command = shlex.split(ssh) + [server, f"{remote_python} -"]
    proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    written = []
    try:
        proc.stdin.write(script.encode())
        proc.stdin.close()
        _skip_to_magic(proc.stdout)
        root = current = None
        while True:
            rec = _read_record(proc.stdout)
            if rec.get("end"):
                break
            if "missing" in rec:
                print(f"not found on {server}: {rec['file']}")
            elif "attrs" in rec:
                current = by_remote[rec["file"]]
                source = f"{server}:{rec['file']}"
                root = _open_local(Path(current["local"]), rec, source)
                if verbose:
                    print(f"{source} -> {current['local']}")
            elif "name" in rec:
                dtype = np.dtype(rec["dtype"])
                size = dtype.itemsize * int(np.prod(rec["shape"]))
                data = np.frombuffer(_read_exact(proc.stdout, size), dtype=dtype)
                data = data.reshape(rec["shape"])
                root[rec["name"]][tuple(slice(a, b) for a, b in rec["sel"])] = data
            elif "end_file" in rec:
                # 全部書けてから記録する。途中で切れたら記録が残らず、次にまた持ってくる
                entries = root.attrs.asdict()[PARTIAL_ATTR]
                entries.append({"keys": current["keys"], "region": current["region"]})
                root.attrs[PARTIAL_ATTR] = entries
                written.append(Path(current["local"]))
                root = current = None
    except EOFError as exc:
        proc.wait()
        raise RuntimeError(
            f"remote extraction on {server} failed (exit {proc.returncode}); "
            "see the messages above"
        ) from exc
    finally:
        proc.stdout.close()
    if proc.wait() != 0:
        raise RuntimeError(f"remote extraction on {server} exited with {proc.returncode}")
    return written
