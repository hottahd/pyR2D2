"""リモートで zarr (ディレクトリまたは .zarr.zip) の一部を読み、標準出力に流すスクリプト。

:py:meth:`pyR2D2.Sync.remap_zarr` などが ``ssh server python -`` の標準入力に
このファイルの中身と ``main('<要求>')`` の 1 行を送って実行する。リモートに要るのは numpy と zarr (v3) だけで、
pyR2D2 は要らない。リモートには何も書かない (一時ファイルを作らない)。

要求 (``main`` に渡す JSON)::

    {"files": [{"path": ..., "keys": [...] or null, "region": [[a, b] or null, ...] or null}],
     "block": 軸 0 方向に一度に送るチャンクの数}

出力: ``MAGIC`` のあと、長さ 8 バイト (little endian) + JSON の記録と、配列のバイト列の並び。

* ``{"file": path, "missing": true}``: ファイルがない
* ``{"file": path, "attrs": ..., "arrays": {name: {"shape", "chunks", "dtype", "whole"}}}``
* ``{"name": name, "sel": [[a, b], ...], "dtype", "shape"}`` のあとにバイト列: 配列 ``name`` の ``sel`` の部分
* ``{"end_file": path}`` / ``{"end": true}``
"""

import json
import struct
import sys
from pathlib import Path

MAGIC = b"PYR2D2-ZARR-REGION-1\n"


def _open(path):
    import zarr

    path = Path(path).expanduser()
    if path.is_dir():
        return zarr.open_group(path, mode="r")
    for p in [path, Path(str(path) + ".zip")]:
        if p.is_file() and p.suffix == ".zip":
            return zarr.open_group(store=zarr.storage.ZipStore(p, mode="r"), mode="r")
    return None


def _record(out, obj):
    data = json.dumps(obj).encode()
    out.write(struct.pack("<Q", len(data)))
    out.write(data)


def _send_file(out, req, block):
    import numpy as np

    root = _open(req["path"])
    if root is None:
        _record(out, {"file": req["path"], "missing": True})
        return

    names = list(root.array_keys())
    ndim = max(root[name].ndim for name in names)
    keys = req.get("keys")
    if keys is not None:
        # slice の ro_yin などのように、keys の ro で ro_<添字> も選ぶ
        def wanted(name, key):
            return name == key or name.startswith(key + "_")

        unknown = [k for k in keys if not any(wanted(name, k) for name in names)]
        if unknown:
            raise KeyError(f"{req['path']}: no array {unknown} (has {names})")

    arrays = {}
    for name in names:
        arr = root[name]
        # 次元の低い配列 (3D の x, y, z など) は座標とみなし、keys によらず丸ごと送る
        whole = arr.ndim < ndim
        if not whole and keys is not None and not any(wanted(name, k) for k in keys):
            continue
        arrays[name] = {
            "shape": list(arr.shape),
            "chunks": list(arr.chunks),
            "dtype": str(arr.dtype),
            "whole": whole,
        }
    _record(out, {"file": req["path"], "attrs": root.attrs.asdict(), "arrays": arrays})

    region = req.get("region") or [None] * ndim
    for name, meta in arrays.items():
        arr = root[name]
        shape = meta["shape"]
        if meta["whole"]:
            sel = [[0, s] for s in shape]
        else:
            sel = [[0, s] if r is None else [max(0, r[0]), min(s, r[1])] for r, s in zip(region, shape)]
        if any(a >= b for a, b in sel):
            raise ValueError(f"{req['path']}:{name}: empty region {sel} for shape {shape}")
        # 軸 0 をチャンクの境目で区切って送り、リモートのメモリを抑える
        step = meta["chunks"][0] * block
        a = sel[0][0]
        while a < sel[0][1]:
            b = min(sel[0][1], (a // step + 1) * step)
            part = [[a, b]] + sel[1:]
            data = np.ascontiguousarray(arr[tuple(slice(p, q) for p, q in part)])
            _record(out, {"name": name, "sel": part, "dtype": data.dtype.str, "shape": list(data.shape)})
            # パイプには numpy の tofile が使えないので、バイト列をそのまま書く
            out.write(memoryview(data).cast("B"))
            a = b
    _record(out, {"end_file": req["path"]})


def main(request):
    request = json.loads(request)
    out = sys.stdout.buffer
    out.write(MAGIC)
    for req in request["files"]:
        _send_file(out, req, request.get("block", 1))
        out.flush()
    _record(out, {"end": True})
    out.flush()


# pyR2D2 は標準入力でこのファイルを送り、その後ろに main('<要求>') の 1 行を足す
if __name__ == "__main__" and len(sys.argv) > 1:
    main(sys.argv[1])
