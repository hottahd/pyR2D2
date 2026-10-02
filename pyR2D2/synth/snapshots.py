"""3D スナップショットが実在するかの確認（ODF-radiation odfgen.r2d2io から移した。R2D2plus DEC-589）。"""
import glob
import os


def check_snapshots(datadir, ns):
    """**3D 出力が実在するか確かめる** (無ければ例外)。

    **なぜ必要か**: R2D2 の `qr.read(n)` は**存在しない番号を渡しても例外を投げず、前回読んだ中身をそのまま残す**。
    間引いた run に対して素朴に `range(n0, n1)` で回すと、**欠けた番号が直前の番号の複製になり、時間変動を測る解析が
    黙って壊れる**。legacy (``qq.dac.*``) と R2D2plus の compressed (``qq.z.*``、2026-10-02 に対応) の両方を探す。

    Parameters
    ----------
    datadir : str
        `.../data/` (末尾の `/` は任意)。
    ns : int の列
        調べたいスナップショット番号。

    Returns
    -------
    ns : list[int]  (全部あればそのまま返す)
    """
    root = os.path.join(str(datadir).rstrip("/"), "remap", "qq")
    ns = [int(n) for n in ns]

    def exists(n):
        return any(glob.glob(os.path.join(root, "*", "*", f"qq.{kind}.{n:08d}.*")) for kind in ("dac", "z"))

    missing = [n for n in ns if not exists(n)]
    if missing:
        have = sorted({int(os.path.basename(f).split(".")[2])
                       for kind in ("dac", "z")
                       for f in glob.glob(os.path.join(root, "*", "00000000", f"qq.{kind}.*.00000000"))})
        raise FileNotFoundError(
            f"3D 出力が無いスナップショット {len(missing)} 個: {missing}\n"
            f"  間引き済みの run である可能性が高い。実在するのは {len(have)} 個で、例えば {have[:8]} ... {have[-4:]}\n"
            f"  **黙って直前の中身を読んでしまうので必ず止める**")
    return ns
