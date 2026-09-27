"""台帳の列の単位を値に合わせて自動で決める。

量は「種類」(時刻・長さ・R_star 比・角度・容量) と、その種類の基準単位での値
(時刻は秒、長さは cm、容量はバイト) で持つ。表示の単位は列ごとに 1 つ選ぶ。

単位の選び方
------------
列の典型値 (0 でない値の絶対値の中央値) を ``v`` とする。

1. 見出しに今の単位があり、同じ種類で ``0.1 <= v/単位 < 1e4`` なら、それを保つ
   (境界付近の値が 1 つ増えただけで列全体を書き直さないためのヒステリシス)。
2. そうでなければ、下の表で ``v < 上限`` となる最初の単位を選ぶ。どの単位でも
   表示値は ``[1, 1000)`` に入り、複数の単位が該当するときは慣用の繰り上がり
   (1000 s までは s、48 h までは h など) で 1 つに決める。

   ======== =============================================
   時刻     s (<1000 s), min (<100 min), h (<48 h), d (<730 d), yr
   長さ     km (<1000 km), Mm
   容量     MB (<1000 MB), GB (<1000 GB), TB  (10 進、1 GB = 1e9 B)
   R_star比 R_star のみ
   角度     deg のみ
   ======== =============================================

   例: 240 s → s、3 日 → d、-5.44e8 cm → Mm、2.4e6 cm → km。
"""

import math
import re
from dataclasses import dataclass

YEAR = 365.25 * 86400.0

# 種類ごとの (単位名, 基準単位での大きさ, この単位を使う典型値の上限)
UNITS = {
    "time": [
        ("s", 1.0, 1.0e3),
        ("min", 60.0, 6.0e3),
        ("h", 3600.0, 48 * 3600.0),
        ("d", 86400.0, 730 * 86400.0),
        ("yr", YEAR, math.inf),
    ],
    "length": [
        ("km", 1.0e5, 1.0e8),
        ("Mm", 1.0e8, math.inf),
    ],
    "ratio": [("R_star", 1.0, math.inf)],
    "angle": [("deg", 1.0, math.inf)],
    "bytes": [
        ("MB", 1.0e6, 1.0e9),
        ("GB", 1.0e9, 1.0e12),
        ("TB", 1.0e12, math.inf),
    ],
}

# 典型値が決まらない (値が全部 0) ときの既定
DEFAULT_UNIT = {
    "time": "s",
    "length": "Mm",
    "ratio": "R_star",
    "angle": "deg",
    "bytes": "GB",
}

# 単位名 → (種類, 大きさ)
UNIT_TABLE = {
    name: (family, factor)
    for family, units in UNITS.items()
    for name, factor, _ in units
}

_HEADER_RE = re.compile(r"^(.*?)\s*\[([^\[\]]*)\]\s*$")


@dataclass(frozen=True)
class Quantity:
    """種類と基準単位での値 (``dx`` のように 2 つ以上の数を持つ列もある)。"""

    family: str
    values: tuple

    def canonical(self):
        """ハッシュ用の表現 (浮動小数点は repr で丸めない)。"""
        return [self.family, [repr(float(v)) for v in self.values]]


def split_header(text):
    """見出し ``"dtout [s]"`` を ``("dtout", "s")`` に分ける。単位が無ければ ``None``。"""
    text = (text or "").strip()
    m = _HEADER_RE.match(text)
    if m:
        return m.group(1).strip(), m.group(2).strip()
    return text, None


def typical(values):
    """0 でない値の絶対値の中央値。全部 0 (または空) なら ``None``。"""
    xs = sorted(abs(v) for v in values if v != 0 and math.isfinite(v))
    if not xs:
        return None
    n = len(xs)
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def choose_unit(family, values, current=None):
    """種類 ``family`` の値の並び (基準単位) に合う表示単位を返す。

    ``current`` は今の見出しの単位。同じ種類で典型値が 0.1〜1e4 に収まるなら保つ。
    """
    units = UNITS[family]
    v = typical(values)
    cur = UNIT_TABLE.get(current) if current else None
    if cur is not None and cur[0] != family:
        cur = None
    if v is None:
        return current if cur is not None else DEFAULT_UNIT[family]
    if cur is not None and 0.1 <= v / cur[1] < 1.0e4:
        return current
    for name, _, upper in units:
        if v < upper:
            return name
    return units[-1][0]


def format_number(x):
    """表の数値の書式 (今のシートと同じ小数 2 桁。極端な値だけ有効数字 3 桁)。"""
    if x == 0:
        return "0.00"
    ax = abs(x)
    if 0.01 <= ax < 1.0e6:
        return f"{x:.2f}"
    return f"{x:.3g}"


def render(q, unit, embed=False):
    """量 ``q`` を単位 ``unit`` で書く。``embed`` なら ``"0.00 [Mm]"`` の形。"""
    factor = UNIT_TABLE[unit][1]
    text = " ".join(format_number(v / factor) for v in q.values)
    if embed:
        text += f" [{unit}]"
    return text


def parse_cell(text, header_unit=None):
    """シートのセルを量として読む。

    Returns
    -------
    (Quantity, unit, embedded) or None
        セルに ``[単位]`` があればそれを、無ければ見出しの単位を使う。
        数として読めない、または単位が分からなければ ``None``。
    """
    text = (text or "").strip()
    if not text:
        return None
    body, unit = split_header(text)
    embedded = unit is not None
    if not embedded:
        unit = header_unit
    if unit not in UNIT_TABLE:
        return None
    try:
        nums = tuple(float(w) for w in body.split())
    except ValueError:
        return None
    if not nums:
        return None
    family, factor = UNIT_TABLE[unit]
    return Quantity(family, tuple(v * factor for v in nums)), unit, embedded
