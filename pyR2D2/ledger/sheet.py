"""シートの今の値と集めたランの値から、書き換えるセルの一覧を作る。

人の列の保護
------------
* 管理列 (:data:`collect.COLUMNS`) は見出しの**単位を除いた名前**で探す
  (``dtout [s]`` も ``dtout [min]`` も ``dtout``)。それ以外の列は人の列で、
  読み書きの対象にしない。
* 書き込みは管理列のセル 1 個ずつの範囲 (``B2`` など) を ``batch_update`` に
  並べる。行全体を読んで書き戻すことはしないので、人の列は送信のあいだに人が
  書き換えても壊れない。
* シートに無い管理列は、見出し行の右端 (最後の空でない見出しの次) から足す。
  見出し行が空なら、管理列の見出しを左から全部書く。
* 行は Case ID ``dNNN`` の番号で決まる (d001 → 2 行目)。その行の ``Server`` が
  送り手と違う、または ``Case ID`` が違うときは警告して飛ばす (``force`` で上書き)。

単位が変わるとき
----------------
量の列は、送る行の値と、**送らない行 (他のサーバのラン・手元に無いラン) の
セルを今の見出しの単位 (またはセル内の ``[単位]``) で読んだ値**を合わせて単位を
選ぶ (:func:`units.choose_unit`)。単位が変わったら見出しを書き換え、送らない行の
セルも読んだ値を新しい単位へ換算して書き直す。元の値は表示の 2 桁までしか
残っていないので、その丸めの分だけ誤差が入る (送る行は手元のファイルから
計算し直すので誤差は無い)。数として読めないセルは触らず警告する。

1 つの列に種類の違う量が混ざる (直交座標の ``ymin`` [Mm] と球の ``ymin`` [deg]
など) ときは、見出しを単位なしにして、各セルに ``"0.00 [Mm]"`` の形で単位を書く。
"""

from collections import Counter
from dataclasses import dataclass, field

from .collect import COLUMNS, QUANTITY_COLUMNS, case_row
from .units import Quantity, UNITS, choose_unit, parse_cell, render, split_header

_ALIASES = {}
for _name, _aliases, _ in COLUMNS:
    _ALIASES[_name] = _name
    for _a in _aliases:
        _ALIASES[_a] = _name


def col_letter(col):
    """1 始まりの列番号を A1 表記の列名にする (1 → A, 27 → AA)。"""
    s = ""
    while col > 0:
        col, r = divmod(col - 1, 26)
        s = chr(ord("A") + r) + s
    return s


def a1(row, col):
    return f"{col_letter(col)}{row}"


def managed_name(header_text):
    """見出しが管理列なら正規の名前、そうでなければ ``None``。"""
    name, _ = split_header(header_text)
    name = " ".join(name.split())
    return _ALIASES.get(name)


@dataclass
class Plan:
    """書き換えの計画。``cells`` は ``(行, 列, 文字列)`` (1 始まり)。"""

    cells: list = field(default_factory=list)
    written: list = field(default_factory=list)  # 送る RunRecord
    skipped: list = field(default_factory=list)  # (caseid, 理由)
    warnings: list = field(default_factory=list)
    n_cols: int = 0  # 必要な列数
    n_rows: int = 0  # 必要な行数
    headers: dict = field(default_factory=dict)  # 列名 → 見出しの文字列


def _cell(values, row, col):
    """``values`` (get_all_values の結果) の (行, 列) (1 始まり)。範囲外は空。"""
    if row - 1 < len(values):
        r = values[row - 1]
        if col - 1 < len(r):
            return r[col - 1]
    return ""


def decide_quantity_column(name, header_text, new, existing):
    """量の列 1 本の見出しとセルの文字列を決める。

    Parameters
    ----------
    name : str
        管理列の名前
    header_text : str or None
        今の見出し (無ければ ``None``)
    new : dict
        キー → 送る値 (:class:`Quantity` / 文字列 / ``None``)
    existing : dict
        キー → 送らない行の今のセルの文字列

    Returns
    -------
    header : str
    new_texts : dict
        キー → 送る行のセルの文字列
    rewrites : dict
        キー → 送らない行のうち書き直すセルの文字列
    warnings : list of str
    """
    _, header_unit = split_header(header_text) if header_text else (name, None)
    parsed = {}
    unparsed = []
    for key, text in existing.items():
        r = parse_cell(text, header_unit)
        if r is None:
            if text.strip():
                unparsed.append(key)
        else:
            parsed[key] = r

    quantities = [q for q in new.values() if isinstance(q, Quantity)]
    quantities += [r[0] for r in parsed.values()]
    families = {q.family for q in quantities}

    def current_unit(family):
        if header_unit and header_unit in {u for u, _, _ in UNITS[family]}:
            return header_unit
        units = Counter(r[1] for r in parsed.values() if r[0].family == family)
        return units.most_common(1)[0][0] if units else None

    chosen = {}
    for family in families:
        vals = [v for q in quantities if q.family == family for v in q.values]
        chosen[family] = choose_unit(family, vals, current_unit(family))

    embed = len(families) > 1
    if embed:
        header = name
    elif families:
        header = f"{name} [{chosen[next(iter(families))]}]"
    else:
        header = header_text if header_text else name

    new_texts = {}
    for key, q in new.items():
        if isinstance(q, Quantity):
            new_texts[key] = render(q, chosen[q.family], embed)
        else:
            new_texts[key] = "" if q is None else str(q)

    rewrites = {}
    for key, (q, unit, was_embedded) in parsed.items():
        if unit != chosen[q.family] or was_embedded != embed:
            rewrites[key] = render(q, chosen[q.family], embed)

    warnings = []
    header_changed = split_header(header)[1] != header_unit
    if unparsed and header_changed:
        warnings.append(
            f"列 {name}: 単位を {header_unit} → {split_header(header)[1]} に変えたが、"
            f"数として読めないセルが {len(unparsed)} 個あり、そのままにした"
        )
    return header, new_texts, rewrites, warnings


def render_records(records):
    """シートを使わずに、ランの並びから見出しと各行の文字列を作る (export 用)。"""
    headers = []
    rows = [[None] * len(COLUMNS) for _ in records]
    for j, (name, _, is_qty) in enumerate(COLUMNS):
        new = {i: r.values.get(name) for i, r in enumerate(records)}
        if is_qty:
            header, texts, _, _ = decide_quantity_column(name, None, new, {})
        else:
            header = name
            texts = {i: "" if v is None else str(v) for i, v in new.items()}
        headers.append(header)
        for i in range(len(records)):
            rows[i][j] = texts[i]
    return headers, rows


def plan_updates(values, records, server, force=False):
    """シートの今の値 ``values`` (``get_all_values()``) に対する書き換えを計画する。"""
    plan = Plan()
    header_row = list(values[0]) if values else []
    while header_row and not header_row[-1].strip():
        header_row.pop()

    # 管理列の位置
    positions = {}
    header_texts = {}
    for col, text in enumerate(header_row, start=1):
        name = managed_name(text)
        if name is not None and name not in positions:
            positions[name] = col
            header_texts[name] = text
    next_col = len(header_row) + 1
    for name, _, _ in COLUMNS:
        if name not in positions:
            positions[name] = next_col
            next_col += 1
    plan.n_cols = next_col - 1

    # 送る行の決定
    rows = {}
    for rec in records:
        row = case_row(rec.caseid)
        if row is None:
            plan.skipped.append((rec.caseid, "Case ID が dNNN の形でない"))
            continue
        cur_case = _cell(values, row, positions["Case ID"]).strip()
        cur_server = _cell(values, row, positions["Server"]).strip()
        if cur_case and cur_case != rec.caseid and not force:
            plan.skipped.append(
                (rec.caseid, f"{row} 行目の Case ID が {cur_case} (--force で上書き)")
            )
            continue
        if cur_server and cur_server != server and not force:
            plan.skipped.append(
                (rec.caseid, f"Server が {cur_server} (送り手は {server}。--force で上書き)")
            )
            continue
        rows[row] = rec
        plan.written.append(rec)
    plan.n_rows = max([len(values)] + list(rows))

    other_rows = [r for r in range(2, len(values) + 1) if r not in rows]

    for name, _, _ in COLUMNS:
        col = positions[name]
        old_header = header_texts.get(name)
        new = {
            row: rec.values[name] for row, rec in rows.items() if name in rec.values
        }
        if name in QUANTITY_COLUMNS:
            # 送らない行と、この列を測っていない送る行 (--no-du の容量) は今の値のまま
            keep = other_rows + [r for r in rows if r not in new]
            existing = {r: _cell(values, r, col) for r in keep}
            header, texts, rewrites, warns = decide_quantity_column(
                name, old_header, new, existing
            )
            plan.warnings.extend(warns)
        else:
            header = old_header if old_header is not None else name
            texts = {r: "" if v is None else str(v) for r, v in new.items()}
            rewrites = {}
        plan.headers[name] = header
        if header != old_header:
            plan.cells.append((1, col, header))
        for row, text in list(texts.items()) + list(rewrites.items()):
            if text != _cell(values, row, col):
                plan.cells.append((row, col, text))
    return plan


def apply_plan(worksheet, plan, chunk=500):
    """計画を worksheet (gspread.Worksheet 互換) に書く。

    必要なら列・行を足してから、管理列のセルだけを ``batch_update`` で書く。
    """
    if plan.n_cols > worksheet.col_count:
        worksheet.add_cols(plan.n_cols - worksheet.col_count)
    if plan.n_rows > worksheet.row_count:
        worksheet.add_rows(plan.n_rows - worksheet.row_count)
    data = [{"range": a1(r, c), "values": [[text]]} for r, c, text in plan.cells]
    for i in range(0, len(data), chunk):
        worksheet.batch_update(data[i : i + chunk])
