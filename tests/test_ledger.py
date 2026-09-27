"""r2d2plus-ledger (pyR2D2.ledger) のテスト。

Google には一切接続しない。シートは gspread.Worksheet の最小限の API を真似た
:class:`FakeWorksheet` で置き換える。
"""

import datetime
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

from pyR2D2.ledger import cli, collect, sheet, unsent_hint
from pyR2D2.ledger.units import Quantity, choose_unit, parse_cell, render

# 利用者の実際のシート (sunspot) の見出し行
REAL_HEADER = [
    "Case ID", "Mstar", "(ix,jx,kx)", "xmin [Mm]", "xmax [Mm]", "ymin", "ymax",
    "zmin", "zmax", "uni", "dx [km]", "m ray", "dtout [s]", "dtout_tau [s]", "al",
    "RSST", "Om", "Gemetry", "origin", "update time", "Note", "Finish", "HPCI",
    "富岳削除", "ISEE",
]  # fmt: skip
HUMAN = ["Note", "Finish", "HPCI", "富岳削除", "ISEE"]

RSUN_R2D2 = 0.6959894528e11


# ---------------------------------------------------------------------------
# 偽のワークシート
# ---------------------------------------------------------------------------


def _parse_a1(ref):
    m = re.fullmatch(r"([A-Z]+)(\d+)", ref)
    assert m, f"単一セルの範囲でない: {ref}"
    col = 0
    for ch in m.group(1):
        col = col * 26 + ord(ch) - ord("A") + 1
    return int(m.group(2)), col


class FakeWorksheet:
    """gspread.Worksheet のうち台帳が使う部分だけを持つ偽物。"""

    def __init__(self, rows, row_count=100, col_count=None):
        self.cells = {}
        for r, row in enumerate(rows, start=1):
            for c, v in enumerate(row, start=1):
                if v != "":
                    self.cells[(r, c)] = v
        self.row_count = row_count
        self.col_count = col_count or max([len(r) for r in rows] + [26])
        self.updated = []  # batch_update で書かれた (行, 列)

    def get_all_values(self):
        if not self.cells:
            return []
        nr = max(r for r, _ in self.cells)
        nc = max(c for _, c in self.cells)
        return [
            [self.cells.get((r, c), "") for c in range(1, nc + 1)]
            for r in range(1, nr + 1)
        ]

    def get(self, row, col):
        return self.cells.get((row, col), "")

    def add_cols(self, n):
        self.col_count += n

    def add_rows(self, n):
        self.row_count += n

    def batch_update(self, data):
        for item in data:
            row, col = _parse_a1(item["range"])
            assert row <= self.row_count and col <= self.col_count
            (value,) = item["values"][0]
            self.updated.append((row, col))
            if value == "":
                self.cells.pop((row, col), None)
            else:
                self.cells[(row, col)] = value

    def header(self):
        return self.get_all_values()[0]

    def col(self, name):
        """見出しの単位を除いた名前で列番号を探す。"""
        for c, text in enumerate(self.header(), start=1):
            if sheet.managed_name(text) == name or text == name:
                return c
        raise KeyError(name)


# ---------------------------------------------------------------------------
# 合成のランディレクトリ
# ---------------------------------------------------------------------------


def _write_params(param_dir, nx=32, ny=16, nz=16, xmin_off=-5.44e8, xmax_off=0.7e8,
                  dtout=60.0, dtout_tau=30.0, tend=3.0e4, geometry="Cartesian"):  # fmt: skip
    lines = [
        "R2D2 version 2.0",
        "       2 xdcheck i", "       2 ydcheck i", "       2 zdcheck i",
        "       2 margin i",
        f"{nx:8d} nx i", f"{ny:8d} ny i", f"{nz:8d} nz i",
        "       2 npe i", "       1 ix0 i", "       2 jx0 i", "       1 kx0 i",
        f"{RSUN_R2D2 + xmax_off:.10E} xmax d", f"{RSUN_R2D2 + xmin_off:.10E} xmin d",
        "0.6144000000E+09 ymax d", "0.0000000000E+00 ymin d",
        "0.6144000000E+09 zmax d", "0.0000000000E+00 zmin d",
        f"{dtout:.10E} dtout d", f"{tend:.10E} tend d",
        "       0 swap i",
        f"{dtout_tau:.10E} dtout_tau d",
        "0.0000000000E+00 omfac d",
        f"{RSUN_R2D2:.10E} rstar d",
        "0.1988999963E+34 mstar d",
        "0.1000000000E+01 potential_alpha d",
        "astana server c",
        "T rte l",
        "       1 rte_num_groups i",
        "F ununiform_flag l",
        f"{geometry} geometry c",
    ]  # fmt: skip
    (param_dir / "params.dac").write_text("\n".join(lines) + "\n")
    # back.dac: head, x, y, z, 1 次元量 20 本, tail (Fortran の順序なし書き出し)
    ixg, jxg, kxg = nx + 4, 2 * ny + 4, nz + 4
    x = RSUN_R2D2 + xmin_off + 2.4e6 * (np.arange(ixg) - 2 + 0.5)
    arrays = [x, np.zeros(jxg), np.zeros(kxg)]
    for name in range(20):
        arrays.append(np.ones(ixg) * (2.0 if name == 19 else 1.0))  # xi != 1 → RSST
    with open(param_dir / "back.dac", "wb") as f:
        np.array([0], "<i4").tofile(f)
        for a in arrays:
            a.astype("<f8").tofile(f)
        np.array([0], "<i4").tofile(f)
    (param_dir / "nd.dac").write_text("       3       0\n")


def make_run(root, caseid, status="normal_end", time=86400.0 * 3, step=1200,
             updated="2026-09-27T08:53:55+0900", parent=None, t_end=100 * 86400.0,
             modifiers=None, params=True, dirty=True):  # fmt: skip
    """R2D2plus のランを模した小さなディレクトリを作る。"""
    run = Path(root) / caseid
    param = run / "data" / "param"
    param.mkdir(parents=True)
    if params:
        _write_params(param)
    (param / "run_summary.toml").write_text(
        "[status]\n"
        f'status = "{status}"\n'
        f"step = {step}\n"
        f"time = {time!r}\n"
        f'updated = "{updated}"\n'
        "world_size = 4\n"
        'commit = "457b4ef2f5396027484768084c81cedfeb8faae1"\n'
        f"dirty = {'true' if dirty else 'false'}\n"
        'backend = "openmp"\n'
        "\n[run]\n"
        'geometry = "cartesian"\nyinyang = false\n'
        "nx = 64\nny = 64\nnz = 128\n"
        f"t_end = {t_end!r}\n"
    )
    if parent:
        (param / "origin.toml").write_text(
            f'[origin]\nparent = "/somewhere/run/{parent}"\nslot = "o"\nstep = 10\n'
        )
    seg = param / "runs" / "20260927-085355"
    seg.mkdir(parents=True)
    (seg / "effective_config.txt").write_text(
        "# 実効設定\n"
        "domain.xmin = 0\ndomain.xmax = 6.144e+08\n"
        "domain.ymin = 0\ndomain.ymax = 6.144e+08\n"
        f"domain.zmin = {RSUN_R2D2 - 5.44e8:.6e}\ndomain.zmax = {RSUN_R2D2 + 0.7e8:.6e}\n"
        f"time.t_end = {t_end:g}\n"
        "output.dtout = 60\noutput.dtout_tau = 30\n"
        'output.format = "compressed"\n'
    )
    (seg / "run.log").write_text(
        "### NORMAL END ###\n"
        "wall time: time loop 1.5 s for 6 steps (250.00 ms/step) [rank 0]\n"
    )
    for slot, st in (("o", step - 10), ("e", step)):
        d = run / "data" / "restart" / slot
        d.mkdir(parents=True)
        mods = modifiers if slot == "e" else None
        text = f"[restart]\nstep = {st}\n"
        if mods:
            text += "\n[modifiers]\napplied = [" + ", ".join(f'"{m}"' for m in mods) + "]\n"
        (d / "meta.toml").write_text(text)
    return run


@pytest.fixture
def project(tmp_path):
    runs = tmp_path / "sunspot" / "run"
    runs.mkdir(parents=True)
    make_run(runs, "d001", modifiers=["add_uniform_field|b_vertical=200"])
    make_run(runs, "d003", parent="d001", status="running",
             updated=(datetime.datetime.now().astimezone()
                      - datetime.timedelta(hours=1)).strftime("%Y-%m-%dT%H:%M:%S%z"))  # fmt: skip
    return runs


def real_sheet():
    """実シートを模した行: d001 に古い値と人のメモ、d009 は人のメモだけ。"""
    rows = [list(REAL_HEADER)]
    d001 = ["d001", "1.00", "256 256 256", "-5.44", "0.70", "0.00 [Mm]", "6.14 [Mm]",
            "0.00 [Mm]", "6.14 [Mm]", "T", "48.00 48.00", "T", "60.00", "30.00", "1.00",
            "F", "0.0", "Cartesian", "N/A", "2026-01-01 00:00:00",
            "初期の対流", "済", "hp200137", "", "共有"]  # fmt: skip
    rows.append(d001)
    rows += [[""] * len(REAL_HEADER) for _ in range(7)]  # 3..9 行目
    d009 = [""] * len(REAL_HEADER)
    d009[0] = "d009"
    d009[REAL_HEADER.index("Note")] = "計算していない"
    rows.append(d009)  # 10 行目
    return rows


def _human_snapshot(ws):
    cols = [ws.col(h) for h in HUMAN]
    return {(r, c): ws.get(r, c) for r in range(1, 40) for c in cols}


# ---------------------------------------------------------------------------
# 単位
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "family, value, unit",
    [
        ("time", 240.0, "s"),
        ("time", 86400.0 * 3, "d"),
        ("time", 3.0 * 3600, "h"),
        ("time", 30 * 86400.0 * 365.25, "yr"),
        ("length", -5.44e8, "Mm"),
        ("length", 2.4e6, "km"),
        ("length", 7.0e7, "km"),
        ("bytes", 777e6, "MB"),
        ("bytes", 167e9, "GB"),
        ("bytes", 3.2e12, "TB"),
    ],
)
def test_choose_unit(family, value, unit):
    assert choose_unit(family, [value]) == unit
    # 表示値は [1, 1000) に入る
    shown = float(render(Quantity(family, (value,)), unit))
    assert 1.0 <= abs(shown) < 1000.0


def test_choose_unit_keeps_current_within_band():
    # 0.70 Mm は km なら 700 だが、見出しが Mm なら保つ (0.1 <= 0.7 < 1e4)
    assert choose_unit("length", [7.0e7], current="Mm") == "Mm"
    # 100 日を s で表すと 8.64e6 で帯を外れるので替える
    assert choose_unit("time", [8.64e6], current="s") == "d"
    # 値が全部 0 なら今の単位
    assert choose_unit("time", [0.0], current="min") == "min"
    assert choose_unit("time", [0.0]) == "s"


def test_parse_cell_header_and_embedded_units():
    q, unit, embedded = parse_cell("0.00 [Mm]")
    assert (q.family, unit, embedded) == ("length", "Mm", True)
    q, unit, embedded = parse_cell("6.14", "Mm")
    assert q.values == (6.14e8,) and not embedded
    q, _, _ = parse_cell("48.00 24.00", "km")
    assert q.values == (4.8e6, 2.4e6)
    assert parse_cell("abc", "s") is None
    assert parse_cell("12", None) is None


# ---------------------------------------------------------------------------
# 収集
# ---------------------------------------------------------------------------


def test_collect_synthetic_run(project):
    rec = collect.collect_run(project / "d001", server="astana")
    v = rec.values
    assert v["(ix,jx,kx)"] == "32 32 16"
    assert v["Geometry"] == "Cartesian"
    assert v["Mstar"] == "1.00"
    assert v["xmin"].family == "length"
    assert v["xmin"].values[0] == pytest.approx(-5.44e8, rel=1e-6)
    assert v["dx"].values == pytest.approx((2.4e6, 2.4e6))
    assert v["RSST"] == "T"
    assert v["状態"] == "正常終了"
    assert v["commit"] == "457b4ef2-dirty"
    assert v["形式"] == "compressed"
    assert v["加工"] == "add_uniform_field|b_vertical=200"
    assert v["ms/step"] == "250.00"
    assert v["ranks"] == "4"
    assert v["origin"] == "N/A"
    assert v["容量"].family == "bytes"
    d003 = collect.collect_run(project / "d003", server="astana")
    assert d003.values["origin"] == "d001 (o)"
    assert d003.values["状態"] == "走行中"


def test_stale_running_is_unknown(tmp_path):
    run = make_run(tmp_path, "d002", status="running",
                   updated="2026-09-20T00:00:00+0900")  # fmt: skip
    now = datetime.datetime(2026, 9, 27, tzinfo=datetime.timezone.utc)
    assert collect.collect_run(run, now=now).values["状態"] == "不明"
    now = datetime.datetime(2026, 9, 20, 5, tzinfo=datetime.timezone.utc)
    assert collect.collect_run(run, now=now).values["状態"] == "走行中"


def test_run_without_params_uses_summary_and_config(tmp_path):
    run = make_run(tmp_path, "d004", params=False)
    v = collect.collect_run(run, server="astana").values
    # C++ の (nx, ny, nz) = (64, 64, 128) を (ix 鉛直, jx, kx) に並べ替える
    assert v["(ix,jx,kx)"] == "128 64 64"
    assert v["dtout"] == Quantity("time", (60.0,))
    assert v["Mstar"] is None and v["dx"] is None
    # rstar が分からないので鉛直範囲は空
    assert v["xmin"] is None


# ---------------------------------------------------------------------------
# シートへの送信 (偽のシート)
# ---------------------------------------------------------------------------


def test_push_protects_human_columns(project):
    ws = FakeWorksheet(real_sheet())
    before = _human_snapshot(ws)
    plan = cli.push(project, worksheet=ws, server="astana", out=lambda *a: None)
    assert {r.caseid for r in plan.written} == {"d001", "d003"}

    # 人の列は 1 セルも書かれていない
    human_cols = {ws.col(h) for h in HUMAN}
    assert not [rc for rc in ws.updated if rc[1] in human_cols]
    assert _human_snapshot(ws) == before
    assert ws.get(10, ws.col("Note")) == "計算していない"

    # 見出しは既存の位置のまま、無い管理列は右端 (ISEE の次) から足された
    header = ws.header()
    assert header[: len(REAL_HEADER)][17] == "Gemetry"
    assert header[len(REAL_HEADER)] == "Server"
    assert header[len(REAL_HEADER) + 1] == "状態"
    assert ws.col_count >= len(header)

    # 行は Case ID の番号 (d001 → 2 行目、d003 → 4 行目)
    assert ws.get(2, 1) == "d001" and ws.get(4, 1) == "d003"
    assert ws.get(2, ws.col("Server")) == "astana"
    assert ws.get(2, ws.col("Geometry")) == "Cartesian"
    assert ws.get(2, ws.col("(ix,jx,kx)")) == "32 32 16"
    assert ws.get(4, ws.col("状態")) == "走行中"
    assert ws.get(4, ws.col("origin")) == "d001 (o)"
    assert ws.get(2, ws.col("加工")) == "add_uniform_field|b_vertical=200"
    # 単位は見出しに入り、セルは数だけ (ymin の "[Mm]" 埋め込みは見出しへ移る)
    assert header[ws.col("ymax") - 1] == "ymax [Mm]"
    assert ws.get(2, ws.col("ymax")) == "6.14"
    assert header[ws.col("dtout") - 1] == "dtout [s]"
    assert header[ws.col("到達時刻") - 1] == "到達時刻 [d]"
    assert ws.get(2, ws.col("到達時刻")) == "3.00"
    # d009 (人のメモだけの行) の Case ID もそのまま
    assert ws.get(10, 1) == "d009"

    # 送信の記録が残り、status は未送信なし
    for c in ("d001", "d003"):
        assert (project / c / "data" / "param" / "ledger_sent.toml").exists()
    lines = []
    cli.cmd_status(type("A", (), {"path": project, "server": "astana"})(), out=lines.append)
    assert all(line.rstrip().endswith("-") for line in lines[1:])

    # 2 回目は何も変わらない (容量は ledger_sent.toml の分だけ増えるので測らない)
    ws.updated.clear()
    cli.push(project, worksheet=ws, server="astana", du=False, out=lambda *a: None)
    assert ws.updated == []


def test_push_skips_other_server_unless_forced(project):
    rows = [REAL_HEADER + ["Server"]]
    rows.append(["d001"] + [""] * (len(REAL_HEADER) - 1) + ["fugaku"])
    ws = FakeWorksheet(rows)
    msgs = []
    plan = cli.push(project, worksheet=ws, server="astana", out=msgs.append)
    assert [c for c, _ in plan.skipped] == ["d001"]
    assert any("fugaku" in m for m in msgs)
    assert not [rc for rc in ws.updated if rc[0] == 2]
    assert not (project / "d001" / "data" / "param" / "ledger_sent.toml").exists()

    plan = cli.push(project, worksheet=ws, server="astana", force=True, out=lambda *a: None)
    assert "d001" in {r.caseid for r in plan.written}
    assert ws.get(2, ws.col("Server")) == "astana"


def test_unit_change_converts_other_rows(project):
    header = ["Case ID", "t_end [s]", "ymin", "Note"]
    rows = [header]
    rows += [[""] * 4 for _ in range(19)]
    rows.append(["d020", "8640000", "1.00 [Mm]", "他のサーバのラン"])  # 21 行目
    ws = FakeWorksheet(rows)
    cli.push(project, worksheet=ws, server="astana", out=lambda *a: None)
    # 送る 2 本の t_end は 100 日 → 列の単位は d
    assert ws.get(1, 2) == "t_end [d]"
    assert ws.get(2, 2) == "100.00"
    # 手元に無い d020 の値も換算される
    assert ws.get(21, 2) == "100.00"
    assert ws.get(1, 3) == "ymin [Mm]"
    assert ws.get(21, 3) == "1.00"
    assert ws.get(21, 4) == "他のサーバのラン"


def test_mixed_families_embed_units(project, tmp_path):
    ws = FakeWorksheet([["Case ID", "ymax"], ["", ""], ["", ""], ["", ""],
                        ["d004", "90.00 [deg]"]])  # fmt: skip
    cli.push(project, worksheet=ws, server="astana", out=lambda *a: None)
    assert ws.get(1, 2) == "ymax"
    assert ws.get(2, 2) == "6.14 [Mm]"
    assert ws.get(5, 2) == "90.00 [deg]"


def test_empty_sheet_gets_full_header(project):
    ws = FakeWorksheet([], col_count=5, row_count=1)
    cli.push(project, worksheet=ws, server="astana", out=lambda *a: None)
    names = [sheet.managed_name(h) for h in ws.header()]
    assert names == collect.COLUMN_NAMES
    assert ws.col_count >= len(collect.COLUMNS) and ws.row_count >= 4


def test_only_and_non_case_dirs(project):
    make_run(project, "test_run")
    ws = FakeWorksheet([list(REAL_HEADER)])
    plan = cli.push(project, worksheet=ws, server="astana", only=["d003", "test_run"],
                    out=lambda *a: None)  # fmt: skip
    assert [r.caseid for r in plan.written] == ["d003"]
    assert [c for c, _ in plan.skipped] == ["test_run"]
    assert ws.get(2, 1) == ""


def test_no_du_keeps_capacity_cell(project):
    ws = FakeWorksheet([list(REAL_HEADER)])
    cli.push(project, worksheet=ws, server="astana", out=lambda *a: None)
    col = ws.col("容量")
    before = ws.get(2, col)
    assert before
    ws.updated.clear()
    cli.push(project, worksheet=ws, server="astana", du=False, out=lambda *a: None)
    assert ws.get(2, col) == before
    assert (2, col) not in ws.updated


def test_dry_run_does_not_touch_google_or_disk(project, monkeypatch):
    monkeypatch.setitem(sys.modules, "gspread", None)  # import すれば失敗する
    out = []
    plan = cli.push(project, dry_run=True, server="astana", out=out.append)
    assert plan.cells
    assert any("接続していない" in line for line in out)
    assert any(line.startswith("A2\td001") for line in out)
    assert not (project / "d001" / "data" / "param" / "ledger_sent.toml").exists()


def test_push_without_gspread_gives_clear_error(project, monkeypatch):
    monkeypatch.setitem(sys.modules, "gspread", None)
    with pytest.raises(SystemExit, match="gspread"):
        cli.push(project, server="astana", out=lambda *a: None)


def test_export_csv(project, tmp_path):
    out = tmp_path / "ledger.csv"
    assert cli.main(["export", str(project), "-o", str(out), "--server", "astana"]) == 0
    import csv

    rows = list(csv.reader(out.open(encoding="utf-8")))
    assert rows[0][0] == "Case ID"
    assert "dtout [s]" in rows[0] and "xmin [Mm]" in rows[0]
    assert [r[0] for r in rows[1:]] == ["d001", "d003"]


def test_set_cells_gspread_uses_header_names(project, monkeypatch):
    """旧来の set_cells_gspread も Server を落とさず、人の列に触らない。"""
    import pyR2D2.write.google.google as g

    ws = FakeWorksheet(real_sheet())
    before = _human_snapshot(ws)

    class Book:
        sheet1 = ws

    class Client:
        def open(self, name):
            assert name == "sunspot"
            return Book()

    monkeypatch.setattr(g, "init_gspread", lambda key, project: Client())
    monkeypatch.setattr(collect, "default_server", lambda: "astana")
    data = type("D", (), {"datadir": project / "d003" / "data"})()
    g.set_cells_gspread(data, json_key="dummy", project="sunspot")
    assert ws.get(4, ws.col("Server")) == "astana"
    assert _human_snapshot(ws) == before


# ---------------------------------------------------------------------------
# pyR2D2.Data を開いたときの 1 行の案内
# ---------------------------------------------------------------------------


def test_unsent_hint(project, tmp_path):
    datadir = project / "d001" / "data"
    msg = unsent_hint(datadir)
    assert msg == f"台帳へ未送信の変更があります: r2d2plus-ledger push {project.resolve()}"
    rec = collect.collect_run(project / "d001", server="astana", du=False)
    collect.write_sent(rec, "sunspot")
    assert unsent_hint(datadir) is None
    # run_summary が更新されれば再び出る
    summary = datadir / "param" / "run_summary.toml"
    summary.write_text(summary.read_text().replace("step = 1200", "step = 1300"))
    assert unsent_hint(datadir) is not None
    # 壊れた ledger_sent.toml でも例外を出さない
    (datadir / "param" / "ledger_sent.toml").write_text("[[[ broken")
    assert unsent_hint(datadir) is None
    # Fortran 版 (run_summary.toml が無い) は黙る
    assert unsent_hint(tmp_path / "nothing" / "data") is None


def test_data_prints_hint_only_for_r2d2plus(capsys):
    """テスト用の Fortran 形式のデータを開いても案内は出ない。"""
    import pyR2D2

    pyR2D2.Data(Path(__file__).parent / "data")
    assert "台帳" not in capsys.readouterr().out


# ---------------------------------------------------------------------------
# 手元の実データ (あれば。小さな param ファイルだけ tmp に写して使う)
# ---------------------------------------------------------------------------

_SCRATCH = Path(
    "/home/sc/c0234hotta/tmp/claude-40234/-scr-a000-c0234hotta-Repository-R2D2plus/"
    "59364ee3-4d14-480c-ac49-4376767e1fd5/scratchpad/proj1/run"
)
_COMP3D = Path("/scr/a000/c0234hotta/io_fixture_20260927/comp_3d")
_FORTRAN = Path("/scr/a000/c0234hotta/odf/run/d001")


def _copy_small(src_run, dst_run):
    """ラン出力のうち台帳が読む小さなファイルだけを写す (元には書かない)。"""
    src, dst = src_run / "data", dst_run / "data"
    (dst / "param").mkdir(parents=True)
    for f in (src / "param").iterdir():
        if f.is_file() and f.suffix in (".dac", ".toml"):
            shutil.copy2(f, dst / "param" / f.name)
    runs = src / "param" / "runs"
    if runs.is_dir():
        for seg in runs.iterdir():
            d = dst / "param" / "runs" / seg.name
            d.mkdir(parents=True)
            for name in ("effective_config.txt", "run.log"):
                if (seg / name).is_file():
                    shutil.copy2(seg / name, d / name)
    rst = src / "restart"
    if rst.is_dir():
        for slot in rst.iterdir():
            if (slot / "meta.toml").is_file():
                (dst / "restart" / slot.name).mkdir(parents=True)
                shutil.copy2(slot / "meta.toml", dst / "restart" / slot.name)
    if (src / "cont_log.txt").is_file():
        shutil.copy2(src / "cont_log.txt", dst / "cont_log.txt")
    return dst_run


@pytest.mark.skipif(not _SCRATCH.is_dir(), reason="手元の R2D2plus ラン群が無い")
def test_real_cpp_runs(tmp_path):
    runs = tmp_path / "proj1" / "run"
    for src in sorted(_SCRATCH.iterdir()):
        if (src / "data" / "param").is_dir():
            _copy_small(src, runs / src.name)
    ws = FakeWorksheet([list(REAL_HEADER)])
    plan = cli.push(runs, worksheet=ws, server="astana", out=lambda *a: None)
    assert len(plan.written) == len(list(runs.iterdir()))
    assert ws.get(3, ws.col("origin")) == "d001 (o)"  # d002 の親
    assert ws.get(6, ws.col("加工")).startswith("add_uniform_field")  # d005
    assert ws.get(6, ws.col("状態")) == "異常終了"
    assert ws.get(2, ws.col("(ix,jx,kx)")) == "24 6 10"


@pytest.mark.skipif(not _COMP3D.is_dir(), reason="comp_3d の fixture が無い")
def test_real_compressed_run(tmp_path):
    run = _copy_small(_COMP3D, tmp_path / "proj" / "run" / "comp_3d")
    v = collect.collect_run(run, server="astana").values
    assert v["形式"] == "compressed"
    assert v["(ix,jx,kx)"] == "128 128 128"
    assert v["ranks"] == "4"
    assert v["xmin"].values[0] == pytest.approx(-5.44e8, rel=1e-3)
    assert collect.project_name(run) == "proj"


@pytest.mark.skipif(not _FORTRAN.is_dir(), reason="Fortran 版のランが無い")
def test_real_fortran_run(tmp_path):
    run = _copy_small(_FORTRAN, tmp_path / "odf" / "run" / "d001")
    rec = collect.collect_run(run, server="astana", du=False)
    v = rec.values
    assert rec.fortran
    assert v["(ix,jx,kx)"] == "256 256 256"
    assert v["状態"] is None and v["commit"] is None
    assert v["dx"].values[0] == pytest.approx(2.4e6, rel=1e-6)
    assert render(v["xmin"], "Mm") == "-5.44"
    assert unsent_hint(run / "data") is None


def test_console_script_entry_point():
    import importlib
    import tomllib

    with open(Path(__file__).resolve().parents[1] / "pyproject.toml", "rb") as f:
        scripts = tomllib.load(f)["project"]["scripts"]
    module, _, func = scripts["r2d2plus-ledger"].partition(":")
    assert callable(getattr(importlib.import_module(module), func))


# ---------------------------------------------------------------------------
# 対になる列は単位を共有する
# ---------------------------------------------------------------------------


def _record(caseid, **values):
    rec = collect.RunRecord(caseid=caseid, run_dir=Path("/nonexistent") / caseid)
    rec.values = {name: None for name in collect.COLUMN_NAMES}
    rec.values["Case ID"] = caseid
    rec.values.update(values)
    return rec


def test_pairs_share_unit_in_export():
    # odf の Fortran ランと同じ: xmin -5.44 Mm、xmax +0.70 Mm (単独なら 700 km)
    L = lambda v: Quantity("length", (v,))  # noqa: E731
    T = lambda v: Quantity("time", (v,))  # noqa: E731
    recs = [
        _record("d001", xmin=L(-5.44e8), xmax=L(7.0e7), 到達時刻=T(240.0),
                t_end=T(3 * 86400.0), dtout=T(1.0e9), dtout_tau=T(30.0)),  # fmt: skip
        _record("d002", xmin=L(-5.44e8), xmax=L(1.0e8), 到達時刻=T(600.0),
                t_end=T(3 * 86400.0), dtout=T(1.0e9), dtout_tau=T(30.0)),  # fmt: skip
    ]
    headers, rows = sheet.render_records(recs)
    col = {sheet.managed_name(h): i for i, h in enumerate(headers)}
    assert headers[col["xmin"]] == "xmin [Mm]"
    assert headers[col["xmax"]] == "xmax [Mm]"
    assert rows[0][col["xmax"]] == "0.70" and rows[1][col["xmax"]] == "1.00"
    # 到達時刻と t_end は同じ単位 (比べられるように)
    u1 = headers[col["到達時刻"]].split("[")[1]
    u2 = headers[col["t_end"]].split("[")[1]
    assert u1 == u2
    # dtout と dtout_tau は独立
    assert headers[col["dtout"]] == "dtout [yr]"
    assert headers[col["dtout_tau"]] == "dtout_tau [s]"


def test_pair_with_different_header_units_is_unified(project):
    header = ["Case ID", "xmin [Mm]", "xmax [km]", "Note"]
    rows = [header] + [[""] * 4 for _ in range(19)]
    rows.append(["d020", "-5.44", "700.00", "他のサーバ"])  # 21 行目 (手元に無い)
    ws = FakeWorksheet(rows)
    cli.push(project, worksheet=ws, server="astana", out=lambda *a: None)
    assert ws.get(1, 2) == "xmin [Mm]"
    assert ws.get(1, 3) == "xmax [Mm]"
    assert ws.get(21, 2) == "-5.44"
    assert ws.get(21, 3) == "0.70"  # 700 km を Mm へ換算
    assert ws.get(2, 3) == "0.70"
    assert ws.get(21, 4) == "他のサーバ"


def test_pair_switches_together(project):
    # 両方 [s] の見出しに 100 日の値 → 帯を外れるので組で d に替える
    header = ["Case ID", "到達時刻 [s]", "t_end [s]"]
    rows = [header] + [[""] * 3 for _ in range(19)]
    rows.append(["d020", "259200", "8640000"])
    ws = FakeWorksheet(rows)
    cli.push(project, worksheet=ws, server="astana", out=lambda *a: None)
    assert ws.get(1, 2) == "到達時刻 [d]" and ws.get(1, 3) == "t_end [d]"
    assert ws.get(21, 2) == "3.00" and ws.get(21, 3) == "100.00"
