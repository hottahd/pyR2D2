"""``r2d2plus-ledger``: ランの一覧 (台帳) を Google スプレッドシートへ送る。

.. code:: console

    $ r2d2plus-ledger status  ~/work/proj/run            # 要約と未送信の有無 (送らない)
    $ r2d2plus-ledger export  ~/work/proj/run -o ledger.csv
    $ r2d2plus-ledger push    ~/work/proj/run --dry-run  # 書き込むセルを表示するだけ
    $ r2d2plus-ledger push    ~/work/proj/run            # シートの管理列だけを更新
    $ r2d2plus-ledger pull    ~/work/proj/run            # シートで書き換えた目的などを 1 つずつ確かめて手元へ

シートはプロジェクトごとに 1 枚 (既定の名前はプロジェクトのディレクトリ名)。
認証は従来の :py:func:`pyR2D2.write.google.init_gspread` と同じサービスアカウント
(``~/json/`` の JSON 鍵)。
"""

import argparse
import csv
import glob
import os
import sys

from . import collect, note, sheet
from .units import Quantity, render, choose_unit


def _records(path, server, du, only=None):
    runs = collect.find_runs(path)
    if only:
        wanted = set(only)
        runs = [r for r in runs if r.name in wanted]
        missing = wanted - {r.name for r in runs}
        if missing:
            print("見つからないラン: " + ", ".join(sorted(missing)), file=sys.stderr)
    return [collect.collect_run(r, server=server, du=du) for r in runs]


def _fmt_time(q):
    if not isinstance(q, Quantity):
        return ""
    unit = choose_unit(q.family, list(q.values))
    return render(q, unit) + " " + unit


def cmd_status(args, out=print):
    """各ランの要約と、台帳へ未送信の変更があるかを表示する。シートには触らない。"""
    records = _records(args.path, args.server, du=False)
    rows = [("Case ID", "状態", "到達時刻", "step", "未送信")]
    for r in records:
        unsent = r.unsent()
        mark = "初回" if unsent == "never" else ("あり" if unsent else "-")
        rows.append(
            (
                r.caseid,
                r.values.get("状態") or "",
                _fmt_time(r.values.get("到達時刻")),
                r.values.get("step") or "",
                mark,
            )
        )
    widths = [max(_width(row[i]) for row in rows) for i in range(len(rows[0]))]
    for row in rows:
        out("  ".join(_pad(c, w) for c, w in zip(row, widths)).rstrip())
    return records


def _width(s):
    # 全角は 2 桁で数える (端末で列を揃えるため)
    return sum(2 if ord(ch) > 0x2E7F else 1 for ch in str(s))


def _pad(s, w):
    return str(s) + " " * (w - _width(s))


def cmd_export(args, out=print):
    """管理列を CSV に書く (単位は見出しに入る)。"""
    records = _records(args.path, args.server, du=not args.no_du)
    headers, rows = sheet.render_records(records)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerows(rows)
    out(f"{len(records)} 本のランを {args.output} に書いた")
    return headers, rows


def open_worksheet(sheet_name=None, sheet_id=None, json_key=None):
    """Google スプレッドシートの最初のワークシートを開く (ここでだけ gspread を読む)。"""
    try:
        import gspread  # noqa: F401
        from google.oauth2.service_account import Credentials  # noqa: F401
    except ImportError as e:
        raise SystemExit(
            "r2d2plus-ledger push には gspread と google-auth が要る "
            f"(pip install gspread google-auth): {e}"
        )
    from ..write.google.google import init_gspread

    if json_key is None:
        keys = sorted(glob.glob(os.path.join(os.environ["HOME"], "json", "*")))
        if not keys:
            raise SystemExit("~/json/ にサービスアカウントの JSON 鍵が無い (--json-key で指定)")
        json_key = keys[0]
    gc = init_gspread(json_key, sheet_name)
    book = gc.open_by_key(sheet_id) if sheet_id else gc.open(sheet_name)
    return book.sheet1


def push(
    path,
    worksheet=None,
    sheet_name=None,
    sheet_id=None,
    json_key=None,
    server=None,
    dry_run=False,
    force=False,
    only=None,
    du=True,
    reorder=True,
    out=print,
):
    """ランの管理列をシートへ送る。

    Parameters
    ----------
    path : str or pathlib.Path
        ``<project>/run`` (またはラン 1 本のディレクトリ)
    worksheet : object, optional
        gspread.Worksheet 互換のもの。省略時は ``sheet_name`` / ``sheet_id`` で開く。
        ``dry_run`` で ``worksheet`` が無ければ、空のシートに対する計画を表示する
        (Google には接続しない)。

    Returns
    -------
    plan : sheet.Plan
    """
    records = _records(path, server, du=du, only=only)
    label = sheet_id or sheet_name or collect.project_name(path)

    if worksheet is None and not dry_run:
        worksheet = open_worksheet(
            sheet_name=None if sheet_id else label, sheet_id=sheet_id, json_key=json_key
        )
    return push_records(
        records, worksheet, server, label, dry_run=dry_run, force=force, reorder=reorder, out=out
    )


def push_records(
    records, worksheet, server, label, dry_run=False, force=False, reorder=True, out=print
):
    """集めた :class:`collect.RunRecord` の並びを worksheet へ送る (:func:`push` の本体)。"""
    values = worksheet.get_all_values() if worksheet is not None else []
    plan = sheet.plan_updates(values, records, server, force=force)

    for caseid, reason in plan.skipped:
        out(f"警告: {caseid} を飛ばした: {reason}")
    for w in plan.warnings:
        out(f"警告: {w}")
    # 管理列を COLUMNS の順に並べ替える (人の列は動かさない)。書いた後の見出しで決める
    moves, final = sheet.column_moves(plan.header_after) if reorder else ([], [])

    if dry_run:
        if moves:
            out("列を並べ替える: " + ", ".join(h for h in final if h))
        if worksheet is None:
            out(f"(dry-run: シート {label} には接続していない。空のシートに対する計画)")
        for r, c, text in sorted(plan.cells):
            out(f"{sheet.a1(r, c)}\t{text}")
        out(f"(dry-run) {len(plan.cells)} セルを書く予定、{len(plan.written)} 本のラン")
        return plan

    sheet.apply_plan(worksheet, plan)
    if moves:
        sheet.apply_column_moves(worksheet, moves)
        out(f"管理列を決まった順に並べ替えた ({len(moves)} 回の移動)")
    for rec in plan.written:
        collect.write_sent(rec, label)
    out(f"{len(plan.written)} 本のランをシート {label} に送った ({len(plan.cells)} セル)")
    return plan


def pull(
    path,
    worksheet=None,
    sheet_name=None,
    sheet_id=None,
    json_key=None,
    server=None,
    dry_run=False,
    force=False,
    only=None,
    ask=None,
    out=print,
):
    """シートで書き換えた目的・プロジェクト・保存・結論を、1 つずつ確かめて手元のメモ (run_note.toml) へ戻す。

    ほかの列 (状態・時刻・格子など) は機械がランから読む値なので戻さない。シートのセルが空なら戻さない。

    Parameters
    ----------
    ask : callable, optional
        質問の文字列を受けて答えを返す (既定は ``input``)。``y`` で反映、``n`` で飛ばす、``q`` で終わる。

    Returns
    -------
    applied : list of sheet.Change
    """
    records = _records(path, server, du=False, only=only)
    label = sheet_id or sheet_name or collect.project_name(path)
    if ask is None:
        if not dry_run and not sys.stdin.isatty():
            raise SystemExit("pull は 1 つずつ y/n で聞くので端末から実行する (--dry-run なら表示だけ)")
        ask = input
    if worksheet is None:
        worksheet = open_worksheet(
            sheet_name=None if sheet_id else label, sheet_id=sheet_id, json_key=json_key
        )
    plan = sheet.plan_pull(worksheet.get_all_values(), records, server, force=force)
    for caseid, reason in plan.skipped:
        out(f"警告: {caseid} を飛ばした: {reason}")
    for w in plan.warnings:
        out(f"警告: {w}")
    if not plan.changes:
        out(f"シート {label} に、手元と違う目的・プロジェクト・保存・結論は無い")
        return []

    applied = []
    for i, c in enumerate(plan.changes, start=1):
        out(f"[{i}/{len(plan.changes)}] {c.record.caseid} の {c.column}")
        out(f"  手元  : {c.local or '(空)'}")
        out(f"  シート: {c.remote}")
        if c.column == "保存":
            cur = note.read(c.record.run_dir) or {}
            if "expires" not in c.fields and cur.get("expires"):
                out(f"  (手元の「{cur['expires']} まで残す」はシートに無いが、消さずに残る。消すなら手元で消す)")
            if "keep_parts" not in c.fields and cur.get("keep_parts"):
                parts = note.KEEP_PARTS_JOIN.join(cur["keep_parts"])
                out(f"  (手元の「{note.KEEP_PARTS_SEP}{parts}」はシートに無いが、消さずに残る。消すなら手元で消す)")
        if dry_run:
            continue
        while True:
            ans = ask("  手元のメモに反映しますか [y/n/q]: ").strip().lower()
            if ans in ("y", "n", "q"):
                break
        if ans == "q":
            break
        if ans == "y":
            p = note.update(c.record.run_dir, c.fields)
            applied.append(c)
            out(f"  反映した: {p}")
    if dry_run:
        out(f"(dry-run) {len(plan.changes)} か所が手元と違う。手元は書き換えていない")
    else:
        out(f"{len(applied)} / {len(plan.changes)} か所を手元のメモに反映した")
    return applied


def cmd_pull(args, out=print):
    only = [s.strip() for s in args.only.split(",")] if args.only else None
    return pull(
        args.path,
        sheet_name=args.sheet,
        sheet_id=args.sheet_id,
        json_key=args.json_key,
        server=args.server,
        dry_run=args.dry_run,
        force=args.force,
        only=only,
        out=out,
    )


def cmd_push(args, out=print):
    only = [s.strip() for s in args.only.split(",")] if args.only else None
    return push(
        args.path,
        sheet_name=args.sheet,
        sheet_id=args.sheet_id,
        json_key=args.json_key,
        server=args.server,
        dry_run=args.dry_run,
        force=args.force,
        only=only,
        du=not args.no_du,
        reorder=not args.keep_column_order,
        out=out,
    )


def build_parser():
    ap = argparse.ArgumentParser(
        prog="r2d2plus-ledger",
        description="ランの一覧 (台帳) を Google スプレッドシートへ送る",
    )
    sub = ap.add_subparsers(dest="command", required=True)

    def common(p):
        p.add_argument("path", help="<project>/run (またはラン 1 本のディレクトリ)")
        p.add_argument("--server", default=None, help="サーバ名 (既定: 前回送ったときの名前、無ければホスト名)")

    p = sub.add_parser("status", help="各ランの要約と未送信の有無 (送らない)")
    common(p)
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("export", help="管理列を CSV に書く")
    common(p)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--no-du", action="store_true", help="容量 (data/ の使用量) を測らない")
    p.set_defaults(func=cmd_export)

    p = sub.add_parser("push", help="シートの管理列を更新する")
    common(p)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--sheet", default=None, help="シート名 (既定: プロジェクト名)")
    g.add_argument("--sheet-id", default=None, help="シートの ID (URL の /d/ の後)")
    p.add_argument("--json-key", default=None, help="サービスアカウントの JSON 鍵")
    p.add_argument("--dry-run", action="store_true", help="書くセルを表示するだけ")
    p.add_argument("--force", action="store_true", help="Server の違う行も上書きする")
    p.add_argument("--only", default=None, help="送るラン (例: d001,d005)")
    p.add_argument("--no-du", action="store_true", help="容量 (data/ の使用量) を測らない")
    p.add_argument(
        "--keep-column-order",
        action="store_true",
        help="管理列を決まった順 (Case ID, プロジェクト, 目的, 保存, 容量, …) に並べ替えない",
    )
    p.set_defaults(func=cmd_push)

    p = sub.add_parser(
        "pull", help="シートで書き換えた目的・プロジェクト・保存・結論を、1 つずつ y/n で確かめて手元のメモへ戻す"
    )
    common(p)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--sheet", default=None, help="シート名 (既定: プロジェクト名)")
    g.add_argument("--sheet-id", default=None, help="シートの ID (URL の /d/ の後)")
    p.add_argument("--json-key", default=None, help="サービスアカウントの JSON 鍵")
    p.add_argument("--dry-run", action="store_true", help="違うところを表示するだけ")
    p.add_argument("--force", action="store_true", help="Server の違う行も戻す")
    p.add_argument("--only", default=None, help="戻すラン (例: d001,d005)")
    p.set_defaults(func=cmd_pull)
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
