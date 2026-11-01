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

from .. import _term as term
from . import collect, note, sheet
from .units import Quantity, render, choose_unit

PROG = "r2d2plus-ledger"

# 画面に出す名前 (シートの列名・値は日本語のまま。R2D2plus DEC-666)
# status の「状態」の値 → R2D2plus の状態名と色
STATUS_NAMES = {label: key for key, label in collect.STATUS_LABELS.items()}
STATUS_NAMES[collect.STALE_LABEL] = "unknown"
STATUS_STYLES = {"running": "green", "normal_end": "green", "wall_limit": "yellow",
                 "error": "red", "unknown": "yellow"}


def _records(path, server, du, only=None):
    runs = collect.find_runs(path)
    if only:
        wanted = set(only)
        runs = [r for r in runs if r.name in wanted]
        missing = wanted - {r.name for r in runs}
        if missing:
            print(term.warning("runs not found: " + ", ".join(sorted(missing)), sys.stderr), file=sys.stderr)
    return [collect.collect_run(r, server=server, du=du) for r in runs]


def _fmt_time(q):
    if not isinstance(q, Quantity):
        return ""
    unit = choose_unit(q.family, list(q.values))
    return render(q, unit) + " " + unit


def cmd_status(args, out=print):
    """各ランの要約と、台帳へ未送信の変更があるかを表示する。シートには触らない。"""
    records = _records(args.path, args.server, du=False)
    header = ("CASE ID", "STATUS", "TIME", "STEP", "UNSENT")
    rows, styles = [], []
    for r in records:
        unsent = r.unsent()
        mark = "never sent" if unsent == "never" else ("yes" if unsent else "-")
        label = r.values.get("状態") or ""
        status = STATUS_NAMES.get(label, label)
        rows.append((r.caseid, status, _fmt_time(r.values.get("到達時刻")), r.values.get("step") or "", mark))
        styles.append((None, STATUS_STYLES.get(status), None, None, "dim" if mark == "-" else "yellow"))
    widths = [max(_width(row[i]) for row in [header] + rows) for i in range(len(header))]
    # 色の制御文字は幅に数えないよう、幅を揃えてから色を付ける
    widths[-1] = 0  # 最後の列は揃えない (行末に空白を残さない)
    out(term.bold("  ".join(_pad(c, w) for c, w in zip(header, widths)).rstrip()))
    for row, style in zip(rows, styles):
        cells = [term.paint(_pad(c, w), s) for c, w, s in zip(row, widths, style)]
        out("  ".join(cells).rstrip())
    return records


def _n(count, noun):
    return f"{count} {noun}" + ("" if count == 1 else "s")


def _width(s):
    # 全角は 2 桁で数える (端末で列を揃えるため)
    return sum(2 if ord(ch) > 0x2E7F else 1 for ch in str(s))


def _pad(s, w):
    return str(s) + " " * max(0, w - _width(s))


def cmd_export(args, out=print):
    """管理列を CSV に書く (単位は見出しに入る)。"""
    records = _records(args.path, args.server, du=not args.no_du)
    headers, rows = sheet.render_records(records)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerows(rows)
    out(term.green(f"wrote {_n(len(records), 'run')} to {args.output}"))
    return headers, rows


def open_worksheet(sheet_name=None, sheet_id=None, json_key=None):
    """Google スプレッドシートの最初のワークシートを開く (ここでだけ gspread を読む)。"""
    try:
        import gspread  # noqa: F401
        from google.oauth2.service_account import Credentials  # noqa: F401
    except ImportError as e:
        raise SystemExit(
            term.error(PROG, "push and pull need gspread and google-auth "
                       f"(pip install gspread google-auth): {e}")
        )
    from ..write.google.google import init_gspread

    if json_key is None:
        keys = sorted(glob.glob(os.path.join(os.environ["HOME"], "json", "*")))
        if not keys:
            raise SystemExit(term.error(PROG, "no service account JSON key in ~/json/ (give one with --json-key)"))
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

    _warn(plan, out)
    # 管理列を COLUMNS の順に並べ替える (人の列は動かさない)。書いた後の見出しで決める
    moves, final = sheet.column_moves(plan.header_after) if reorder else ([], [])

    if dry_run:
        if moves:
            out("would reorder columns: " + ", ".join(h for h in final if h))
        if worksheet is None:
            out(term.dim(f"(dry-run: not connected to sheet {label}; plan against an empty sheet)"))
        for r, c, text in sorted(plan.cells):
            out(f"{term.dim(sheet.a1(r, c))}\t{text}")
        out(f"(dry-run) would write {_n(len(plan.cells), 'cell')} for {_n(len(plan.written), 'run')}")
        return plan

    sheet.apply_plan(worksheet, plan)
    if moves:
        sheet.apply_column_moves(worksheet, moves)
        out(f"reordered the managed columns ({_n(len(moves), 'move')})")
    for rec in plan.written:
        collect.write_sent(rec, label)
    out(term.green(f"sent {_n(len(plan.written), 'run')} to sheet {label} ({_n(len(plan.cells), 'cell')})"))
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
            raise SystemExit(term.error(PROG, "pull asks y/n for each change; run it from a terminal "
                                        "(or use --dry-run to only show the changes)"))
        ask = input
    if worksheet is None:
        worksheet = open_worksheet(
            sheet_name=None if sheet_id else label, sheet_id=sheet_id, json_key=json_key
        )
    plan = sheet.plan_pull(worksheet.get_all_values(), records, server, force=force)
    _warn(plan, out)
    if not plan.changes:
        out(term.green(f"sheet {label} has no purpose, project, retention or cited_in that differs from the local notes"))
        return []

    applied = []
    for i, c in enumerate(plan.changes, start=1):
        out(term.bold(f"[{i}/{len(plan.changes)}] {c.record.caseid} {note.PULL_KEYS.get(c.column, c.column)}"))
        out(f"  local: {c.local or term.dim('(empty)')}")
        out(f"  sheet: {c.remote}")
        if c.column == "保存":
            cur = note.read(c.record.run_dir) or {}
            if "expires" not in c.fields and cur.get("expires"):
                out(term.dim(f"  (local expires = {cur['expires']} is not in the sheet; it is kept. "
                             "Remove it in the local note to drop it)"))
            if "keep_parts" not in c.fields and cur.get("keep_parts"):
                parts = ", ".join(cur["keep_parts"])
                out(term.dim(f"  (local keep_parts = [{parts}] is not in the sheet; it is kept. "
                             "Remove it in the local note to drop it)"))
        if dry_run:
            continue
        while True:
            ans = ask("  apply to the local note? [y/n/q]: ").strip().lower()
            if ans in ("y", "n", "q"):
                break
        if ans == "q":
            break
        if ans == "y":
            p = note.update(c.record.run_dir, c.fields)
            applied.append(c)
            out(term.green(f"  applied: {p}"))
    if dry_run:
        out(f"(dry-run) {_n(len(plan.changes), 'change')} from the local notes; nothing was written")
    else:
        out(term.green(f"applied {len(applied)} of {len(plan.changes)} changes to the local notes"))
    return applied


def _warn(plan, out):
    for caseid, reason in plan.skipped:
        out(term.warning(f"skipped {caseid}: {reason}"))
    for w in plan.warnings:
        out(term.warning(w))


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
        description="Send the list of runs (the ledger) to a Google spreadsheet",
    )
    sub = ap.add_subparsers(dest="command", required=True)

    def common(p):
        p.add_argument("path", help="<project>/run (or the directory of one run)")
        p.add_argument("--server", default=None, help="server name (default: the name used last time, else the host name)")

    p = sub.add_parser("status", help="show each run and whether it has unsent changes (sends nothing)")
    common(p)
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("export", help="write the managed columns to a CSV file")
    common(p)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--no-du", action="store_true", help="do not measure the size (disk usage of data/)")
    p.set_defaults(func=cmd_export)

    p = sub.add_parser("push", help="update the managed columns of the sheet")
    common(p)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--sheet", default=None, help="sheet name (default: the project name)")
    g.add_argument("--sheet-id", default=None, help="sheet ID (the part after /d/ in the URL)")
    p.add_argument("--json-key", default=None, help="JSON key of the service account")
    p.add_argument("--dry-run", action="store_true", help="only show the cells that would be written")
    p.add_argument("--force", action="store_true", help="also overwrite rows whose Server differs")
    p.add_argument("--only", default=None, help="runs to send (e.g. d001,d005)")
    p.add_argument("--no-du", action="store_true", help="do not measure the size (disk usage of data/)")
    p.add_argument(
        "--keep-column-order",
        action="store_true",
        help="do not reorder the managed columns into the fixed order (Case ID, project, purpose, retention, size, ...)",
    )
    p.set_defaults(func=cmd_push)

    p = sub.add_parser(
        "pull", help="bring purpose, project, retention and cited_in edited in the sheet back to the local notes, asking y/n for each"
    )
    common(p)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--sheet", default=None, help="sheet name (default: the project name)")
    g.add_argument("--sheet-id", default=None, help="sheet ID (the part after /d/ in the URL)")
    p.add_argument("--json-key", default=None, help="JSON key of the service account")
    p.add_argument("--dry-run", action="store_true", help="only show the differences")
    p.add_argument("--force", action="store_true", help="also bring back rows whose Server differs")
    p.add_argument("--only", default=None, help="runs to bring back (e.g. d001,d005)")
    p.set_defaults(func=cmd_pull)
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
