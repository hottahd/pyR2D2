"""ランのメモ ``data/param/run_note.toml`` (R2D2plus DEC-600) の書き換えと、シートの値からの読み戻し。

メモは人が書く。``r2d2plus-ledger pull`` はシートで書き換えた目的・プロジェクト・保存・結論を、
1 か所ずつ利用者に確かめてからメモへ戻す (2026-10-28、利用者の指示)。

書き換えは R2D2plus ``tools/run_note.py`` の ``update`` と同じく、該当する行だけを変え、
ほかのバイト (見出し・コメント・並び・末尾の改行) は保つ。

    [note]
    purpose   = "..."            # 目的
    project   = "..."            # プロジェクト
    retention = "keep"           # 保存 (keep | regen_only | auto。thin は古いメモだけ。R2D2plus DEC-666)
    expires   = "2026-12-31"     # 保存の「(… まで残す)」
    cited_in  = ["DEC-533"]      # 結論
"""

import re
import tomllib
from pathlib import Path

NOTE_FILE = "run_note.toml"

# シートの「保存」の表示 (collect.RETENTION_LABELS と同じ)。R2D2plus DEC-666 (2026-10-31) で意味を改めた:
# auto は未定で消さない、thin は新しく使わず keep と同じに扱う
RETENTION_LABELS = {
    "keep": "残す",
    "thin": "間引き済み (残す)",
    "regen_only": "再生成の情報だけ",
    "auto": "未定 (消さない)",
}
# 2026-10-31 までの表示。シートに残っていても pull で読めるようにする
OLD_RETENTION_LABELS = {
    "全部残す": "keep",
    "間引いて残す": "thin",
    "自動 (DEC-599)": "auto",
}

# シートから戻してよい列 (列名 → メモの項目)。ほかの列は機械が書く値なので戻さない。
PULL_COLUMNS = ("プロジェクト", "目的", "保存", "結論")

_EXPIRES_RE = re.compile(r"^(.*?)\s*\((\d{4}-\d{2}-\d{2})\s*まで残す\)\s*$")
_DEC_RE = re.compile(r"^DEC-\d+$")


def note_path(run_dir):
    return Path(run_dir) / "data" / "param" / NOTE_FILE


def read(run_dir):
    p = note_path(run_dir)
    if not p.exists():
        return None
    return tomllib.loads(p.read_text(encoding="utf-8")).get("note", {})


def parse_cell(column, text):
    """シートのセルをメモの項目へ直す。

    Returns
    -------
    fields : dict
        書き換える項目 (``{"retention": "thin", "expires": "2026-12-31"}`` など)
    Raises
    ------
    ValueError
        読めない値 (保存の表示が既知のどれでもない、結論が DEC-NNN でない)
    """
    text = text.strip()
    if column == "目的":
        return {"purpose": text}
    if column == "プロジェクト":
        return {"project": text}
    if column == "保存":
        expires = None
        m = _EXPIRES_RE.match(text)
        if m:
            text, expires = m.group(1).strip(), m.group(2)
        keys = {label: key for key, label in RETENTION_LABELS.items()}
        keys.update(OLD_RETENTION_LABELS)
        keys.update({key: key for key in RETENTION_LABELS})
        if text not in keys:
            raise ValueError(
                f"保存の値 {text!r} が読めない ("
                + "・".join(RETENTION_LABELS.values())
                + " のどれか、後ろに「(YYYY-MM-DD まで残す)」を付けられる)"
            )
        fields = {"retention": keys[text]}
        if expires:
            fields["expires"] = expires
        return fields
    if column == "結論":
        decs = [s for s in re.split(r"[,、\s]+", text) if s]
        bad = [s for s in decs if not _DEC_RE.match(s)]
        if bad:
            raise ValueError(f"結論は DEC-NNN を並べる ({', '.join(bad)} が読めない)")
        return {"cited_in": decs}
    raise KeyError(column)


def _quote(s):
    return '"' + str(s).replace("\\", "\\\\").replace('"', '\\"') + '"'


def _toml_value(value):
    if isinstance(value, list):
        return "[" + ", ".join(_quote(c) for c in value) + "]"
    return _quote(value)


def update(run_dir, fields):
    """メモの ``fields`` の行だけを書き換える (無い項目は ``[note]`` の最後に足す)。ほかのバイトは保つ。

    メモが無ければ新しく作る。書いた値が読み戻せることを確かめてから置き換える。
    """
    p = note_path(run_dir)
    if p.exists():
        lines = p.read_text(encoding="utf-8").splitlines(keepends=True)
    else:
        lines = ["# ランのメモ（人が書く。DEC-600）。機械が書く状態は run_summary.toml。\n"]
    section, header, note_end = None, None, None
    found = {}
    for i, line in enumerate(lines):
        m = re.match(r"\s*\[([^\]]+)\]\s*(#.*)?$", line)
        if m:
            section = m.group(1).strip()
            if section == "note":
                header = i
            continue
        if section == "note":
            note_end = i
            m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_-]*)\s*=", line)
            if m and m.group(1) in fields:
                found[m.group(1)] = i
    if header is None:  # [note] が無い: 最後に足す
        if lines and not lines[-1].endswith("\n"):
            lines[-1] += "\n"
        lines.append("[note]\n")
        header = len(lines) - 1
    if note_end is None:
        note_end = header
    for key, i in found.items():
        eol = "\n" if lines[i].endswith("\n") else ""
        lines[i] = f"{key} = {_toml_value(fields[key])}" + eol
    extra = [f"{k} = {_toml_value(v)}\n" for k, v in fields.items() if k not in found]
    if extra:
        if not lines[note_end].endswith("\n"):
            lines[note_end] += "\n"
        lines[note_end + 1 : note_end + 1] = extra
    new_text = "".join(lines)
    new_note = tomllib.loads(new_text).get("note", {})
    for key, value in fields.items():
        if new_note.get(key) != value:
            raise RuntimeError(f"{p}: {key} を書けなかった")
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".toml.tmp")
    tmp.write_text(new_text, encoding="utf-8")
    tmp.replace(p)
    return p
