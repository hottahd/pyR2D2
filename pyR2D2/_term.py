"""コマンドの出力の色 (r2d2plus-ledger・pack・archive・restore で共通)。

configuration の ``claude/claude-sessions`` と同じ決まり:

* 色を付けるのは出し先が端末のときだけ (``isatty``)。``NO_COLOR`` が設定されていれば付けない。
  パイプ・ファイル・cron・tmux run-shell には制御文字を出さない。
* 緑 = 正常・済み、赤 = 失敗・止まっている、黄 = 注意、薄い色 = ID などの補足、太字 = 見出し。
* 表の幅は色を付ける前の文字列で揃える (制御文字を幅に数えない)。
"""

import os
import sys

_CODES = {"green": "32", "red": "31", "yellow": "33", "dim": "2", "bold": "1"}


def enabled(stream=None):
    """``stream`` (既定は標準出力) に色を付けてよいか。呼ぶたびに調べる。"""
    if "NO_COLOR" in os.environ:
        return False
    stream = sys.stdout if stream is None else stream
    try:
        return stream.isatty()
    except (AttributeError, ValueError):
        return False


def paint(text, style, stream=None):
    """``text`` を ``style`` (green・red・yellow・dim・bold) の色にする。色を付けないときはそのまま。"""
    text = str(text)
    if not style or not text or not enabled(stream):
        return text
    return f"\033[{_CODES[style]}m{text}\033[0m"


def green(text, stream=None):
    return paint(text, "green", stream)


def red(text, stream=None):
    return paint(text, "red", stream)


def yellow(text, stream=None):
    return paint(text, "yellow", stream)


def dim(text, stream=None):
    return paint(text, "dim", stream)


def bold(text, stream=None):
    return paint(text, "bold", stream)


def warning(message, stream=None):
    """``warning: …`` の 1 行 (``warning:`` を黄色に)。"""
    return f"{yellow('warning:', stream)} {message}"


def error(prog, message, stream=None):
    """``<prog>: error: …`` の 1 行 (標準エラー向け、``error:`` を赤に)。"""
    stream = sys.stderr if stream is None else stream
    return f"{prog}: {red('error:', stream)} {message}"
