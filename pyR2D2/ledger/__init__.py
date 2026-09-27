"""ランの一覧 (台帳) を Google スプレッドシートへ送る道具 ``r2d2plus-ledger``。

:py:class:`pyR2D2.Data` を使わずに、各ランの小さなファイル (``params.dac``、
``run_summary.toml``、``origin.toml``、``effective_config.txt`` など) だけを読む。
使い方は :py:mod:`pyR2D2.ledger.cli`、列と単位の決め方は
:py:mod:`pyR2D2.ledger.sheet` と :py:mod:`pyR2D2.ledger.units` を参照。
"""

import hashlib
import os
import tomllib
from pathlib import Path

__all__ = ["unsent_hint"]


def unsent_hint(datadir):
    """``pyR2D2.Data`` で開いたランに台帳へ未送信の変更があれば 1 行の案内を返す。

    ``param/run_summary.toml`` (R2D2plus) が無いラン (Fortran 版) では ``None``。
    ``param/ledger_sent.toml`` に残した ``run_summary.toml`` の SHA-256 と比べるだけ
    なので軽い。どんな失敗でも例外を出さず ``None`` を返す。
    """
    try:
        datadir = Path(datadir)
        summary = datadir / "param" / "run_summary.toml"
        if not summary.is_file():
            return None
        sha = hashlib.sha256(summary.read_bytes()).hexdigest()
        sent_sha = None
        sent = datadir / "param" / "ledger_sent.toml"
        if sent.is_file():
            with open(sent, "rb") as f:
                sent_sha = tomllib.load(f).get("sent", {}).get("run_summary_sha256")
        if sent_sha == sha:
            return None
        run_parent = Path(os.path.abspath(os.fspath(datadir))).parent.parent
        return f"台帳へ未送信の変更があります: r2d2plus-ledger push {run_parent}"
    except Exception:
        return None
