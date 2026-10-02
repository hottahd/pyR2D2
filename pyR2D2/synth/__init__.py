"""
3D スナップショットからのスペクトル合成（2026-10-02 に ODF-radiation から移した。R2D2plus DEC-589）。

不透明度表（高分解能の ``.h5``）は R2D2plus-input の ``opacity/``（``scripts/build_opacity_table.py``）で作り、ここで読む。
表を読むのに ``h5py`` が要る（``pip install "pyR2D2[synth]"``）。

- :class:`OpacityTable`: 高分解能の不透明度表
- :func:`synth_column` / :func:`synth_ray`: 鉛直・斜めの視線の形式解（ドップラー偏移を含む）
- :func:`ray_contribution`: 寄与関数、:func:`instrument_profile`・:func:`bisector`: 観測との比較
- :func:`check_snapshots`: 3D 出力が実在するか（間引いたランで黙って前の中身を読まないため）

道具は ``python -m pyR2D2.synth.spectrum`` / ``.disk`` / ``.mu_images`` / ``.contribution``。
"""
from .core import (OpacityTable, bisector, doppler_shift_rows, formal_solution, instrument_profile, ray_contribution,
                   ray_geometry, substeps_needed, synth_column, synth_ray)
from .snapshots import check_snapshots

__all__ = ["OpacityTable", "bisector", "check_snapshots", "doppler_shift_rows", "formal_solution", "instrument_profile",
           "ray_contribution", "ray_geometry", "substeps_needed", "synth_column", "synth_ray"]
