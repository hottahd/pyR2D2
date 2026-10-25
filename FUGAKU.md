# 富岳での pyR2D2 導入について (2026-10-05 時点のメモ)

> **2026-10-25 追記**: その後、富岳のログインノードの `~/.venv` (uv の環境) に pyR2D2 を入れて使っている。
> spack の GCC 15 で組み、`libstdc++` を静的に取り込む (41a2083)。手順と罠の正本は R2D2plus の
> `docs/process/environment.md` の富岳の項。以下の「ログを残していない」「未整備」は 10/05 時点の記述で、
> 計算ノード (aarch64) 向けのビルドは今も未整備。

富岳で pyR2D2 をどう入れるか、という問いに対して、現時点でリポジトリから言えることと
言えないことを分けて記録したもの。手順書ではなく、整備されていない箇所の洗い出しが目的。

## 現状: 富岳向けのビルド手順はリポジトリに存在しない

富岳への言及は以下の2箇所だけで、どちらも「ホイールの glibc を富岳のログインノードに
合わせた」という話であり、導入手順ではない。

- `PACKAGING.md:76-78` — Linux ホイールを **manylinux_2_28** (glibc 2.28) にしているのは、
  いちばん古い行き先である富岳のログインノード (RHEL 8) に合わせるため。同じホイールが
  研究室の5台 (glibc 2.39-2.43) でも動く。
- `pyproject.toml` の `[tool.cibuildwheel.linux]` — 同趣旨のコメントと
  `manylinux-x86_64-image = "manylinux_2_28"` / `archs = ["x86_64"]`。

`doc/source/quick_start.rst` や `doc/source/environment.rst` に出てくる富岳
(`SERVER=FUGAKU`、ジョブスクリプト `fg.sh`) は Fortran 本体の話で、pyR2D2 とは無関係。

## ログインノード (x86_64) で使うだけなら、ビルドは不要 (検証済み: 設定ファイルの確認のみ)

`[tool.cibuildwheel.linux]` が `archs = ["x86_64"]` / `manylinux_2_28` なので、タグ push で
GitHub Release に貼られた Linux ホイールがそのまま入る。C++ コンパイラも OpenMP も要らない。

```bash
python -m pip install https://github.com/hottahd/pyR2D2/releases/download/vX.Y.Z/pyr2d2-...-manylinux_2_28_x86_64.whl
```

注意点:

- ホイールは **CPython 3.12 / 3.13 の2枚のみ** (`[tool.cibuildwheel] build = "cp312-* cp313-*"`)。
  pybind11 は stable ABI に対応しないため版ごとに必要で、3.14 は numpy/scipy のホイールが
  揃ってから足す方針。miniforge 側の Python をこの2つのどちらかに合わせる。
- 富岳はアウトバウンドにプロキシ設定が必要。`pip` が外に出られない環境なら、ホイールを
  手元から転送して `python -m pip install ./pyr2d2-*.whl` とするのが確実。

ここまでは `pyproject.toml` / `PACKAGING.md` の記述から確実に言えることだが、
**実際に富岳のログインノードで install して import まで通したログは残っていない。**

## 計算ノード (A64FX = aarch64) は未整備 (未検証: 富岳一般の知識からの推定)

aarch64 のホイールは作っていないので、ソースから `python -m pip install .` になる。
ここで問題になるのが、富岳のログインノードは x86_64、計算ノードは aarch64 という点。

- Fortran 本体は富士通のクロスコンパイラ (`frtpx` / `mpiFCCpx`) でログインノードから
  A64FX 向けバイナリを吐かせているが、setuptools + pybind11 で Python 拡張を
  クロスビルドさせるのは素直にはいかない。
- 現実的には計算ノードを対話確保して (`pjsub --interact`)、そこで aarch64 の
  Python・numpy・scipy を揃えたうえでビルドする形になると思われる。
- ビルド要件は `PACKAGING.md` の通り: **C++20** と、Linux では `setup.py` が `-fopenmp` を
  無条件に付けるので **libgomp 相当のある処理系** (富士通コンパイラの clang モードか
  spack の gcc)。
- `PYR2D2_NATIVE=1` (`-march=native`) は、そのノード専用と割り切るときだけ付ける。

この節は富岳の構成についての一般的な理解に基づく推定であり、リポジトリ内に裏付けはなく、
実際に通した記録もない。

## 残っている判断

1. 用途はログインノードでの解析だけか、計算ノード上でも pyR2D2 を使うのか。
2. 前者だけなら、上のホイール導入手順を一度実機で確認して `PACKAGING.md` に節として足せば済む。
3. 後者が必要なら、計算ノードでのソースビルド手順を確定させる作業が別途必要
   (aarch64 向けホイールを CI で作るか、手順を文書化するだけにするかも含めて)。
