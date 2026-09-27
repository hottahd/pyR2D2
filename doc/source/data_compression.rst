データ圧縮
================================

概要
--------------------------------

R2D2の計算では、1タイムステップごとにMPIランクの数だけバイナリファイルが出力される(``data/remap/qq/*/*/qq.dac.NNNNNNNN.MMMMMMMM`` のような形式)。このままでは解析のたびに大量のファイルを読み込む必要があり、ディスク容量も圧迫する。pyR2D2は、これらを `zarr <https://zarr.readthedocs.io/>`_ ストアにまとめる ``compress()`` によって、

* **データサイズを小さくする**: 2つの方法を組み合わせている。

  1. zarr側でBlosc/zstdによる圧縮をかける
  2. ``pr`` (圧力) や ``te`` (温度)、``op`` (不透明度) のように、保存されている ``ro`` (密度) と ``se`` (エントロピー) からEOS(状態方程式)経由で後から計算できる変数は、デフォルトでは保存しない

* **データアクセスを最適化する**\ (チャンク単位で必要な部分だけを読み込める)

という2つの目的を実現する。

ただし、zarrストア自体は多数の小さいファイルからなるディレクトリなので、そのままだと今度はディレクトリ内のファイル数が非常に多くなる。Gfarm(HPCI共用ストレージなど)のようにファイル数に厳しい制限があるファイルシステムでは、これが新たな問題になる。そこでpyR2D2は、zarrストアを1つの ``.zarr.zip`` ファイルにまとめる機能(``zip()``)も用意している。

このページでは、次の3つの機能を使ったデータ圧縮のワークフローを説明する。

1. ``compress()``: バラバラのバイナリファイルを1つのzarrストアにまとめる(サイズ削減・アクセス最適化)
2. ``check()`` / ``delete()``: zarrストアが元データと一致することを確認してから、元のバイナリファイルを削除する
3. ``zip()`` (``pyR2D2.zarr_util.zip_zarr``): zarrストアを1つの ``.zarr.zip`` ファイルにまとめる(ファイル数制限への対策。後述の通り、目的が異なるので無条件には使わない)

対象となるデータの種類ごとに、 :py:class:`pyR2D2.FullData` (``d.qf``)、:py:class:`pyR2D2.OpticalDepth` (``d.qt``)、:py:class:`pyR2D2.Slice` (``d.qs``) が、それぞれ同じ名前のメソッドを提供している。

.. warning::

    ``delete()`` は元のバイナリファイルを消去する **不可逆な操作** である。必ず ``check()`` (内部で自動的に呼ばれる)がパスすることを確認してから実行すること。``force=True`` を指定するとこの確認をスキップして強制的に削除するので、通常は使用しない方が良い。

compress(): バイナリファイルをzarrにまとめる
------------------------------------------------

:py:class:`pyR2D2.Data` インスタンスを作成した後、3次元フルデータであれば ``d.qf.compress()`` を呼ぶ。

.. code:: python

    import pyR2D2
    d = pyR2D2.Data('../run/d001/data')

    n = 10
    d.qf.compress(n)

デフォルトでは、``data/remap/qq/zarr/qq.NNNNNNNN.zarr`` にzarrストアが作成される。すでに存在する場合は何もせずメッセージを表示するだけなので、上書きしたい場合は ``overwrite=True`` を指定する。

.. code:: python

    d.qf.compress(n, overwrite=True)

デフォルトで保存される変数は ``ro``, ``vx``, ``vy``, ``vz``, ``bx``, ``by``, ``bz``, ``se`` の8つ(``zarr_keys``)である。``pr``, ``te``, ``op`` はデフォルトでは含まれない。これらが必要な場合は ``keys`` 引数で明示的に指定する。

.. warning::

    ``pr``, ``te``, ``op`` を後から ``ro``, ``se`` と :py:class:`pyR2D2.cpp_util.EOS` (あるいは :py:func:`pyR2D2.util.eos_table`) で計算しても、出力ファイルの値とは一致しない。これらは常に表を引くが、R2D2 本体は相対振幅 ``ct = max(|ro1|/ro0, |se1|/se0)`` が ``3e-3`` を下回る点では線形化した EOS を使っている(硬い切り替え)。深部の計算ではほぼ全域が線形側である。本体と同じ値が要るなら :py:func:`pyR2D2.util.eos_switch` を使うこと(単精度の丸めの位置まで同じではないので最終桁は違いうる)。解析で ``pr``, ``te``, ``op`` を使うなら ``keys`` に入れて保存しておくのが確実である。

.. code:: python

    d.qf.compress(n, keys=["ro", "se", "pr", "te", "op"])

計算領域全体ではなく一部分だけを圧縮したい場合は、``i_start`` / ``i_size`` / ``j_start`` / ``j_size`` / ``k_start`` / ``k_size`` で範囲を指定できる。

.. code:: python

    d.qf.compress(n, i_start=0, i_size=64)

データが大きくメモリに乗らない場合は、``lightweight=True`` を指定すると変数を1つずつ読み込んで書き出すので、メモリ使用量を抑えられる。

.. code:: python

    d.qf.compress(n, lightweight=True)

光学的厚さ一定面のデータ(``d.qt``)やスライスデータ(``d.qs``)でも同様に使える。スライスデータは方向(``'x'``, ``'y'``, ``'z'``)の指定が必要である。

.. code:: python

    d.qt.compress(n)
    d.qs.compress(n, direc='x')

check() と delete(): 元データの安全な削除
------------------------------------------------

``compress()`` でzarrストアを作成しても、元のバイナリファイルは自動的には消えない。作成したzarrストアが元のバイナリファイルと一致するかを確認するのが ``check()`` である。

.. code:: python

    d.qf.check(n)
    # -> "Check passed at n=10" のように表示され、True が返る

``check()`` は、zarrファイルの存在確認・変数の存在確認をした上で、実際にzarr経由とバイナリ経由の両方でデータを読み込み、数値が一致するかどうかを比較する。一致しない場合は ``False`` を返す。

``delete()`` は、内部で ``check()`` を呼び、パスした場合のみ元のバイナリファイルを削除する。

.. code:: python

    d.qf.delete(n)

``check()`` を毎回実行したくない場合(すでに確認済みの場合など)は、``force=True`` を指定すると確認をスキップしてそのまま削除する。前述の通り、事故を防ぐためこのオプションは慎重に使うこと。

.. code:: python

    d.qf.delete(n, force=True)  # checkをスキップして削除。取り扱い注意

zip(): ファイル数制限への対策
------------------------------------------------

``compress()`` によってデータサイズ・アクセス速度の問題は解決するが、zarrストアの実体は多数の小さいファイルからなるディレクトリである。Gfarmなどファイル数に制限のあるファイルシステムでは、この点が新たな問題になる。``zip()`` はこれを1つの ``.zarr.zip`` ファイルにまとめることで、ファイル数の問題を解消する。

.. code:: python

    d.qf.zip(n)
    # -> data/remap/qq/zarr/qq.NNNNNNNN.zarr.zip が作成される

zip化の際は圧縮を行わない(``ZIP_STORED``)。zarr側ですでにBlosc/zstdによる圧縮がかかっているため、二重に圧縮する必要がないからである。

元のzarrディレクトリを消して良い場合は ``remove_original=True`` を指定する。この場合、zip化した内容が元のディレクトリと(配列の形・dtypeについて)一致することを確認してから削除するので、``delete()`` と同様に安全に使える。

.. code:: python

    d.qf.zip(n, remove_original=True)

``.zarr.zip`` になったデータは、読み込み時に特別な指定をしなくても自動的に認識される。``zarr_flag=True`` で読み込む際、``.zarr`` ディレクトリと ``.zarr.zip`` ファイルのどちらが存在してもそのまま読み込める。

.. code:: python

    d.qf.read(n=n, zarr_flag=True)

低レベルAPIとして、任意の ``.zarr`` ディレクトリに対して直接 ``pyR2D2.zarr_util.zip_zarr`` を呼ぶこともできる。

.. code:: python

    import pyR2D2
    pyR2D2.zarr_util.zip_zarr('path/to/data.zarr', remove_original=True)

.. important::

    **zipすべきタイミングについて**

    zarrストアはディレクトリのまま(zip化しない)であれば、複数プロセスからの並列アクセスが可能である(ただしpyR2D2側でこれを積極的に活用する実装は現時点では用意されていない)。一方、``.zarr.zip`` にまとめると、そうした並列アクセスの余地は失われる。

    したがって、``zip()`` は無条件に実行するのではなく、**Gfarm(HPCI共用ストレージなど)のようにファイル数制限が厳しいシステムに保存する場合に限定して実行する**\ のが良い。ローカルディスクや通常のPOSIXファイルシステムに置く場合は、zarrストアをディレクトリのまま(zip化せず)残しておく方が、将来的な並列アクセスの余地を残せる。

time.zip: 時刻ファイルの圧縮
------------------------------------------------

各タイムステップの出力時刻を記録した ``data/time/`` 以下のファイル群も、数が多くなりがちである。:py:meth:`pyR2D2.Data.zip_time` を使うと、これらを ``data/time.zip`` にまとめられる。この操作もファイル数対策としての位置づけであり、``zip()`` と同様の考え方(Gfarmなどファイル数制限の厳しい環境で使う)が当てはまる。

.. code:: python

    d.zip_time()

すでに ``time.zip`` に含まれているファイルはスキップされ、新しく追加されたファイルだけが追記される。元ファイルを消してよければ ``remove_original=True`` を指定する。

.. code:: python

    d.zip_time(remove_original=True)

``time_read()`` は、デフォルトでは元ファイルを優先して読み込み、存在しない場合に ``time.zip`` から読み込む。``use_zip=True`` を指定すると、常に ``time.zip`` から読み込む。

.. code:: python

    d.time_read(n, use_zip=True)

典型的なワークフロー
------------------------------------------------

計算が完了したタイムステップについて、まとめて圧縮・確認・削除を行う例を示す。ここでは、通常のディスクに保存する場合(``zip()`` は行わない)と、Gfarmのようなファイル数制限の厳しいストレージに保存する場合(``zip()`` まで行う)の2通りを示す。

**通常のディスクに保存する場合**

.. code:: python

    import pyR2D2

    d = pyR2D2.Data('../run/d001/data')

    for n in range(d.p.nd + 1):
        d.qf.compress(n)
        if d.qf.check(n):
            d.qf.delete(n)
        else:
            print(f"n={n}: check failed, keeping original binary files")

**Gfarmなど、ファイル数制限の厳しいストレージに保存する場合**

.. code:: python

    import pyR2D2

    d = pyR2D2.Data('../run/d001/data')

    for n in range(d.p.nd + 1):
        d.qf.compress(n)
        if d.qf.check(n):
            d.qf.delete(n)
            d.qf.zip(n, remove_original=True)
        else:
            print(f"n={n}: check failed, keeping original binary files")

    d.zip_time(remove_original=True)

.. note::

    計算がまだ進行中のディレクトリに対してこの処理を行う場合、計算プロセスが書き込み中のファイルを誤って対象にしないよう注意すること。基本的には、計算が完全に終了した(あるいは十分に古い)タイムステップに対してのみ実行するのが安全である。

R2D2plus の新出力形式 (``format = "compressed"``)
------------------------------------------------------

R2D2plus (C++版) は ``[output] format = "compressed"`` で、計算の時点で圧縮した出力を書ける(R2D2plus DEC-576)。上の ``compress()`` とは別物で、pyR2D2 側の操作は要らない。``param/format.toml`` があり ``format = "compressed"`` と書かれていれば :py:class:`pyR2D2.Data` が自動で見分け、``d.qf.read(n)``, ``d.qx``, ``d.qz``, ``d.qm``, ``d.qr``, ``d.qt``, ``d.qs``, ``d.q2``, ``d.qp``/``d.qa`` をこれまでと同じ引数・同じ属性名・同じ配列の形で使える(``d.vc`` は従来形式と同じファイル)。従来形式(Fortran R2D2 と R2D2plus の ``format = "legacy"``)の読み方は変わらない。

* 3D の ``remap/qq`` は 1 ランク 1 ファイルの ``r2d2plus-z`` コンテナ(変数ごとに byte shuffle + zstd)。中身は従来形式とビット単位で同一。読むには zstd の復号器が要る(``pip install numcodecs``。``pyR2D2[zarr]`` を入れていれば入っている)。
* ``pr``, ``te``, ``op`` は書かれない。pyR2D2 は要求されたとき(``keys="all"`` の既定も含む)に、出力に同梱された ``param/eos_table.tbl`` を使って ``ro``, ``se`` から計算する。式は R2D2plus の従来形式の書き出しと同じで、単精度の丸めの位置まで写してあるので、3D の ``remap/qq`` と tau の ``pr``, ``te``, ``op`` は従来形式のファイルとビット単位で一致する(対照ランで確認)。slice と 2D の ``pr``, ``te`` は R2D2plus が倍精度の状態から作っていたので、単精度の出力から計算するこちらとは最終桁が違いうる(対照ランの slice で背景量に対し 2e-8 以内)。
* ``ph`` (発散クリーニングの変数) は既定では書かれない。``keys="all"`` では属性が ``None`` になり、明示的に要求すると例外になる。
* tau の高さは ``r - rstar`` [cm] で書かれる。``d.qt.height`` (と ``height01``, ``height001``) でそのまま読める。従来の ``d.qt.he`` (中心からの半径) は ``rstar + height`` を倍精度で返す(従来の単精度の半径は約 82 m 刻みだった)。
* ``d.qp``, ``d.qa`` の配列は単精度(従来形式は倍精度)。値は同じ。

``r2d2plus-z`` のファイル 1 個や、R2D2plus のリスタートを直接読む関数もある(リスタートは実験的)。

.. code:: python

    from pyR2D2.data_io import compressed

    header, qq = compressed.read_r2d2plus_z(
        "data/remap/qq/00000/00000000/qq.z.00000001.00000000")
    state, meta = compressed.read_restart_state("data/restart/e", rank=0)

前の出力を R2D2plus の初期条件にする
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`pyR2D2.write_initial_state` は R2D2plus の初期状態ファイル(``r2d2plus-z``、``kind = "initial_state"``)を書く。R2D2plus 側では ``[initial_condition] type = "file"``、``[initial_condition.params] path = "..."`` で読む。変数名と軸は pyR2D2 の 3D 配列と同じ(``vx``, ``bx`` が鉛直成分、形は ``(ix, jx, kx)``、領域全体の物理セルのみ)なので、``d.qf.read(n)`` で読んだ配列をそのまま渡せばよい(従来形式・新形式のどちらのランでも同じ)。

.. warning::

    ``ro`` と ``se`` は **背景からの摂動** (``d.qf.ro``, ``d.qf.se`` と同じ量) を渡す。全量 (``ro0 + ro``) を渡してはならない。

.. code:: python

    import pyR2D2

    d = pyR2D2.Data('../run/d001/data')
    n = 120
    keys = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]
    d.qf.read(n, keys=keys)
    pyR2D2.write_initial_state(
        "initial_state.z",
        {k: d.qf.__dict__[k] for k in keys},
        attributes={"source": str(d.datadir), "nd": n,
                    "time": float(d.time_read(n, verbose=False))},
    )

書き出しには ``xxhash`` が要る(``pip install xxhash``、または ``pyR2D2[compressed]``)。R2D2plus は変数ごとの XXH3-64 を照合して読むので、チェックサムを省略できない。

境界値の時系列を R2D2plus の境界条件に与える
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`pyR2D2.write_boundary_series` は R2D2plus の汎用境界条件の規則 ``"file"`` が読む時系列(``kind = "boundary_series"``)を書く。``values`` の形は ``(n1, n2, nt)`` (``n1`` = 第1水平 ``jx``、``n2`` = 第2水平 ``kx``、**時刻が最後の軸**)、``times`` は狭義単調増加。値はその壁の ghost 層すべてに使われ、ステップ開始時刻で線形内挿される(範囲外は端の値)。

.. warning::

    ``[boundary.params]`` の変数名は **C++ の成分名** (``ro, vx, vy, vz, bx, by, bz, se, ps``、C++ では **z が鉛直**) である。pyR2D2 (Fortran の名前) の ``vx`` (鉛直) は ``vz``、``vy`` は ``vx``、``vz`` は ``vy`` になる(磁場も同様、``ph`` は ``ps``)。``ro``, ``se`` は背景からの摂動を与える。

.. code:: python

    import numpy as np
    import pyR2D2

    d = pyR2D2.Data('../run/d001/data')
    steps = range(100, 121)
    planes, times = [], []
    for n in steps:
        d.qf.read(n, keys=["vx"])
        planes.append(d.qf.vx[0].copy())   # 下端の面の鉛直速度 (jx, kx)
        times.append(d.time_read(n, verbose=False))
    pyR2D2.write_boundary_series("bottom_vz.z", times, np.stack(planes, axis=-1))
    # R2D2plus 側: [boundary.params] bottom.vz = "file", bottom.vz_file = "bottom_vz.z"

最終更新日：|today|
