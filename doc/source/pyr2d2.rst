=================
pyR2D2
=================

Introduction
----------------------

Python module for dealing with calculation results from R2D2 code is available.

R2D2で実施した計算を取り扱うためのPythonモジュールを公開している。

https://github.com/hottahd/pyR2D2

You can install it by cloning the repository and running pip install as

以下のようにしてインストールできる。データの読み込みには :code:`zarr` パッケージが必要なので、通常は :code:`[zarr]` を付けてインストールする。

.. code-block:: bash

    git clone https://github.com/hottahd/pyR2D2
    cd pyR2D2
    pip install .[zarr]

matplotlibを使ったプロット機能もあわせて使いたい場合は :code:`[all]` を指定する。

.. code-block:: bash

    pip install .[all]

Here is a simple example of how to use the pyR2D2 module.

PythonでpyR2D2で定義された関数を使うには以下のようにする。

.. code:: python

    import pyR2D2


All APIs are shown in :doc:`/api_reference`.

すべてのAPIは :doc:`/api_reference` に示されている。

:class:`pyR2D2.Data`
--------------------------

pyR2D2D2 provides a class :py:class:`pyR2D2.Data` to handle data in an object-oriented manner. Instance can be created as follows:

pyR2D2には :py:class:`pyR2D2.Data` クラスが定義してあり、これをオブジェクト指向的に用いてデータを取り扱う。以下のようにしてインスタンスを作成する。

.. code:: python

    import pyR2D2
    datadir = '../run/d002/data'
    d = pyR2D2.Data(datadir)

.. code:: python

    help(pyR2D2)
    help(pyR2D2.Data)

などとすると実行環境で、モジュール全体や各関数の簡単な説明を見ることができる。            

最終更新日：|today|
出力の番号と時刻
--------------------------

ランの出力は 4 つの族に分かれ、それぞれ番号と時刻を持つ。

.. list-table::
   :header-rows: 1

   * - 族
     - 読み手
     - 番号
     - 間隔
     - 時刻
   * - 3D (``remap/qq``)
     - ``d.qf`` など
     - ``d.nd``
     - ``d.dtout``
     - ``time/mhd``
   * - 平均量など (``remap/vl``)
     - ``d.vc``
     - ``d.nd_vl``
     - ``d.dtout_vl``
     - ``time/vl``
   * - 光学的深さ (``tau``)
     - ``d.qt``
     - ``d.nd_tau``
     - ``d.dtout_tau``
     - ``time/tau``
   * - 断面 (``slice``)
     - ``d.qs``
     - ``d.nd_slice``
     - ``d.dtout_slice``
     - ``time/slice``

R2D2plus の ``output_cadence = 2`` のラン (2026-10-31 以降に始めたもの) では、4 つが別々の間隔と番号を持つ。
それより前のラン (``d.output_cadence == 1``。延長しても変わらない) では、vl は 3D と、slice は tau と同じ番号・
時刻なので、``d.nd_vl == d.nd``、``d.nd_slice == d.nd_tau`` になる。どちらのランでも、上の名前で書けば同じように動く。

時刻は ``dtout * n`` で求めず、``time/`` から読む (途中で間隔を変えたランもある)。

.. code:: python

    t = d.time_read(n, kind="vl")          # remap/vl の番号 n の時刻 (kind は "mhd"・"vl"・"tau"・"slice")
    steps, t = d.times("slice")            # slice の全部の番号と時刻
    m = steps[abs(t - t0).argmin()]        # 時刻 t0 にいちばん近い slice の番号

今までの数え方の R2D2plus のランでは、3D を ``three_d_every`` 回に 1 回しか書かないことがある
(``time/mhd`` と ``remap/vl`` は毎回書く)。``d.qq_steps()`` は 3D のある番号の一覧を返す
(``param/runs/*/run.toml`` の ``three_d_every`` と、手元にある 3D のファイルから求める)。
``Data(..., verbose=True)`` の表示も、4 つの族のランでは族ごとに 1 行、``three_d_every`` を使ったランでは 3D のある番号の数を出す。

EUV・軟 X 線の強度 (emission measure)
---------------------------------------

R2D2plus の量の並び 2 のラン (restart の ``meta.toml`` に ``vl_layout = 2``。2026-10-31 以降に初期条件・``r2d2plus_run new``・
``r2d2plus_upgrade`` から始めた、出力の数え方 2 のラン。R2D2plus DEC-671) では、``remap/vl`` に次のものが増える。
``d.vc.read(n)`` が一緒に読み、ファイルの無いランでは ``None`` になる。

.. list-table::
   :header-rows: 1

   * - 属性
     - 形
     - 中身
   * - ``bxcxco`` 〜 ``bzczco``、``qrk``・``qrte``・``qart``・``qcond``・``qthin``・``qoth``、``keart``・``meart``
     - (ix, jx)
     - B と電流の相関、項目別の加熱率 [erg cm^-3 s^-1]、人工拡散が変えた運動・磁気エネルギー (``vl_xy`` の末尾)
   * - ``fc``
     - (ix + 1, jx)
     - 熱伝導の鉛直流束 [erg cm^-2 s^-1] (``vl_flux`` の末尾)
   * - ``hist_te``・``hist_ro``
     - (ix, 100)・(ix, 80)
     - 高さごとの log10 T・log10(ρ/ρ0) の面積の割合 (各高さで和 1)。箱は ``d.vc.hist_logte_edges``・``hist_logro_edges``
   * - ``em_x``・``em_y``・``em_z``
     - (jx, kx, 60)・(ix, kx, 60)・(ix, jx, 60)
     - 視線ごとの EM = Σ n_e² dl [cm^-5] を log10 T の箱 ([5, 8)、0.05 dex) に分けたもの。最後の文字が視線の向き
       (``em_x`` は鉛直の視線で上から見た図)。Cartesian で熱伝導があり EOS 表に ``ne`` があるランだけ

AIA・XRT の強度は、EM に応答関数を掛けて作る (:py:mod:`pyR2D2.synth.euv`)。

.. code:: python

    from pyR2D2.synth import euv

    d.vc.read(n)                                   # n は vl の番号 (d.nd_vl)
    aia = euv.intensity(d.vc.em_x, "aia")          # {"A94": (jx, kx), "A131": ..., ...} [DN s^-1 pix^-1]
    xrt = euv.intensity(d.vc.em_x, "xrt", channels=["Al-poly", "Be-thin"])
    dem = d.vc.em_x / np.diff(10 ** d.vc.dem_logte_edges)   # DEM(T) [cm^-5 K^-1]

同梱の応答関数 (``pyR2D2/synth/data/``) は SSW の公開ファイルから作ったもの。表の頭に出どころが書いてある
(``euv.load_response("aia").meta``)。

- ``"aia"``: ``aia_get_response(/temp, /dn, /evenorm, /chiantifix)`` (V10、時間変化の補正なし) と同じ手順。
  94・131・171・193・211・335 Å (304 Å は光学的に厚いので入れない)。
- ``"xrt"``: ``make_xrt_temp_resp`` (CHIANTI 10、コロナの組成) と同じ手順。CCD とフィルターの汚染の厚さは表を作った日付のもの。

応答関数は差し替えられる。``intensity(em, "my_table.txt")`` のようにファイルを渡すか、:py:class:`pyR2D2.synth.euv.ResponseTable`
を作って渡す。表は ``# logT ch1 ch2 ...`` の見出しと数の列のテキスト (``ResponseTable.save``)。
観測した日付の AIA の劣化や XRT の汚染に合わせるときは、SSW のファイルを手元に落として作り直す。

.. code:: bash

    # https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/response/ の aia_V9_all_fullinst.genx, aia_V9_fullemiss.genx (790 MB),
    # aia_V9_chiantifix.genx, aia_V10_20201119_190000_response_table.txt を DIR に置く (読むのに sunpy が要る)
    python -m pyR2D2.synth.euv aia --ssw DIR --out aia_20240101.txt --date 2024-01-01
    # https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/response/ の channels/xrt_channels_v0017.geny,
    # emiss/XRT_emiss_model.default_CHIANTI_20240306.geny, contam/ の xrt_contam_on_ccd.geny,
    # xrt_contam_on_filter.geny, n_DEHP.txt を DIR に置く
    python -m pyR2D2.synth.euv xrt --ssw DIR --out xrt_20240101.txt --date 2024-01-01

EM の無いラン (古いラン、球座標など) や、箱を変えたいときは、3D 出力から作り直せる。EOS 表に ``ne`` が要る
(R2D2plus-input の ``stratification_eos/add_electron_density.py`` が足す。``d.p.read_eos_table`` が読むと
``d.p.log_ne_e`` ができる)。

.. code:: python

    em = euv.em_from_state(d, nd, los="x")         # nd は 3D の番号。R2D2plus の vl_dem と float32 の精度で一致する
