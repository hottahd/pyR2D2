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
