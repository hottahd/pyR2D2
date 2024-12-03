出力と読込
=================

出力
----------------------

Fortranコード
::::::::::::::::::::::

Sliceデータ
^^^^^^^^^^^^^^^^^^^^^^^^^^^
R2D2では、スライスデータも高いケーデンスで出力できるようになっている。多くのデータを出力するのには効率の悪い方法になっているので、3-4枚のスライスを出力するのに留めておくことが推奨される。 ``slice_def.F90`` 内でどの部分のスライスを出力するかを定義している。

.. code:: fortran

    integer, parameter :: nx_slice = 3 ! number of slice in x direction
    integer, parameter :: ny_slice = 2 ! number of slice in y direction
    integer, parameter :: nz_slice = 2 ! number of slice in z direction

この部分で、それぞれの方向に何枚のスライスを出力するかを定義。例えば ``nx_slice`` はスライスするy-z平面の数となる。実際にどの部分を出力するかは続く部分で指定している。

.. code:: fortran

    real(KIND(0.d0)), dimension(nx_slice), save :: x_slice = (/xmin,rsun,xmax/)
    real(KIND(0.d0)), dimension(ny_slice), save :: y_slice = (/ymin,0.5d0*(ymin + ymax)/)
    real(KIND(0.d0)), dimension(nz_slice), save :: z_slice = (/zmin,0.5d0*(zmin + zmax)/)

``nx_slice`` で定義した数と合うように個数を指定しなければならない。

読込
----------------------
読み込みは、PythonコードとIDLコードを用意しているが、開発の頻度が高いPythonコードの利用を推奨している。

Pythonコード
::::::::::::::::::::::

PythonでR2D2で定義された関数を使うには

.. code:: python

    import R2D2

として、モジュールを読み込む。R2D2には :py:class:`R2D2_data` クラスが定義してあり、これをオブジェクト指向的に用いてデータを取り扱う。

以下にそれぞれの関数を示すが、docstringは記入してあるので

.. code:: python

    help(R2D2)
    help(R2D2.R2D2_data)

などとすると実行環境で、モジュール全体や各関数の簡単な説明を見ることができる。


データの読み込みには ``R2D2`` モジュールの中で定義されている ``R2D2_data`` クラスを使う必要がある。

.. code:: python

    import R2D2
    datadir = '../run/d002/data'
    d = R2D2.R2D2_data(datadir)

などとしてインスタンスを生成する。

データダウンロード
||||||||||||||||||||||||

スパコンなどで計算した後に、ローカルの環境にデータをダウンロードするメソッドも提供している。堀田と全く同じようにディレクトリ構造を作ってないといけないので注意。

解像度・計算領域変更
::::::::::::::::::::::

R2D2のPythonの機能を用いて, 解像度や計算領域を変更することができる.

以下の手順に従う

1. fortranコードで何らかの計算を実行
2. pythonで読み込み。解像度変換を実行
3. fortranで再度, 計算を実行

pythonでの解像度変換には :py:meth:`R2D2_data.upgrade_resolution` メソッドを用いる.


IDLコード
::::::::::::::::::::::

`GitHubの公開レポジトリ <https://github.com/hottahd/R2D2_idl>`_ に簡単な説明あり

バージョン履歴
----------------------

* ver. 1.0: バージョン制を導入
* ver. 1.1: 光学的厚さが0.1, 0.01の部分も出力することにした。qq_in, vcをconfigのグローバル変数として取扱うことにした。
* ver. 1.2: データ構造を変更。

最終更新日：|today|