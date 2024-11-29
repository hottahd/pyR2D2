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
読み込みについては、PythonコードとIDLコードを用意しているが、開発の頻度が高いPythonコードの利用を推奨している。

Pythonコード
::::::::::::::::::::::

.. py:module:: R2D2

PythonでR2D2で定義された関数を使うには

.. code:: python

    import R2D2

として、モジュールを読み込む。R2D2には :py:class:`R2D2_data` クラスが定義してあり、これをオブジェクト指向的に用いてデータを取り扱う。

以下にそれぞれの関数を示すが、docstringは記入してあるので

.. code:: python

    help(R2D2)
    help(R2D2.R2D2_data)

などとすると実行環境で、モジュール全体や各関数の簡単な説明を見ることができる。

クラス
^^^^^^^^^^^^^^^^^^^^^^^
.. py:class:: R2D2_data(datadir)

データの読み込みには ``R2D2`` モジュールの中で定義されている ``R2D2_data`` クラスを使う必要がある。

.. code:: python

    import R2D2
    datadir = '../run/d002/data'
    d = R2D2.R2D2_data(datadir)

などとしてインスタンスを生成する。

Attribute
^^^^^^^^^^^^^^^^^^^^^^^

.. py:attribute:: R2D2_data.p

    基本的なパラメタ。格子点数 ``ix`` や領域サイズ ``xmax`` など。
    インスタンス生成時に同時に読み込まれる。

.. py:attribute:: R2D2_data.t

    時間を格納するデータ。　:py:meth:`R2D2_data.read_time` では戻り値として同じ値を返すために :py:attr:`R2D2_Data.t` はあまり使われない。

.. py:attribute:: R2D2_data.qs

    ある高さの2次元のndarrayが含まれる辞書型。 :py:meth:`R2D2_data.read_qq_select` で読み込んだ結果。

.. py:attribute:: R2D2_data.qq
    
    3次元のnumpy array。計算領域全体のデータ。:py:meth:`R2D2_data.read_qq` で読み込んだ結果。

.. py:attribute:: R2D2_data.qt

    2次元のnumpy array。ある光学的厚さの面でのデータ。現在は光学的厚さ1, 0.1, 0.01でのデータを出力している。 :py:meth:`R2D2_data.read_qq_tau` で読み込んだ結果。

.. py:attribute:: R2D2_data.vc

    Fortranの計算の中で解析した結果。 :py:meth:`R2D2_data.read_vc` で読み込んだ結果。
.. py:attribute:: R2D2_data.qc

    3次元のnumpy array。計算領域全体のデータ。Fortranの計算でチェックポイントのために出力しているデータを読み込む。主に解像度をあげたいときのために使う :py:meth:`R2D2_data.read_qq_check` で読み込んだ結果。

.. py:attribute:: R2D2_data.ql

    2次元のnumpy array。Fortranで定義したスライスデータ :py:meth:`R2D2_data.read_qq_slice` で読み込んだ結果。
    実際にどの位置のスライスを読み込んでいるかは ``R2D2.p['x_slice']``, ``R2D2.p['y_slice']``, ``R2D2.p['z_slice']`` を参照すること。スライスの位置の配列が保存されている。

:py:attr:`R2D2_data.p` については、``init.py`` などで

.. code:: python

    import R2D2
    d = R2D2.R2D2_data(datadir)
    for key in R2D2.p:
        exec('%s = %s%s%s' % (key, 'R2D2.p["',key,'"]'))

などとしているために、辞書型の ``key`` を名前にする変数に値が代入されている。例えば、 ``R2D2_data.p['ix']`` と ``ix`` には同じ値が入っている。

Method
^^^^^^^^^^^^^^^^^^^^^^

データ読み込み
||||||||||||||||||||||||

メソッドで指定する ``datadir`` はデータの場所を示す変数。R2D2の計算を実行すると ``data`` ディレクトリが生成されて、その中にデータが保存される。この場所を指定すれば良い。

.. py:method:: R2D2_data.__init__(datadir)
    
    インスタンス生成時に実行されるメソッド。計算設定などのパラメタが読み込まれる。 :py:attr:`R2D2_data.p` にデータが保存される。
    
    :argument str datadir: データの場所

.. py:method:: R2D2_data.read_qq_select(xs,n,silent=False)
    
    ある高さのデータのスライスを読み込む。戻り値を返さない時も :py:attr:`R2D2_data.qs` にデータが保存される。

    :param float xs: 読み込みたいデータの高さ
    :param int n: 読み込みたい時間ステップ
    :param bool silent: 読み込み時にどの変数に格納されたかメッセージの表示を抑制するフラグ。デフォルトはFalseで、Trueだとメッセージは表示されない。

.. py:method:: R2D2_data.read_qq(n,silent=False)
    
    3次元のデータを読み込む。戻り値を返さない時も :py:attr:`R2D2_data.qq` にデータが保存される。

    :param int n: 読み込みたい時間ステップ
    :param bool silent: 読み込み時にどの変数に格納されたかメッセージの表示を抑制するフラグ。デフォルトはFalseで、Trueだとメッセージは表示されない。

.. py:method:: R2D2_data.read_qq_tau(n,silent=False)
    
    光学的厚さが一定の2次元のデータを読み込む。:py:attr:`R2D2_data.qt` にデータが保存される。

    :param int n: 読み込みたい時間ステップ
    :param bool silent: 読み込み時にどの変数に格納されたかメッセージの表示を抑制するフラグ。デフォルトはFalseで、Trueだとメッセージは表示されない。

.. py:method:: R2D2_data.read_time(n,tau=False,silent=False)
    
    時間を読み込む。 :py:attr:`R2D2.t` にもデータは格納されるが戻り値としても使うことができる。

    :param int n: 読み込みたい時間ステップ
    :param bool tau: 光学的厚さ一定のデータ(高ケーデンス)に対する時間を読むかのフラグ。デフォルトはFalse。
    :param bool silent: 読み込み時にどの変数に格納されたかメッセージの表示を抑制するフラグ。デフォルトはFalseで、Trueだとメッセージは表示されない。
    :return: 時間ステップでの時間

.. py:method:: R2D2_data.read_vc(n,silent=False)
    
    Fortranコードの中で解析した計算結果を読み込む。戻り値を返さない時も :py:attr:`R2D2_data.vc` にデータが保存される。

    :param int n: 読み込みたい時間ステップ
    :param bool silent: 読み込み時にどの変数に格納されたかメッセージの表示を抑制するフラグ。デフォルトはFalseで、Trueだとメッセージは表示されない。

.. py:method:: R2D2_data.read_qq_check(n,silent=False,end_step=False)
    
    チェックポイントのための3次元データを読み込む。主に解像度をあげるときに使う。 :py:attr:`R2D2_data.qc` にデータが保存される。

    :param int n: 読み込みたい時間ステップ
    :param bool silent: 読み込み時にどの変数に格納されたかメッセージの表示を抑制するフラグ。デフォルトはFalseで、Trueだとメッセージは表示されない。
    :param book end_step: Falseの時は、 ``n`` で指定された時間ステップのデータを読み込むが3次元データはそれほど高頻度ではな出力していない。Trueの時は、 ``qq.dac.e`` もしくは　``qq.dac.o`` という最後の1ステップの出力データを読み込む。こちらは常に上書きされてしまっているが、毎ステップ必ず書き込むので最後のステップのデータを読み込みたい時はこちらをTrueにする。

.. py:method:: R2D2_data.read_qq_slice(n,n_slice,direc,silent=False)
    
    ``slice_def.F90`` で指定したスライスデータを読み込む。:py:attr:`R2D2_data.ql` にデータが保存される。

    :param int n: 読み込みたい時間ステップ
    :param int n_slice: 何枚目のスライスを読み込むか
    :param str direc: スライスの方向 'x', 'y', 'z'のどれか
    :param bool silent: 読み込み時にどの変数に格納されたかメッセージの表示を抑制するフラグ。デフォルトはFalseで、Trueだとメッセージは表示されない。    

データダウンロード
||||||||||||||||||||||||

スパコンなどで計算した後に、ローカルの環境にデータをダウンロードするメソッドも提供している。堀田と全く同じようにディレクトリ構造を作ってないといけないので注意。

.. py:function:: R2D2.sync.set(server,caseid,project=os.getcwd().split('/')[-2],dist='../run/')

    設定のみをダウンロードするメソッド。ひとまずGoogleスプレッドシートに書き込みたい時などに有用。

    :param char server: ダウンロード先のサーバー名。sshで使うサーバー名を用いれば良い。
    :param char server: ダウンロードしたいcaseid。'd001'などとする。
    :param char project: プロジェクト名。'R2D2'など。何も入力しなければ一個上のディレクトリの名前をプロジェクト名と判断する。
    :param char dist: データダウンロード先。特別な用途がなければデフォルトのままとする。
    
.. py:method:: R2D2_data.sync_tau(server,project=os.getcwd().split('/')[-2])

    光学的厚さ一定の面上でのデータをダウンロードするメソッド。

    :param char server: ダウンロード先のサーバー名。sshで使うサーバー名を用いれば良い。
    :param char project: プロジェクト名。'R2D2'など。何も入力しなければ一個上のディレクトリの名前をプロジェクト名と判断する。

.. py:method:: R2D2_data.sync_select(xs,server,project=os.getcwd().split('/')[-2])

    ２次元データをダウンロードするメソッド

    :param float xs: ダウンロードする高さ。
    :param char server: ダウンロード先のサーバー名。sshで使うサーバー名を用いれば良い。
    :param char project: プロジェクト名。'R2D2'など。何も入力しなければ一個上のディレクトリの名前をプロジェクト名と判断する。

.. py:method:: R2D2_data.sync_vc(server,project=os.getcwd().split('/')[-2])

    計算実行中に解析した物理量をダウンロードするメソッド

    :param char server: ダウンロード先のサーバー名。sshで使うサーバー名を用いれば良い。
    :param char project: プロジェクト名。'R2D2'など。何も入力しなければ一個上のディレクトリの名前をプロジェクト名と判断する。

.. py:method:: R2D2_data.sync_check(n,server,project=os.getcwd().split('/')[-2],end_step=False)

    チェックポイントデータをダウンロードするメソッド

    :param int n: ダウンロードしたい時間ステップ。
    :param char server: ダウンロード先のサーバー名。sshで使うサーバー名を用いれば良い。
    :param char project: プロジェクト名。'R2D2'など。何も入力しなければ一個上のディレクトリの名前をプロジェクト名と判断する。

解像度・計算領域変更
::::::::::::::::::::::

R2D2のPythonの機能を用いて, 解像度や計算領域を変更することができる.

以下の手順に従う

1. fortranコードで何らかの計算を実行
2. pythonで読み込み。解像度変換を実行
3. fortranで再度, 計算を実行

pythonでの解像度変換には :py:meth:`R2D2_data.upgrade_resolution` メソッドを用いる.

.. py:method:: R2D2_data.upgrade_resolution(caseid, n, xmin, xmax, ymin, ymax, zmin, zmax, ix, jx, kx, ix_ununi=32, dx00=48e5, x_ununif=False, endian='<', end_step=False, memory_saving=False)

    データの解像度や計算領域を変更するためのメソッド

    :param char caseid: 出力先のcaseid e.g. 'd002'
    :param int n: 何番目のデータの解像度・計算領域を変換するか. データが必ずしもあるとは限らないので `end_step=True` が推奨される.
    :param float xmax: max location in x direction
    :param float xmin: min location in x direction
    :param float ymax: max location in y direction
    :param float ymin: min location in y direction
    :param float zmax: max location in z direction
    :param float zmin: min location in z direction
    
    :param char endian: endian, "<" もしくは, ">"
    
    :param int ix: updated grid point in x direction
    :param int jx: updated grid point in y direction
    :param int kx: updated grid point in z direction
    :param bool memory_saving: If true, upgraded variables are saved in memory separately for saving memory
    :param bool end_step: `end_step=True` のときは, パラメータ `n` は無視されて持っている一番最後のステップのデータの解像度・計算領域変更がされる. `end_step=False` のときは　`n` ステップのデータが利用される.cd .
        
    これより下のパラメタは `x_ununif=True` を用いたときのみ有効となる.    
    
    :param int ix_ununi: number of grid in uniform grid region
    :param float dx00: grid spacing in uniform grid region
    :param bool x_ununif: whethere ununiform grid is used
    


例えば, `caseid='d001'` のデータの解像度を変更して `caseid='d002'` へと出力する時は

.. code:: python

    d = R2D2.R2D2_data('../run/d001')
    d.upgrade_resolution('d002',...)

として, 出力された結果を参考にd002のプログラムを変更する.

IDLコード
::::::::::::::::::::::

`GitHubの公開レポジトリ <https://github.com/hottahd/R2D2_idl>`_ に簡単な説明あり

バージョン履歴
----------------------

* ver. 1.0: バージョン制を導入
* ver. 1.1: 光学的厚さが0.1, 0.01の部分も出力することにした。qq_in, vcをconfigのグローバル変数として取扱うことにした。
* ver. 1.2: データ構造を変更。

最終更新日：|today|