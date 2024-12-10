R2D2 pythonでのキーワードの説明
=========================================

以下では、R2D2 pythonで使われている辞書型に含まれるキーの説明を行う

* キーの名前 (型) -- 説明 [単位]

というフォーマットを採用する。

R2D2では、:code:`R2D2_data` を提供している。

self.p [dictionary]
--------------------------------

.. code::

    import pyR2D2
    self = pyR2D2.Data(datadir)

とすると、初期設定が読み込まれる。 :code:`self` は :code:`pyR2D2.Data` のオブジェクトであり、名前は任意である。 :code:`init.py` や :code:`mov.py` では、オブジェクト名は :code:`d` としてある。

出力・時間に関する量
::::::::::::::::::::::::::::::::

* datadir (str) -- データの保存場所
* nd (int) -- 現在までのアウトプット時間ステップ数(3次元データ)
* nd_tau (int) -- 現在までのアウトプット時間ステップ数(光学的厚さ一定のデータ)
* dtout (float) --　出力ケーデンス [s]
* dtout_tau (float) --　光学的厚さ一定のデータの出力ケーデンス [s]
* ifac (int) -- dtout/dtout_tau
* tend (float) -- 計算終了時間。大きく取ってあるためにこの時間まで計算することはあまりない [s]
* swap (int) --　エンディアン指定。big endianは1、little endianは0。IDLの定義に従っている。
* endian (char) -- エンディアン指定。big endianは > 、little endianは < 。pythonの定義に従っている。
* m_in (int) -- 光学的厚さ一定のデータを出力する変数の数
* m_tu (int) -- 光学的厚さ一定のデータの層の数


座標に関する量
::::::::::::::::::::::::::::::::
* xdcheck (int) -- x軸方向に解いているか。解いていたら2、解いていなかったら1
* ydcheck (int) -- y軸方向に解いているか。解いていたら2、解いていなかったら1
* zdcheck (int) -- z軸方向に解いているか。解いていたら2、解いていなかったら1    
* margin (int) -- マージン(ゴーストセル)の数
* nx (int) -- 1 MPIスレッドあたりのx方向の格子点の数
* ny (int) -- 1 MPIスレッドあたりのy方向の格子点の数
* nz (int) -- 1 MPIスレッドあたりのz方向の格子点の数
* ix0 (int) -- x方向のMPI領域分割の数
* jx0 (int) -- y方向のMPI領域分割の数
* kx0 (int) -- z方向のMPI領域分割の数
* ix (int) -- x方向の格子点数 ix0*nx
* jx (int) -- y方向の格子点数 jx0*ny
* kx (int) -- z方向の格子点数 kx0*nz
* npe (int) -- 全MPIスレッドの数 ``npe = ix0*jx0*kx0``
* mtype (int) -- 変数の数
* xmax (float) -- x方向境界の位置(上限値)　[cm]
* xmin (float) -- x方向境界の位置(下限値)　[cm]
* ymax (float) -- y方向境界の位置(上限値)　[cm]
* ymin (float) -- y方向境界の位置(下限値)　[cm]
* zmax (float) -- z方向境界の位置(上限値)　[cm]
* zmin (float) -- z方向境界の位置(下限値)　[cm]
* x (float) [ix] -- x方向の座標 [cm]
* y (float) [jx] -- y方向の座標 [cm]
* z (float) [kx] -- z方向の座標 [cm]
* xr (float) [ix] -- x/rsun
* xn (float) [ix] -- ``(x-rsun)*1.e-8``
* deep_top_flag (int) --
* ib_excluded_top (int) --
* rsun (float) [ix] -- 太陽半径 [cm]

背景場に関する量
::::::::::::::::::::::::::::::::
* pr0 (float) [ix] -- 背景場の圧力 [dyn cm `-2`:sup:]
* te0 (float) [ix] -- 背景場の温度 [K]
* ro0 (float) [ix] -- 背景場の密度 [g cm `-3`:sup:]
* se0 (float) [ix] -- 背景場のエントロピー [erg g `-1`:sup: K `-1`:sup:]
* en0 (float) [ix] -- 背景場の内部エネルギー [erg cm `-3`:sup:]
* op0 (float) [ix] -- 背景場のオパシティー [g `-1`:sup: cm `-2`:sup:]
* tu0 (float) [ix] -- 背景場の光学的厚さ
* dsedr0 (float) [ix] -- 背景場の鉛直エントロピー勾配 [erg g `-1`:sup: K `-1`:sup: cm `-1`:sup:]
* dtedr0 (float) [ix] -- 背景場の鉛直温度勾配 [K cm `-1`:sup:]
* dprdro (float) [ix] -- 背景場の :math:`(\partial p/\partial \rho)_s` 
* dprdse (float) [ix] -- 背景場の :math:`(\partial p/\partial s)_\rho` 
* dtedro (float) [ix] -- 背景場の :math:`(\partial T/\partial \rho)_s` 
* dtedse (float) [ix] -- 背景場の :math:`(\partial T/\partial s)_\rho`
* dendro (float) [ix] -- 背景場の :math:`(\partial e/\partial \rho)_s` 
* dendse (float) [ix] -- 背景場の :math:`(\partial e/\partial s)_\rho` 
* gx (float) [ix] -- 重力加速度 [cm s `-2`:sup:]
* kp (float) [ix] -- 放射拡散係数 [cm `2`:sup: s `-1`:sup:]
* cp (float) [ix] -- 定圧比熱　[erg g `-1`:sup: K `-1`:sup:]
* fa (float) [ix] -- 対流層の底付近の輻射によるエネルギーフラックス。光球付近では輻射輸送を直に解くために含まれないが、上部境界が光球にない場合は、上部境界付近の人工的なエネルギーフラックス(冷却が含まれる) [erg cm `-2`:sup:]
* sa (float) [ix] -- 上記faによる加熱率 [erg cm `-3`:sup:]
* xi (float) [ix] -- 音速抑制率
* ix_e (int) -- 状態方程式の密度の格子点数
* jx_e (int) -- 状態方程式のエントロピーの格子点数

解析のためのデータ再配置(remap)に関する量
::::::::::::::::::::::::::::::::::::::::::::

* m2da (int) -- remapで出力した解析量の数
* cl (char) [m2da] --　remapで出力した解析量の名前
* jc (int) -- ``self.vc['vxp']`` などで出力するスライスのy方向の位置
* kc (int) -- 浮上磁場の中心と思っている場所を出力(あまり使わない)
* ixr (int) -- remap後のx方向分割の数
* jxr (int) -- remap後のy方向分割の数
* iss (int) [npe] -- remap後配列のそれぞれのMPIプロセスのx方向の初めの位置
* iee (int) [npe] -- remap後配列のそれぞれのMPIプロセスのx方向の終わりの位置
* jss (int) [npe] -- remap後配列のそれぞれのMPIプロセスのy方向の初めの位置
* jee (int) [npe] -- remap後配列のそれぞれのMPIプロセスのy方向の終わりの位置
* iixl (int) [npe] -- remap後配列のそれぞれのMPIプロセスのx方向の格子点数
* jjxl (int) [npe] -- remap後配列のそれぞれのMPIプロセスのy方向の格子点数
* np_ijr (int) [npe] -- x, y方向のMPIプロセスの位置を入力するとMPIプロセス番号を返す配列
* ir (int) [npe] -- MPIプロセス番号を入れるとx方向のMPIプロセスの位置を返す配列
* jr (int) [npe] -- MPIプロセス番号を入れるとy方向のMPIプロセスの位置を返す配列
* i2ir (int) [ix] -- x方向の格子点の位置を入れるとx方向のMPIプロセスの位置を返す配列
* j2jr (int) [jx] -- y方向の格子点の位置を入れるとy方向のMPIプロセスの位置を返す配列

スライスデータに関する量
::::::::::::::::::::::::::::::::::::::::::::

* nx_slice [int] -- x一定面のスライスの数
* ny_slice [int] -- y一定面のスライスの数
* nz_slice [int] -- z一定面のスライスの数
* x_slice [float] -- x一定面のスライスの位置 [cm]
* y_slice [float] -- y一定面のスライスの位置 [cm]
* z_slice [float] -- z一定面のスライスの位置 [cm]

self.qs [dictionary]
--------------------------------

.. code::
    
    xs = 0.99*rsun
    ns = 10
    self.read_qq_select(xs,ns)

として高さ :code:`xs` での二次元スライスを読み込む

* ro (float) [jx,kx] -- 密度の擾乱 :math:`\rho_1` [g cm `-3`:sup:]
* vx (float) [jx,kx] -- x方向の速度 :math:`v_x` [cm s `-1`:sup:]
* vy (float) [jx,kx] -- y方向の速度 :math:`v_y` [cm s `-1`:sup:]
* vz (float) [jx,kx] -- z方向の速度 :math:`v_z` [cm s `-1`:sup:]
* bx (float) [jx,kx] -- x方向の磁場 :math:`B_x` [G]
* by (float) [jx,kx] -- y方向の磁場 :math:`B_y` [G]
* bz (float) [jx,kx] -- z方向の磁場 :math:`B_z` [G]
* se (float) [jx,kx] -- エントロピーの擾乱 :math:`s_1` [erg g `-1`:sup: K `-1`:sup:]
* pr (float) [jx,kx] -- 圧力の擾乱 :math:`p_1` [dyn cm `-2`:sup:]
* te (float) [jx,kx] -- 温度の擾乱 :math:`T_1` [K]
* op (float) [jx,kx] -- 不透明度(オパシティー) :math:`\kappa` [g `-1`:sup: cm `-2`:sup:]

self.qq [dictionary]
--------------------------------

:code:`self.qs` と同様

self.qt [dictionary]
--------------------------------

ほぼself.qsと同様だが、以下の追加量が保存してある。

self.vc [dictionary]
--------------------------------

数値計算実行時に解析・出力している統計量。しばしばバグがあるので注意すること。

.. code::

    ns = 10
    self.read_vc(ns)

として統計量を読み込む.

* su, sd (float) [ix,jx] -- ある動径位置 :math:`r`, 余緯度 :math:`\theta` における上昇流(su), 下降流(sd)の格子点数 [個]
* :code:`**m` と表されるものは経度方向平均。以下の物理量がある
    * rom (float) [ix,jx] -- 密度　[g cm `-3`:sup:]
    * vxm (float) [ix,jx] -- x方向の速度 [cm s `-1`:sup:]
    * vym (float) [ix,jx] -- y方向の速度 [cm s `-1`:sup:]
    * vzm (float) [ix,jx] -- z方向の速度 [cm s `-1`:sup:]
    * rxm (float) [ix,jx] -- x方向の運動量 ::math:`\rho v_x`  [g cm `-2`:sup: s `-1`:sup:]
    * rym (float) [ix,jx] -- y方向の運動量 ::math:`\rho v_x`  [g cm `-2`:sup: s `-1`:sup:]
    * rzm (float) [ix,jx] -- z方向の運動量 ::math:`\rho v_x`  [g cm `-2`:sup: s `-1`:sup:]
    * bxm (float) [ix,jx] -- x方向の磁場 [G]
    * bym (float) [ix,jx] -- y方向の磁場 [G]
    * bzm (float) [ix,jx] -- z方向の磁場 [G]
    
.. * ``**rms``と表されるものは、経度に対するRMS量

最終更新日：|today|