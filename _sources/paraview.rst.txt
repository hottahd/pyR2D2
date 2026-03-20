Paraviewを用いた3Dデータ表示
================================

ここでは、Paraviewを用いてR2D2の計算結果を三次元表示する方法を説明する。

データ準備
--------------------------------
データをParaviewで扱うためにVTKフォーマットに変換する。
変換後に用意するべきデータは

- ある物理量の三次元データ
- ある物理量の ::math:`\tau=1` でのデータ

それぞれのデータのためにpyR2D2では以下の関数が用意してある。

- :meth:`pyR2D2.write.vtk.write_3D`
- :meth:`pyR2D2.write.vtk.write_3D_vector`
- :meth:`pyR2D2.write.vtk.write_optical_surface`

例えば、以下のようにして実行する

.. code:: python

    d = pyR2D2.read.R2D2Data('data/') # dataディレクトリにあるデータを読込
    d.qf.read(100, ['bx', 'by', 'bz']) # 100番目の3次元データを読込
    bb = sqrt(d.qq['bx']**2 + d.qq['by']**2 + d.qq['bz']**2) # 磁場の強さを計算
    d.qt.read(100) # 100番目のtau=1の2次元データを読込

    pyR2D2.write.vtk.write_3D(bb, d.x, d.y, d.z,'bb.vtk','bb')
    # 変数名をbbとしてファイル名bb.vtkに3次元データを出力
    pyR2D2.write.write_optical_surface(d.qt.rt, d.qt.he, d.y, d.z, 'in.vtk', 'in')
    # 変数名をinとしてファイル名in.vtkに2次元データを出力

Paraviewを用いて3次元表示
--------------------------------
`Paraviewのサイト <https://www.paraview.org/download/>`_ からParaviewをダウンロード。Windows, Linux, macOSのそれぞれのソフトウェアがあるのでインストール方法は各自確認すること。ここでは、macOSでの利用方法を示すが、確認している限りは、Linuxでもほとんど同じ。ここでは非常に簡単にParaviewの使い方を説明する。詳しくはParaviewの公式マニュアルなどを読むこと。

1. まず右上のファイルアイコンをクリック

    .. image:: _static/figs/paraview/paraview01.png
        :width: 600 px

2. Pythonで生成したファイルを選択。2次元, ３次元ファイルどちらも選択する。一個一個選択しても良いし、一度に選択しても良い。時系列データの時は、すべてを一度に選択するとアニメーションを作りやすい。

    .. image:: _static/figs/paraview/paraview02.png
        :width: 600 px

3. Applyをクリック。選択した二つのデータが表示される。

    .. image:: _static/figs/paraview/paraview03.png
        :width: 600 px

4. 2次元データの方は、すぐに面として表示されるが、三次元データは表示方法を選ぶ必要がある。

    .. image:: _static/figs/paraview/paraview04.png
        :width: 600 px

5. 三次元データのボリュームレンダリングが行いたいので、Volumeを選ぶ。

    .. image:: _static/figs/paraview/paraview05.png
        :width: 600 px

6. 三次元データのボリュームレンダリングが表示されるので、便宜描画を回転させるなどして、解析する。

    .. image:: _static/figs/paraview/paraview06.png
        :width: 600 px

最終更新日：|today|