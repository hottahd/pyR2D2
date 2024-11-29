座標生成
=================

R2D2では中央差分法を用いているが、そのほとんどは数値フラックスを用いて書き直すことでき、提供される ``x`` , ``y``　, ``z`` などは **セル中心** で定義される。よって計算領域内の最初のグリッドは、計算境界から半グリッド進んだところにある。

また、R2D2では一様グリッドと非一様グリッドどちらでも計算できるようにしている。

一様グリッド
--------------------------------

一様グリッドを用いるときは

* 格子間隔を計算する

.. math::

    \Delta x = \frac{x_\mathrm{max} - x_\mathrm{min}}{N_x}

ここで、コードでは、配列の要素数には ``margin`` も含むので
:math:`N_x` を計算するには ``margin`` の部分を引く必要があることに注意。

* :math:`x_1` を設定。``margin`` の分も考慮して計算する。
* ``do loop`` で順次足していく

コードは以下のようになる

.. code:: fortran

    dx_unif = (xmax-xmin)/real(ix00-2*marginx)
    x00(1) = xmin + (0.5d0-dble(marginx))*dx_unif
    if(xdcheck == 2) then
        do i = 1+i1,ix00
            x00(i) = x00(i-i1) + dx_unif
        enddo
    endif    

非一様グリッド
--------------------------------

非一様グリッドを用いるときは、太陽光球付近は、輻射輸送のために一様なグリッド、ある程度の深さから格子間隔が線形に増加する非一様グリッドを使うことにしている。光球近くは、光球をちゃんと解像するために一様グリッド、光球からある程度進むと、非一様グリッドを採用することにしている。
実際の構造は以下のようになっている。非一様グリッド領域の両端２つのグリッド間隔は一様グリッドをとるようにしている。

fortranのコードの中では

* :math:`\Delta x_0` → ``dx00`` : 一様グリッドでの格子点間隔
* :math:`i_\mathrm{x\left(uni\right)}` → ``ix_ununi``: 一様グリッドの格子点数
* :math:`x_\mathrm{ran}` → ``xrange``: 領域サイズ
* :math:`x_\mathrm{ran0}` → ``xrange0``: 一様グリッドの領域サイズ
* :math:`x_\mathrm{ran1}` → ``xrange1``: 非一様グリッドの領域サイズ
* :math:`n_x` → ``nxx`` : 非一様グリッドの格子点数

.. image:: figs/ununiform_grid.png
   :width: 700 px

一様グリッドでのグリッド間隔は :math:`\Delta x_0` として、非一様グリッドでは :math:`\delta x` ずつグリッド間隔が大きくなっていくとする。

.. math::

    x_\mathrm{tran}&={\color{red} \frac{1}{2} \Delta x_0} + {\color{blue} \Delta x_0} + \Delta x_0
    + \left(\Delta x_0 + \delta x\right)
    + \left(\Delta x_0 + 2\delta x\right) + [...] \\
    &+ \left[\Delta x_0 + \left(n_x - 4\right)\delta x\right]
    + {\color{blue}\left[\Delta x_0 + \left(n_x - 4\right)\delta x\right]}
    + {\color{red}\frac{1}{2}\left[\Delta x_0 + \left(n_x - 4\right)\delta x\right]} \\
    &= {\color{red} \Delta x_0 + \frac{1}{2}\left(n_x-4\right)\delta x}
    +{\color{blue} 2\Delta x_0 + \left(n_x - 4\right)\delta x} 
    + \sum_{n=0}^{n_x - 4}\left(\Delta x_0 + n\delta x\right) \\
    &= 3\Delta x_0 + \frac{3\left(n_x-4\right)\delta x}{2}
    + \frac{\left[2\Delta x_0 + \left(n_x-4\right)\delta x\right]\left(n_x - 3\right)}{2} \\
    &= n_x \Delta x_0 + \frac{1}{2} n_x\left(n_x - 4\right)\delta x 

この関係式より、グリッド間隔の増分 :math:`\delta x` を以下のように求めることができる。

.. math::

    \delta x = \frac{2\left(x_\mathrm{tran} - n_x\Delta x_0\right)}{\left(n_x - 4\right)n_x}

最終更新日：|today|