方程式
=================

R2D2で解く方程式は以下である。現状では、デカルト座標 :math:`(x,y,z)`  と球座標 :math:`(r,\theta,\phi)` を提供している。球座標を使う場合は、Yin-Yang格子を用いることもできる。数値計算コードの中では、:math:`x` を重力方向(鉛直方向、動径方向)に取っているが、論文を書く際は各自適切に判断されたい。

磁気流体力学
-----------------
R2D2で解いている磁気流体力学の方程式は以下である。


.. math::

    \frac{\partial \rho_1}{\partial t} &= - \frac{1}{\xi^2}\nabla\cdot
    \left(\rho \boldsymbol{v}\right) \\
    \frac{\partial}{\partial t}\left(\rho \boldsymbol{v}\right) &=
    -\nabla\cdot\left(\rho\boldsymbol{vv}\right)
    - \nabla p_1 - \rho_1 g\boldsymbol{e}_x
    +\frac{1}{4\pi}\left(\nabla\times\boldsymbol{B}\right)
    \times\boldsymbol{B} \\
    \frac{\partial \boldsymbol{B}}{\partial t} &= 
    \nabla\times\left(\boldsymbol{v\times B}\right)
    \\
    \rho T \frac{\partial s_1}{\partial t} &= -\rho T 
    \left(\boldsymbol{v}\cdot\nabla\right) s + Q_\mathrm{rad} \\
    p_1 &= p_1(\rho_1,s_1,x)

ここで :math:`\rho` は密度、:math:`\boldsymbol{v}` は流体速度、:math:`\boldsymbol{B}` は磁場、 :math:`s` はエントロピー、:math:`p` は圧力、 :math:`T` は温度、 :math:`g` は重力加速度、 :math:`Q_\mathrm{rad}` は輻射による加熱率である。

R2D2では熱力学量を以下のように時間的に一定で :math:`x` 方向の依存性のみを持つ0次の量とそこから擾乱の1次の量に分けている。

.. math::

    \rho &= \rho_0 + \rho_1 \\
    p &= p_0 + p_1 \\
    s &= s_0 + s_1 \\
    T &= T_0 + T_1 \\

太陽内部では、:math:`\rho_1 \ll \rho_0` などが成り立っているが、太陽表面では熱対流による擾乱と背景場は同程度となるので、R2D2の中では :math:`\rho_1 \ll \rho_0` などは仮定しない。0次の量はModel Sを参考にして計算をしている。`Hotta & Iijima2020 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.2523H/abstract>`_ や `Hotta, Kusano and Shimada, 2022 <https://ui.adsabs.harvard.edu/abs/2022ApJ...933..199H/abstract>`_ を参照されたい。

輻射輸送
-----------------

.. todo:: 輻射輸送の方程式

最終更新日：|today|