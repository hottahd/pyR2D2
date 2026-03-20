.. R2D2 documentation master file, created by
   sphinx-quickstart on Fri Dec 20 15:59:33 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

R2D2 Document
================================

このページは太陽と恒星のための輻射磁気流体コードR2D2 (RSST and Radiation for Deep Dynamics)のマニュアルである。 

R2D2では、輻射磁気流体力学の方程式を中央差分と非線形の人工粘性を用いて高精度かつ安定に解く。R2D2では、太陽表面付近では現実的に輻射輸送を解き、太陽深部では微小擾乱を正確に解く手法を取ることで、太陽内部から表面までを一貫して解くことができる。

.. toctree::
   :maxdepth: 2

   quick_start
   environment
   typical_case
   pyr2d2
   api_reference
   sphinx.rst
   paraview.rst
   geometry.rst
	     
ライセンス
----------------------------------

R2D2は公開ソフトウェアではなく、再配布も禁じている。
共同研究者のみが使って良いというルールになっており、R2D2の使用には、以下の規約を守る必要がある。

* 再配布しない
* 改変は許されるが、その時の実行結果について堀田は責任を持たない
* R2D2で実行する計算は、堀田と議論する必要がある。パラメタ変更などの細かい変更には相談する必要はないが、新しいプロジェクトを開始するときはその都度相談すること。堀田自身のプロジェクト、堀田の指導学生のプロジェクトとの重複を避けるためである。
* R2D2を用いた論文を出版するときは `Hotta et al., 2019 <https://ui.adsabs.harvard.edu/abs/2019SciA....5.2307H/abstract>`_, `Hotta and Iijima, 2020 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.2523H/abstract>`_ を引用すること。  
* R2D2を用いた研究を発表するときは、`R2D2のロゴ <https://hottahd.github.io/R2D2-manual/_images/R2D2_logo_red.png>`_ の使用が推奨される(強制ではない)。

出版論文
----------------------------------
R2D2を用いた研究で出版された論文は以下です。

1. `Hotta, Iijima, and Kusano, 2019, Science Advances, 5, eaau2307 <https://ui.adsabs.harvard.edu/abs/2019SciA....5.2307H/abstract>`_
2. `Toriumi and Hotta, 2019, ApJ, 886, L21 <https://ui.adsabs.harvard.edu/abs/2019ApJ...886L..21T/abstract>`_
3. `Hotta and Iijima, 2020, MNRAS, 494, 2523 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.2523H/abstract>`_
4. `Hotta and Toriumi, 2020, MNRAS, 498, 2925 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.498.2925H/abstract>`_
5. `Takahata, Hotta, Iida, and Oba, MNRAS, 503, 3610 <https://ui.adsabs.harvard.edu/abs/2021MNRAS.503.3610T/abstract>`_
6. `Hotta and Kusano, 2021, Nature Astronomy, 5, 1100 <https://ui.adsabs.harvard.edu/abs/2021NatAs...5.1100H/abstract>`_
7. `Hotta, Kusano, and Shimada, 2022, ApJ, 933, 2 <https://ui.adsabs.harvard.edu/abs/2022ApJ...933..199H/abstract>`_
8. `Mori and Hotta, 2023, MNRAS, 519, 3091 <https://ui.adsabs.harvard.edu/abs/2023MNRAS.519.3091M/abstract>`_
9. `Kaneko, Hotta, Toriumi and Kusano, 2022, MNRAS, 517, 2775 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.517.2775K/abstract>`_
10. `Mori and Hotta, 2023, MNRAS, 524, 4746 <https://ui.adsabs.harvard.edu/abs/2023MNRAS.524.4746M/abstract>`_
11. `Toriumi, Hotta and Kusano, 2023, Scientific Report, 13, 8994 <https://ui.adsabs.harvard.edu/abs/2023NatSR..13.8994T/abstract>`_
12. `Masaki, Hotta, Katsukawa and Ishikawa, 2023, PASJ, <https://ui.adsabs.harvard.edu/abs/2023PASJ..tmp...76M/abstract>`_
13. `Silva et al., 2023, ApJ, 948, L24 <https://ui.adsabs.harvard.edu/abs/2023ApJ...948L..24S/abstract>`_
14. `Toriumi, Hotta and Kusano, 2024, ApJ, 975, 209, <https://ui.adsabs.harvard.edu/abs/2024ApJ...975..209T/abstract>`_

賞
----------------------------------
* `HPCI ソフトウェア賞 <https://www.hpci-c.jp/hrdevelop/award.html>`_

最終更新日：|today|