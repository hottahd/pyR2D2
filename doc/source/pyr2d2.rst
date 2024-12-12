=================
pyR2D2
=================

Introduction
----------------------

R2D2で実施した計算を取り扱うためのPythonモジュールを公開している。

https://github.com/hottahd/R2D2_py

.. code-block:: bash

    git clone https://github.com/hottahd/pyR2D2
    cd pyR2D2
    pip install .

としてインストールする。
PythonでR2D2で定義された関数を使うには

.. code:: python

    import pyR2D2

とする。すべてのAPIは :doc:`/api_reference` に示されている。

:class:`pyR2D2.Data`
--------------------------

R2D2には :py:class:`pyR2D2.Data` クラスが定義してあり、これをオブジェクト指向的に用いてデータを取り扱う。

.. code:: python

    import pyR2D2
    datadir = '../run/d002/data'
    d = pyR2D2.Data(datadir)

などとしてインスタンスを生成する。

.. code:: python

    help(pyR2D2)
    help(pyR2D2.Data)

などとすると実行環境で、モジュール全体や各関数の簡単な説明を見ることができる。            

最終更新日：|today|