=================
pyR2D2
=================

Introduction
----------------------

Python module for dealing with calculation results from R2D2 code is available.

R2D2で実施した計算を取り扱うためのPythonモジュールを公開している。

https://github.com/hottahd/R2D2_py

You can install it by cloning the repository and running pip install as

以下のようにしてインストールできる。

.. code-block:: bash

    git clone https://github.com/hottahd/pyR2D2
    cd pyR2D2
    pip install .

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