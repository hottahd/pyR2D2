=================
pyR2D2
=================

Introduction
----------------------

R2D2で実施した計算を取り扱うためのPythonモジュールを公開している。

https://github.com/hottahd/R2D2_py

以下でダウンロードできる。

.. code-block:: bash

    git clone https://github.com/hottahd/pyR2D2

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

:class:`pyR2D2.Data` classは :class:`pyR2D2.Read` と :class:`pyR2D2.Sync` の2つのクラスを以下のようにコンポジットしている。

.. code:: python

    class Data:
        (some codes...)        

        def __init__(self, datadir):
            self.read = Read(...)
            self.sync = Sync(...)

            (some codes...)        


つまり、 :class:`pyR2D2.Data` クラスは、 :class:`pyR2D2.Read` クラスと :class:`pyR2D2.Sync` クラスのメソッドをすべて持っており、それぞれ :class:`pyR2D2.Data.read` と :class:`pyR2D2.Data.sync` でアクセスできる。

一方、 :class:`pyR2D2.Read` の属性(attribute)は、 :class:`pyR2D2.Data` の属性としてもアクセスできる。
:class:`pyR2D2.Read` の属性は頻繁に使うので、このような仕様としている。一方、 :class:`pyR2D2.Read` のメソッドは、 :class:`pyR2D2.Data` のメソッドとしてはアクセスできない。
            

最終更新日：|today|