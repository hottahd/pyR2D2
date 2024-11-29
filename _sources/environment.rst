R2D2を使うための環境設定
===============================

R2D2コードを使って計算するだけならば任意のfortranコンパイラ, FFTW, MPIのみがあれば良い。
Pythonコードを使って解析する場合は、いくつかのモジュールが必要なので、そのインストールの方法もここで説明する。

Fortranコードの環境設定
----------------------------------------
Mac
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Homebrewを用いて、必要なコンパイラ・ライブラリをインストールすることを推奨している。
コンパイラとFFTWのインクルードファイルとライブラリの位置だけ指定すれば良いので、
任意の方法でインストールして構わない。Homebrew以外を用いる場合は、便宜make/Makefileを編集すること。

Homebrewのインストール

.. code::

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

gfortranのインストール

.. code::

    brew install gcc

OpenMPIのインストール

.. code::

    brew install openmpi


FFTWのインストール

.. code:: 

    brew install fftw

Linux (Ubuntu 22.04)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ここでは、Ubuntu 22.04の場合のみを説明する。

gfortranのインストール

.. code::

    sudo apt-get install gfortran

OpenMPIのインストール

.. code::

    sudo apt-get install openmpi-doc openmpi-bin libopenmpi-dev

FFTWのインストール

.. code:: 
    
     sudo apt-get install libfftw3-dev

Pythonコードの環境設定
----------------------------------------
Anacondaをインストールし、以下に示すモジュール群をインストールする。
MacとLinuxで共通する部分が多いのでまとめて説明を記す。

`Anacondaのウェブサイト <https://www.anaconda.com/>`_ から対応するインストーラーをダウンロードする。

- Mac
    dmgファイルをダウンロードして、インストール。インストールされるPATHが変わることが多いが、探してPATHを通す。 :code:`/anaconda/bin` や :code:`~/opt/anaconda/bin` など

- Linux
    ダウンロードしてきたシェルスクリプトファイルのあるディレクトリで
    .. code::

        bash ~/Anaconda***.sh

    インストールするディレクトリは :code:`/ホームディレクトリ/anaconda3` とする。
    :code:`/ホームディレクトリ/anaconda3` にPATHを通す。
    スパコンのログインノードなどでもインストール方法は共通である。

ipythonの初期設定
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
以下は必須ではないが、ipythonを使う時の初期設定ファイルである。
:code:`~/.ipython/profile_default/startup/00_first.py`
というファイルを作りそこに以下のように記す。

.. code::

    import sys, os
    import matplotlib.pyplot as plt
    import scipy as sp
    import numpy as np                                                                                                                                                            
    from matplotlib.pyplot import pcolormesh,plot,clf,close
    from numpy import sin,cos,tan,arcsin,arccos,arctan,exp,log,log2,log10,mod,sqrt,absolute,sinh,cosh,tanh,pi,arange
    plt.ion()
    from IPython.core.magic import register_line_magic
    @register_line_magic
    def r(line):
        get_ipython().run_line_magic('run', ' -i ' + line)
    del r                                                                                  
                                      
                                  
最後に記した設定によって、

.. code::

    r (Pythonスクリプト名)

とするだけで、スクリプトを実行できるようになる。

Googleスプレッドシート利用
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

計算設定などをGoogleスプレッドシートにまとめておくと便利である。
R2D2では、Pythonから直接Googleスプレッドシートに送信する方法を提供しているので、利用したい方は検討されたい。

手順については、 `こちら <https://qiita.com/akabei/items/0eac37cb852ad476c6b9>`_ を参考にしたが、少し手順が違うのでこのページでも解説する。

まずは関連するモジュールのインストール。

.. code:: shell

    pip install gspread
    pip install oauth2client

プロキシなどの影響でpipが使えない時は以下のようにする

gspreadのインストール

.. code:: shell

    git clone git@github.com:burnash/gspread.git
    cd gspread
    ipython setup.py install

oauth2clientのインストール

.. code:: shell

    git clone git@github.com:googleapis/oauth2client.git
    cd oauth2client
    ipython setup.py install


プロジェクト作成
||||||||||||||||||||||||||||||

ウェブブラウザで https://console.developers.google.com/cloud-resource-manager?pli=1 にアクセス。

.. image:: figs/google/gen_project1.png
    :width: 350 px

「プロジェクトを作成」として、プロジェクトを作成

.. image:: figs/google/gen_project2.png
    :width: 400 px

プロジェクト名はR2D2, 場所は組織なしとする。

API有効化
||||||||||||||||||||||||||||||

.. image:: figs/google/google_drive1.png
    :width: 400 px

次に検索窓にGoogle Driveと打ち込んで、Google DriveのAPIを検索

.. image:: figs/google/google_drive2.png
    :width: 400 px

Google Drive APIを有効にする。

.. image:: figs/google/google_sheet1.png
    :width: 400 px

同様にGoogle Sheetsと検索

.. image:: figs/google/google_sheet2.png
    :width: 400 px

Google Sheets APIを有効化

サービスアカウント作成
||||||||||||||||||||||||||||||

.. image:: figs/google/service_account1.png
    :width: 400 px

Google APIロゴ → 認証情報 → サービスアカウントとたどる。

.. image:: figs/google/service_account2.png
    :width: 400 px

サービスアカウント名はR2D2とする

.. image:: figs/google/service_account3.png
    :width: 400 px

役割は編集者を選択

.. image:: figs/google/service_account4.png
    :width: 400 px

キーの生成ではJSONを選択し、キーを生成する。
ダウンロードしたファイルは、使用する計算機のホームディレクトリにjsonというディレクトリを作成し、その下に配置する。そのディレクトリには、このjsonファイル以外には何も置かないこと。

スプレッドシート作成
||||||||||||||||||||||||||||||

以下のウェブサイトからGoogleスプレッドシートを作成
https://docs.google.com/spreadsheets/create

名前はプロジェクト名とする。R2D2では、pyディレクトリの上のディレクトリ名を読みそれをスプレッドシートの名前として情報を送るので、ディレクトリと同じ名前にする。

.. image:: figs/google/spread_sheet1.png
    :width: 400 px

講習会ではR2D2としておく。

.. image:: figs/google/spread_sheet2.png
    :width: 400 px

共有をクリックし、ダウンロードしたjsonファイルの中のclient_email行のEメールアドレスをコピーして、貼り付け。ここまでで、R2D2からGoogleスプレッドシートにアクセスできるようになる。

IDLコードの環境設定
----------------------------------------

システムにIDLをインストールすれば、それのみで使える。ここでは説明しない。

最終更新日：|today|