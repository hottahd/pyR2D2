# pyR2D2 既知の重大な問題 (2026-07-27時点)

Codexによるsetup.py/パッケージング整備をきっかけに、コードベース全体を軽くレビューして見つかった問題をまとめたもの。
「検証済み」は実際にファイルを開いて内容を確認したもの、「未検証」はレビュー担当エージェントの報告のみで、こちらではまだファイルを開いて確認していないもの。

## 解決済み

### 1. 【最重要】メモリ破壊: `pyR2D2/cpp_util/rte.hpp:152` — 2026-07-27 修正済み

```cpp
rt(j, k) = rt1d[ro.i_size - 2];
```

- `rt_np`は3次元(`i_size, j_size, k_size`)で確保しているが、ここは2引数の`rt(j, k)`を呼んでいる。
- `view_array`の2引数版は`data[i * j_size + j]`という2次元配列用の式を使うため、3次元バッファの想定外の場所を読み書きする。
- **発生条件**: `j_size`が`k_size`より大きいグリッド(例: 2次元スライスで`k_size=1`)。確保領域外への書き込みが起こり得る。
- **症状**: クラッシュせずに他の変数の値が壊れる、サイレントな計算結果の破損。
- **対応内容**: `rt_np`の宣言を`{ro.i_size, ro.j_size, ro.k_size}`から`{ro.j_size, ro.k_size}`(2次元)に変更。`vertical_upward_rte`のdocstringも「2D array」と明記されており、`view_array`は次元数を可変に扱える設計のため、この修正で意図通りの2次元配置になった。`python setup.py clean --all` → `python setup.py build_ext --inplace`でリビルドし、コンパイル成功・importも正常に確認済み。

## 検証済み・未修正

### 2. データ消失リスク: `Data.zip_time` (`pyR2D2/data.py:130-150`)

```python
with zipfile.ZipFile(dst, "a", ...) as zf:
    for p in src.rglob("*"):
        ...
        zf.write(p, arcname)
        if remove_original:
            os.remove(p)   # withブロックを抜ける(=zipの目次を確定させる)前に元ファイルを消している
```

- zipファイルの目次(central directory)は`with`ブロックを抜けて`close()`されるまで確定しない。
- **発生条件**: 大量ファイルをアーカイブ中に、HPCジョブのタイムアウトやOOM Killなどでプロセスが強制終了。
- **症状**: zipが壊れて読めない状態のまま、元ファイルはすでに削除済み。再生成不可能なシミュレーションデータの消失。
- **直し方の方向性**: 全ファイルの`zf.write()`が終わって`with`ブロックを抜けた後にまとめて`os.remove()`する(書き込みと削除のループを分離する)。

### 3. データ消失リスク: `check()`が「未確認」を「安全」と誤判定 (`pyR2D2/data_io/read.py:818-823`, 同様のパターンが`Slice.check`約1969行目, `_BasePrevAftr.check`約2503行目にもあり)

```python
for np0 in range(self.npe):
    ...
    if not filepath.exists():
        print(f"File {filepath} does not exist. Anyway you can delete it.")
        return True   # 他のrankのファイルは1つも検証せずに「削除OK」を返す
```

- 複数MPIランクのうち1つでもファイルが欠けていると、残り全ランクのzarrとの整合性チェックを一切せずに`True`(削除して良い)を返す。
- **発生条件**: 過去の途中終了や部分的なクリーンアップなどで、あるランクのファイルだけ先に無くなっている状態。
- **症状**: `delete(force=False)`がこの`True`を信じて呼ばれると、実際には未検証の他ランクのファイルまで一緒に削除されてしまう。
- **直し方の方向性**: ファイルが存在しないrankをスキップしつつ、存在する残りのrankについては通常通り検証を継続する。全rank欠損時のみ「削除可」とするなど、意図を明確にする。

## 保留中(判断待ち)

### オフバイワンの疑い: `pyR2D2/cpp_util/rte.hpp:132-152`

```cpp
rt1d[1] = so_mid[1];
for (size_t i = 1; i < ro.i_size - 1; i++)
{
    ...
    rt1d[i + 1] = ...;
}
rt(j, k) = rt1d[ro.i_size - 2];
```

- ループは`i`が`ro.i_size - 2`まで進み、最後に書き込まれるのは`rt1d[ro.i_size - 1]`のはず。しかし取り出しているのは`rt1d[ro.i_size - 2]`で、最後の1つ手前の値になっている。
- 意図通り(最上層を境界条件として別扱いにしたいなど)なのか、単なるオフバイワンのミスなのか、物理的な意味が分からないと判断できないため保留。Hottaさんの判断待ち。

## 未検証(エージェント報告のみ、要確認)

### C++側

- `pyR2D2/cpp_util/field_line.hpp:24` (`interpolate3d`): 磁力線追跡でグリッド外に出た座標を範囲チェックなしにインデックスとして使っている疑い。磁力線がグリッド外に出るのは通常の終了条件のため、通常の呼び出しで発生し得るとの報告。
- `pyR2D2/cpp_util/eos.hpp:48`: EOSテーブルの要素数が2未満の場合に境界外読み取りが起きる疑い。テーブル間のサイズ不一致も未チェック。
- `pyR2D2/cpp_util/eos.hpp:124`: `eval`の配列版で`ro_val_np`/`se_val_np`のサイズ一致を未チェック。
- `pyR2D2/cpp_util/yin_yang_convert.hpp:130`: Yang格子側の変換で座標の範囲チェックがYin側と非対称(Yin側はチェックあり、Yang側はなし)。負値からの`size_t`アンダーフローの疑いも。
- `pyR2D2/cpp_util/rte.hpp:56`: `eval_tau`で`x_np`の長さと`ro`の第1軸サイズの一致を未チェック。

### Python側

- `pyR2D2/sync/sync.py:49`ほか計6メソッド:
  ```python
  def some_method(self, ..., project=os.getcwd().split("/")[-2]):
  ```
  デフォルト引数が**import時に1回だけ**評価される。別ディレクトリから`project=`を省略して呼ぶと、古いプロジェクト名のまま同期されてしまう。

## 優先度メモ

1. ~~`rte.hpp:152` (メモリ破壊)~~ — 2026-07-27 修正済み
2. `data.py`の`zip_time` (データ消失) — 優先
3. `read.py`の`check()` (データ消失) — 優先
4. `rte.hpp`のオフバイワン疑い — 判断待ち
5. その他の未検証項目 — 時間があれば個別に確認・修正
