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

## 解決済み(続き)

### 2. データ消失リスク: `Data.zip_time` (`pyR2D2/data.py:130-150`) — 2026-07-27 修正済み

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
- **対応内容**: 削除対象を`pending_removal`リストに集約し、`os.remove()`は`with`ブロックを抜けて(zipの目次が確定して)から実行するように変更。書き込みと削除のタイミングを分離した。一時ディレクトリを使った手動テスト(新規書き込み+削除、既存エントリのスキップ、2回目呼び出しでの追記)で動作確認済み。既存のpytestスイート(114件)も全てパス。

## レビューの結果、バグではないと判断したもの

### `check()`がrankファイル欠損時に即`True`を返す (`pyR2D2/data_io/read.py:818-823`, 同様のパターンが`Slice.check`約1969行目, `_BasePrevAftr.check`約2503行目にもあり)

```python
for np0 in range(self.npe):
    ...
    if not filepath.exists():
        print(f"File {filepath} does not exist. Anyway you can delete it.")
        return True   # 他のrankのファイルは1つも検証せずに「削除OK」を返す
```

- 当初は「未検証のまま安全側と誤判定しているのでは」という懸念でレビュー担当エージェントが報告した項目。
- Hottaさんに確認したところ、**これは意図通りの挙動**とのこと: 複数MPIランクのうち一部のファイルだけが欠けている状態は、削除処理が途中で中断された場合以外に起こり得ない。そのため中途半端な状態を残さないよう、検証をスキップしてでも残りのランクファイルを削除しきる必要がある、という設計。
- 対応不要。

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

## 解決済み(続き2) — 2026-07-28、未検証項目をすべて検証・修正

エージェント報告のみだった6件を実際にコードを読んで確認したところ、すべて実在するバグだった(`sync.py`は報告では6メソッドとされていたが、実際には8メソッド該当)。防御的なチェック・クランプ・打ち切りロジックを追加し、`python setup.py clean --all` → `build_ext --inplace`でのリビルド成功、既存pytest 114件全パス、さらに各修正について個別の手動スモークテスト(正常系・異常系とも)で動作を確認済み。

### C++側

- **`pyR2D2/cpp_util/field_line.hpp` (`interpolate3d` / `trace_field_line`)**: 磁力線がグリッド外に出ても`n_steps`まで追跡を続け、範囲外インデックスで読み取っていた。`trace_field_line`にグリッド範囲判定(`in_domain`)を追加し、磁力線が領域を出た時点でトレースを打ち切るように変更。開始点が最初から範囲外なら`RuntimeError`。戻り値の配列は実際に有効な長さに切り詰められる(＝`n_steps`より短くなり得る、という戻り値の仕様変更を伴う)。手動テストで「領域内に留まる短いトレースは指定`n_steps`分そのまま返る」「境界に到達すると打ち切られる」「開始点が範囲外だとエラーになる」の3パターンを確認。
- **`pyR2D2/cpp_util/eos.hpp` (`EOS`コンストラクタ / `eval`配列版)**: テーブル要素数が2未満、テーブル間の形状不一致、`ro_val_np`/`se_val_np`のサイズ不一致、をそれぞれ`RuntimeError`で検出するように変更。正常系(2次元テーブルでの構築・評価)は従来通り動作することを確認。
- **`pyR2D2/cpp_util/yin_yang_convert.hpp`**: `lagrange_interpolation_3rd`内のインデックス計算(`ic`, `jc`)を、`double`のままクランプしてから`size_t`にキャストするよう変更(Yin側・Yang側どちらの呼び出しも保護)。`convert_scalar`の`i_size`/`j_size`計算でも、`qq_yin.i_size`がmargin(既定2)の2倍以下だと符号なし整数演算でアンダーフローするため、事前に`RuntimeError`を投げるチェックを追加。
- **`pyR2D2/cpp_util/rte.hpp` (`eval_tau` / `vertical_upward_rte`)**: `x_np`の長さが`ro_np`の第1軸サイズと一致しない場合に`RuntimeError`を投げるチェックを両関数に追加。

### Python側

- **`pyR2D2/sync/sync.py`**: `project=os.getcwd().split("/")[-2]`という、import時に1回だけ評価されるデフォルト引数が8メソッド(`setup`, `tau`, `remap_qq`, `xselect`, `vc`, `check`, `slice`, `all`)にあった。全て`project=None`に変更し、関数本体の先頭で`if project is None: project = os.getcwd().split("/")[-2]`と解決するように修正。呼び出し時点のカレントディレクトリが反映されることを確認。

## 優先度メモ

1. ~~`rte.hpp:152` (メモリ破壊)~~ — 2026-07-27 修正済み
2. ~~`data.py`の`zip_time` (データ消失)~~ — 2026-07-27 修正済み
3. ~~`read.py`の`check()`~~ — 2026-07-27 レビューの結果、意図通りの挙動と判明。対応不要
4. `rte.hpp`のオフバイワン疑い — 判断待ち(唯一の残項目)
5. ~~未検証項目6件~~ — 2026-07-28 すべて検証・修正済み

## データ定義の変更の周知 — 2026-08-25 追記

### `vl_spex.dac` の degree スペクトル (`fvxl`〜`fsel`) の定義が R2D2plus (C++) で変わった

R2D2plus 側の修正 DEC-258 (2026-08-25) で、球面調和関数スペクトル出力の
**degree スペクトル (`fvxl, fvyl, fvzl, fbxl, fbyl, fbzl, fsel`) の
$m \ne 0$ 成分のパワー2倍過剰**が修正された。

- 旧 Fortran (`-DSPEX`): degree = $|a_0|^2 + 4S$ ($S$ は $m \ne 0$ のパワー和。
  正負の振動数を「和」で畳んでいたため 2 倍過剰)
- R2D2plus (修正後): degree = $|a_0|^2 + 2S$ (正しい定義。独立参照実装
  sht_py と float32 精度 7.3e-8 で一致することを確認済み)
- order スペクトル (`fvxm`〜`fsem`) = $|a_0|^2$ は**両者で不変**

**換算式** (仮定なしの厳密な関係):

```text
新degree = (旧degree + 旧order) / 2
```

pyR2D2 は `vl_spex.dac` の生の値をそのまま読むだけなので**コード変更は不要**。
ただし旧 Fortran の出力と R2D2plus の出力を混ぜて比較する場合は、上の換算を
挟むこと。なお旧 Fortran では `-DSPEX` は本番で一度も有効化されておらず
(しかも立てないと `vl_spex.dac` の中身は未初期化のゴミ)、この定義変更に
影響される既存データ・解析は無いはずである。

あわせて DEC-258 では「YinYang では `xyz.dac` を書かない」(旧方針) も撤回
された。pyR2D2 の `Parameters` 初期化が `xyz.dac` を無条件に開くため、
R2D2plus の YinYang 出力もネイティブに開けることを確認済み。
