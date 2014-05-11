# bed_mod.py
---

author  : @shibasay 

date    : 2014/5/10 (last update)

---

## abstruct 
@aki426 による3Dモデリングソフト be! のセーブファイル（.bed）を読み込み、モデルを修正するスクリプトです。

## usage
`python ./bed_mod.py -i [filename.bed] -o [outname.bed] -m [mode] -l [level] -r [mirrormode] -s [scale] -x [mvx] -y [mvy] -z [mvz] --xdegreee [xdegree] --ydegree [ydegree] --zdegree [zdegree]`

Pythonのバージョンは、Python2系にのみ対応です。Python3系には対応していません。

CPython2.7.3、および PyPy 2.3.0 にて動作確認済み。

## 引数説明

* [filename.bed]: 入力ファイル名（必須オプション）。これはbe!で作ってね
* [outname.bed] : 出力ファイル名。指定が無い場合、標準出力にbedデータを出力します。
* [mode]        : 処理モード（必須オプション）
    * 0 中空モード（hollow）
    * 1 埋めモード（fill）
	* 2 鏡像モード (mirror)
    * 3 スケールモード (scale)
    * 4 平行移動モード (move)
    * 5 回転モード (rotate)
    * デフォルト値は 0
* [level]       : 中空にするレベル。mode 0 の場合にのみ使用。
     * 0 何も消えません。
     * 1 線接触
     * 2以上 面接触
     * デフォルト値は１
* [mirrormode]  : 鏡像モードのモード。mode 2 の場合にのみ使用。
     * 0 左側を使った鏡像
     * 1 右側を使った鏡像
     * 2 もう一個隣に作る鏡像
     * デフォルト値は0
* [scale] : スケールモードの、スケール率。mode 3 の場合にのみ使用
     * デフォルト値は 1 (no scale)
     * 整数値のみ対応
* [mvx] [mvy] [mvz] : 平行移動モードでの移動量。 mode 4 の場合にのみ使用
     * デフォルト値はどれも0（移動なし）
     * 整数値のみ対応
* [xdegree] [ydegree] [zdegree]: X/Y/Z軸まわりの回転量
     * デフォルト値は 0
     * 整数で、度数（0度、90度など）で指定

## アルゴリズム詳細（へぼい）
### 中空モード
内部（ボクセルにより閉じた領域内）にあるボクセルを消します。閉じた空間であるかどうか（厳密には開いた空間であるかどうか）は、以下のアルゴリズムで判定します（厳密にはソースで描いていることと違ったりするかもしれませんが、アイデアレベルでは同じはず）。

1. bedファイルで定義されたボクセル群を包含する直方体を定義する
2. 直方体の、x,y,zすべてが最小値になる頂点に、外側を示す仮想ボクセルをつくる
3. 外側を示す仮想ボクセルの周辺（x,y,zをそれぞれ+1 or -1 した場所（計６箇所）を見る。
4. 見た箇所にボクセルがなければ、まだ外側とみなし、そこを中心として繰り返し外側探索を行う。ボクセルがあれば内側とみなして再帰的探索しない。
5. 直方体内部で探索していない外側ボクセルがなくなったら探索終了。
6. 外側と接する箇所をlevel0、それより内側を順にlevel1、level2とし、中空レベル（level）にあわせて削除する
7. 外側を示す仮想ボクセルを削除する

### 埋めモード
中空モードと同様の閉じた空間判定を行い、閉じた空間を固定ボクセルで埋めます。色は周囲ボクセルの平均色が使われます。
なお、実装上、もともと周囲に1個もボクセルのない箇所については、白で埋められます。
外側から中側へじわじわと色を伝播させていくアルゴリズムではないのです。

### 鏡像モード
左右反転したモデルを生成します。
作ってみたけど左右非対称になっちゃった…。というときに使用することを想定して作りました。
副作用として、左右反転したモデルを隣にもう一個作るオプションもできました。

* 0: モデルの向かって右側半分を消し、左側を反転させたものを代わりに埋めます。左側を使った鏡像になります。
* 1: モデルの向かって左側半分を消し、右側を反転させたものを代わりに埋めます。右側を使った鏡像になります。
* 2: モデルを丸ごと左右反転させたものを、向かってすぐ右隣に配置します。モデルが簡単に2個になるよ！

### スケールモード
ボクセル数をX/Y/Z軸方向すべてに整数倍するモード。ただボクセル数を増やすだけだと意味が無いので、滑らかにするために角をおとします。
ただしアルゴリズムが現状へぼいので、明らかな角しか落とせません。また、ドット絵の感じからして落とさなくていいところまで落としたりします。
気に入らないときは手打ちで直してください。

### 平行移動モード
なんの難しさも無い平行移動。なぜ最初に実装しなかったのか…。

### 回転モード
回転行列を用いて実装。x/y/z軸周りの回転量をそれぞれ指定できます。
ただし最後の浮動小数点→整数変換の扱いをまじめにしていないため、90の倍数以外の度数に対してはボクセル抜けが発生したりします。
あらかじめご了承ください。

