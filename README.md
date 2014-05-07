# bed_mod.py
---

author  : @shibasay 

date    : 2014/5/5

---

## abstruct 
@aki426 による3Dモデリングソフト be! のセーブファイル（.bed）を読み込み、モデルを修正するスクリプトです。

## usage
`python ./bed_mod.py -i [filename.bed] -o [outname.bed] -m [mode] -l [level] -r [mirrormode]`

## 引数説明

* [filename.bed]: 入力ファイル名。これはbe!で作ってね
* [outname.bed] : 出力ファイル名
* [mode]        : 処理モード 0 or 1
    * 0 中空モード（hollow）
    * 1 埋めモード（fill）
	* 2 鏡像モード (mirror)
    * デフォルト値は 0
* [level]       : 中空にするレベル。mode 0 の場合にのみ使用。
     * 0 何も消えません。
     * 1 線接触
     * 2以上 面接触
     * デフォルト値は１
* [mirrormode]  : 鏡像モードのモード。mode 2 の場合にのみ使用。
     * 0 左側を使った鏡像
     * 1 右側を使った鏡像
     * 2 もう一個隣に作る胸像
     * デフォルト値は0

## アルゴリズム詳細（へぼい）
### 中空モード
内部（ボクセルにより閉じた領域内）にあるボクセルを消します。閉じた空間であるかどうか（厳密には開いた空間であるかどうか）は、以下のアルゴリズムで判定します（厳密にはソースで描いていることと違ったりするかもしれませんが、アイデアレベルでは同じはず）。

1. bedファイルで定義されたボクセル群を包含する直方体を定義する
2. 直方体の、x,y,zすべてが最小値になる頂点に、外側を示す仮想ボクセルをつくる
3. 外側を示す仮想ボクセルから、x,y,zをそれぞれ1足した場所（計三箇所）を見る。
4. 見た箇所にボクセルがなければ、まだ外側とみなし、そこを中心として再帰的に外側探索を行う。ボクセルがあれば内側とみなして再帰的探索しない。
5. 再帰処理が終わった後、外側を示す仮想ボクセルの無い部分は、ボクセルの有無にかかわらず、外側ではないとみなせる
6. 外側と接する箇所をlevel0、それより内側を順にlevel1、level2とし、中空レベル（level）にあわせて削除する
7. 外側を示す仮想ボクセルを削除する


### 埋めモード
中空モードと同様の閉じた空間判定を行い、閉じた空間を固定ボクセルで埋めます。色はハードコードされており、現在は白（r,g,b=255,255,255）です。
### 鏡像モード
左右反転したモデルを生成します。
作ってみたけど左右非対称になっちゃった…。というときに使用することを想定して作りました。
副作用として、左右反転したモデルを隣にもう一個作るオプションもできました。
* 0: モデルの向かって右側半分を消し、左側を反転させたものを代わりに埋めます。左側を使った鏡像になります。
* 1: モデルの向かって左側半分を消し、右側を反転させたものを代わりに埋めます。右側を使った鏡像になります。
* 2: モデルを丸ごと左右反転させたものを、向かってすぐ右隣に配置します。モデルが簡単に2個になるよ！
