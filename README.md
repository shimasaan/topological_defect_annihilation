# topological_defect_annihilation
 整数トポロジカル欠陥の対消滅スクリプト

## defect.mについて
このスクリプトは、±1トポロジカル欠陥の対消滅を可視化する為に書かれたものです。
（より詳しく言うと、垂直配向液晶セルのFréedericksz転移後の、トポロジカル欠陥の緩和ダイナミクスを説明する為。某大学学生実験で使用。）
### 初期条件・境界条件
±1欠陥がx軸上に並ぶ状態を初期条件として用意。それ以外に欠陥が生じないよう、配向場を連続的に接続。周期境界条件を設定。
### 相互作用
簡単のため、各サイトは最近接相互作用で結ばれている。polar相互作用の場合には±1欠陥が対消滅。nematic相互作用の場合、半整数欠陥が対生成される。