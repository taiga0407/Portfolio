# Two-Phase Simplex Method in C

C言語で実装した2段階シンプレックス法による線形計画問題の解を求めるプログラムです．

制約条件に等式制約を持つ線形計画問題に対して，人為変数を導入する2段階シンプレックス法を用いて最適解を探索します．

## Overview

本プログラムでは以下の流れで線形計画問題を解きます．

1. 初期タブロの作成
2. Phase 1（第一段階）
   - 人為変数を導入
   - 補助問題を解き，実行可能解の存在を確認
3. Phase 2（第二段階）
   - 元の目的関数に戻す
   - 最適解を探索

対象とする問題形式：

```
最大化

z = c1x1 + c2x2 + ... + cnxn

制約条件

a11x1 + a12x2 + ... + a1nxn = b1
a21x1 + a22x2 + ... + a2nxn = b2
...

xi >= 0
```

## Features

- C言語によるシンプレックス法の実装
- 2段階シンプレックス法（Two-Phase Simplex Method）
- 人為変数の自動導入
- ピボット操作による基底変数交換
- 最小比率検査による離脱変数選択
- 実行可能性判定
- 最適解・最適値の表示
- シンプレックスタブロの途中経過表示

## Algorithm

### Phase 1

通常のシンプレックス法では初期基底解が必要ですが，制約条件によっては初期基底解を作成できない場合があります．
その場合，人為変数を導入して補助問題を作成します．

補助目的関数：

```
-w = -(人工変数の合計)
```

を最大化することで，実行可能解を探索します．

Phase 1終了時に

```
w = 0
```

となれば実行可能解が存在すると判断し，Phase 2へ移行します．

## Phase 2

元の目的関数を用いて通常のシンプレックス法を実行し，最適解を求めます．

目的関数の係数に負の値が存在する間，ピボット操作を繰り返します．

## Implementation Details

### Data Structure

問題情報とシンプレックスタブロを以下の構造体で管理しています．

```c
typedef struct {
    double objective_function[MAX_VARIABLES + 1];
    double constraints[MAX_CONSTRAINTS][MAX_VARIABLES + 1];
    double tableau[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS];
    int basis[MAX_CONSTRAINTS];

    int variable_count;
    int constraint_count;
} SimplexTableau_t;
```

### Main Functions

| Function | Description |
| --- | --- |
| `in_data()` | 問題データの入力 |
| `create_initial_tableau()` | 初期タブロ作成 |
| `phase1()` | 実行可能性判定 |
| `phase2()` | 最適解探索 |
| `pivot()` | ピボット操作 |
| `determine_entering_variable()` | 進入変数の決定 |
| `determine_leaving_variable()` | 離脱変数の決定 |
| `print_answer()` | 最適解の表示 |

## Environment

### Language

- C

### Compiler

動作確認環境：

- GCC

## Compile

以下のコマンドでコンパイルできます．

```bash
gcc simplex.c
```

## Execute

```bash
./simplex < in1.txt
```

実行すると，変数数，目的関数，制約条件を入力できます．

入力例：

```
変数の数を入力してください．-> 5

目的関数の係数の入力
x_1の係数：-2
x_2の係数：3
x_3の係数：-6
x_4の係数：-1
x_5の係数：2
定数:0

制約条件の式の数を入力してください．-> 2

制約条件 1
x_1の係数：-2
x_2の係数：3
x_1の係数：1
x_2の係数：3
x_1の係数：-1
右辺の定数：3

制約条件 2
x_1の係数：1
x_2の係数：1
x_1の係数：-2
x_2の係数：9
x_1の係数：0
右辺の定数：4
```

## Example Result

問題：

```
maximize

z = -2x_1 + 3x_2 + -6x_3 + -1x_4 + 2x_5 + 0  (最大化)

subject to

2x_1 + -3x_2 + 1x_3 + 3x_4 + -1x_5 = 3
1x_1 + 1x_2 + -2x_3 + 9x_4 + 0x_5 = 4
x_i>= 0  (i = 1,2,...,5)

x1,x2,x3,x4,x5 >= 0
```

実行結果：

```
最適解は
x_1= 4.000000 x_2= 0 x_3= 0 x_4= 0 x_5= 5.000000
最適値は 2.000000
```

## Limitations

現在の実装では以下の制限があります．

- 最大変数数：10
- 最大制約数：20
- 等式制約を対象

## Future Improvements

今後の改善予定：

- 不等式制約（<=, >=）への対応
- 無限解・非有界問題の検出
- 数値誤差へのさらなる対応

## Author

C言語によるアルゴリズム実装学習の一環として作成しました．