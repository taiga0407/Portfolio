#include <stdio.h>

#define MAX_VARIABLES 10  //変数の最大値
#define MAX_CONSTRAINTS 20  //制約条件の式の数の最大値
#define MAX_TABLEAU_COLUMNS 50

double objective_function[MAX_VARIABLES + 1];  //目的関数
double constraints[MAX_CONSTRAINTS][MAX_VARIABLES + 1];  //制約条件

double phase1_tableau[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS] ;  //第１段階シンプレックスタブロ

int basis[MAX_CONSTRAINTS];  //基底変数を保持

int var;  //変数の数
int con;  //制約条件の数

// void print_number(double value)
// {
//     if (value == (int)value) {
//         printf("%d", (int)value);
//     } else {
//         printf("%f", value);
//     }
// }

void check_data(){
    fprintf(stderr, "目的関数  z = ");
    for(int i = 1; i <= var; i ++){
        fprintf(stderr, "%2fx_%d", objective_function[i], i);
        fprintf(stderr, " + ");
    }
    fprintf(stderr, "%f", objective_function[0]);

    fprintf(stderr, "  (最大化)\n");

    fprintf(stderr, "制約条件\n");
    for(int i = 0; i < con; i ++){
        for(int j = 1; j <= var; j ++){
            fprintf(stderr, "%2fx_%d", constraints[i][j], j);
                if(j != var){
                    fprintf(stderr, " + ");
            }
        }
        fprintf(stderr, " = %2f\n", constraints[i][0]);
    }

    fprintf(stderr, " x_i>= 0  (i = 1,2,...,%d)\n", var);

}

void in_data(){
    fprintf(stderr, "変数の数を入力してください．-> ");
    scanf("%d", &var);

    fprintf(stderr, "目的関数の係数の入力\n");
    for(int i = 1; i <= var; i ++){
        fprintf(stderr, "x_%dの係数：", i);
        scanf("%lf", &objective_function[i]);
    }
    fprintf(stderr, "定数：");
    scanf("%lf", &objective_function[0]);

    fprintf(stderr, "制約条件の式の数の数を入力してください．-> ");
    scanf("%d", &con);

    fprintf(stderr, "制約条件の入力\n");
    for(int i = 0; i < con; i ++){
            fprintf(stderr, "制約条件 %d\n", i);
        for(int j = 1; j <= var; j ++){
            fprintf(stderr, "x_%dの係数：", j);
            scanf("%lf", &constraints[i][j]);
        }
        fprintf(stderr, "右辺の定数：");
        scanf("%lf", &constraints[i][0]);
    }
    check_data();
}

void print_tableau(){  //タブロの表示
    fprintf(stderr, "\n            定数    ");
    for(int i = 1; i <= var; i ++){
        fprintf(stderr, "  -x_%d    ", i);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, " z   = ");
    for(int i = 0; i <= var; i ++){
        fprintf(stderr," %8.2f ", phase1_tableau[0][i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, " -w  = ");
    for(int i = 0; i <= var; i ++){
        fprintf(stderr," %8.2f ", phase1_tableau[1][i]);
    }
    fprintf(stderr, "\n");

    for(int i = 0; i < con; i ++){
        fprintf(stderr, " x_%d = ", basis[i]);
        for(int j = 0; j <= var; j ++){
            fprintf(stderr," %8.2f ", phase1_tableau[2 + i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void creat_inintial_tableau(){
    //配列basisの初期化と決定
    for(int i = 0; i < MAX_CONSTRAINTS; i ++){
        basis[i] = -1;
    }
    for(int i = 0; i < con; i ++){
        basis[i] = var + 1 + i;
    }

    //元の目的関数z
    for(int i = 0; i <= var; i ++){
        phase1_tableau[0][i] = -objective_function[i];
    }
    //人為変数x_var+1からx_var+conまで
    for(int i = 0; i < con; i ++){
        for(int j = 0; j <= var; j ++){
            phase1_tableau[2 + i][j] = constraints[i][j];
        }
    }
    //人為変数導入後の目的関数-w
    for(int i = 0; i <= var; i ++){
        for(int j = 0; j < con; j ++){
            phase1_tableau[1][i] += phase1_tableau[2 + j][i];
        }
        phase1_tableau[1][i] *= -1;
    }

    print_tableau();
}

int main(void){

    in_data();  //データの入力

    creat_inintial_tableau();  //初期タブロの作成

    return 0;
}