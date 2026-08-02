#include <stdio.h>
#include <math.h>
#include <float.h>

#define MAX_VARIABLES 10  //変数の最大値
#define MAX_CONSTRAINTS 20  //制約条件の式の数の最大値
#define MAX_TABLEAU_COLUMNS 50

#define EPSILON 1e-9

double objective_function[MAX_VARIABLES + 1];  //目的関数
double constraints[MAX_CONSTRAINTS][MAX_VARIABLES + 1];  //制約条件

double phase1_tableau[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS];  //第１段階シンプレックスタブロ

double phase2_tableau[MAX_CONSTRAINTS + 1][MAX_TABLEAU_COLUMNS];  //第２段階シンプレックスタブロ

int basis[MAX_CONSTRAINTS];  //基底変数を保持

int var;  //変数の数
int con;  //制約条件の数

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

void print_phase1_tableau(){  //phase1タブロの表示
    fprintf(stderr, "\n            定数    ");
    for(int i = 1; i <= var + con; i ++){
        fprintf(stderr, "  -x_%d    ", i);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, " z   = ");
    for(int i = 0; i <= var + con; i ++){
        fprintf(stderr," %8.2f ", phase1_tableau[0][i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, " -w  = ");
    for(int i = 0; i <= var + con; i ++){
        fprintf(stderr," %8.2f ", phase1_tableau[1][i]);
    }
    fprintf(stderr, "\n");

    for(int i = 0; i < con; i ++){
        fprintf(stderr, " x_%d = ", basis[i]);
        for(int j = 0; j <= var + con; j ++){
            fprintf(stderr," %8.2f ", phase1_tableau[2 + i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void create_initial_tableau(){
    //配列basisの初期化と決定
    for(int i = 0; i < MAX_CONSTRAINTS; i ++){
        basis[i] = -1;
    }
    for(int i = 0; i < con; i ++){
        basis[i] = var + 1 + i;
    }

    for(int i = 0; i < MAX_CONSTRAINTS + 2; i ++){  //配列phase1_tableauの初期化
        for(int j = 0; j < MAX_TABLEAU_COLUMNS; j ++){
            phase1_tableau[i][j] = NAN;
        }
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
        phase1_tableau[1][i] = 0;
        for(int j = 0; j < con; j ++){
            phase1_tableau[1][i] += phase1_tableau[2 + j][i];
        }
        phase1_tableau[1][i] *= -1;
    }

    print_phase1_tableau();
}

void phase1_pivot(int entering, int leaving){
    double tmp[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS];

    for(int i = 0; i < MAX_CONSTRAINTS + 2; i ++){
        for(int j = 0; j < MAX_TABLEAU_COLUMNS; j ++){
            tmp[i][j] = phase1_tableau[i][j];
        }
    }

    //交換する基底変数の行の設定
    int leaving_line;
    for(int i = 0; i < con; i ++){
        if(basis[i] == leaving){
            leaving_line = 2 + i;
        }
    }

    basis[leaving_line - 2] = entering;

    for(int i = 0; i <= var + con; i ++){
        phase1_tableau[leaving_line][i] /= tmp[leaving_line][entering];
    }
    phase1_tableau[leaving_line][leaving] = 1 / tmp[leaving_line][entering];
    phase1_tableau[leaving_line][entering] = NAN;

    //残りの規定変数の行の操作
    for(int i = 0; i < con; i ++){
        if(basis[i] != entering){
            for(int j = 0; j <= var + con; j ++){
                phase1_tableau[2 + i][j] -= phase1_tableau[leaving_line][j] * tmp[2 + i][entering];
                if(j == leaving){
                    phase1_tableau[2 + i][j] = - 1 / tmp[2 + i][entering];
                }
            }
        }
    }

    //zの行と-wの行の操作
    for(int j = 0; j <= var + con; j ++){
        phase1_tableau[0][j] -= phase1_tableau[leaving_line][j] * tmp[0][entering];
        phase1_tableau[1][j] -= phase1_tableau[leaving_line][j] * tmp[1][entering];
        if(j == leaving){
            phase1_tableau[0][j] = - tmp[0][entering] * phase1_tableau[leaving_line][leaving];
            phase1_tableau[1][j] = - tmp[1][entering] * phase1_tableau[leaving_line][leaving];
        }
    }

    print_phase1_tableau();
}

void phase1_determine_entering_and_leaving_variable(){
    int entering_variable;  //交換する非基底変数
    int leaving_variable;  //交換する基底変数

    //entering_variableの決定
    entering_variable = -1;
    for(int i= 1; i <= var + con; i ++){
        for(int j = 0; j < con; j ++){
            if(i == basis[j]){
                break;
            }
            if(j == con -1){
                entering_variable = i;
            }
        }
        if(entering_variable != -1){
            break;
        }
    }

    for(int i = 1; i <= var + con; i ++){
        if(!isnan(phase1_tableau[1][i]) && phase1_tableau[1][entering_variable] > phase1_tableau[1][i]){
            entering_variable = i;
        }
    }

    //leaving_variableの決定
    double min = phase1_tableau[2][0] / phase1_tableau[2][entering_variable];
    int min_num = 2;

    for(int i = 0; i < con; i ++){  //minの初期値
        if(phase1_tableau[2 + i][entering_variable] > 0){
            min = phase1_tableau[2 + i][0] / phase1_tableau[2 + i][entering_variable];
            min_num = 2 + i;
            break;
        }
    }

    for(int i = 0; i < con; i ++){
        if(phase1_tableau[2 + i][entering_variable] > 0 && min > phase1_tableau[2 + i][0] / phase1_tableau[2 + i][entering_variable]){  //最小添え字規則
            min = phase1_tableau[2 + i][0] / phase1_tableau[2 + i][entering_variable];
            min_num = 2 + i;
        }
    }
    leaving_variable = basis[min_num  - 2];

    fprintf(stderr, "\nx_%dとx_%dの交換\n", entering_variable, leaving_variable);

    phase1_pivot(entering_variable, leaving_variable);
}

int phase1(){
    fprintf(stderr, "\nphase1\n");
    double previous_value_of_negative_w = phase1_tableau[1][0];  //前回の-wの目的関数値

    while(fabs(previous_value_of_negative_w) > EPSILON){
        phase1_determine_entering_and_leaving_variable();
        if(previous_value_of_negative_w > phase1_tableau[1][0]){
            fprintf(stderr, "\n実行不可能です\n");
            return 1;
        }
        previous_value_of_negative_w = phase1_tableau[1][0];  //更新
    }
    fprintf(stderr, "\n実行可能です\n");
    return 0;
}

void print_phase2_tableau(){  //phase2タブロの表示
    fprintf(stderr, "\n            定数    ");
    for(int i = 1; i <= var; i ++){
        fprintf(stderr, "  -x_%d    ", i);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, " z   = ");
    for(int i = 0; i <= var; i ++){
        fprintf(stderr," %8.2f ", phase2_tableau[0][i]);
    }
    fprintf(stderr, "\n");

    for(int i = 0; i < con; i ++){
        fprintf(stderr, " x_%d = ", basis[i]);
        for(int j = 0; j <= var; j ++){
            fprintf(stderr," %8.2f ", phase2_tableau[1 + i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void phase2_pivot(int entering, int leaving){
    double tmp[MAX_CONSTRAINTS + 1][MAX_TABLEAU_COLUMNS];

    for(int i = 0; i < MAX_CONSTRAINTS + 1; i ++){
        for(int j = 0; j < MAX_TABLEAU_COLUMNS; j ++){
            tmp[i][j] = phase2_tableau[i][j];
        }
    }

    //交換する基底変数の行の設定
    int leaving_line;
    for(int i = 0; i < con; i ++){
        if(basis[i] == leaving){
            leaving_line = 1 + i;
        }
    }

    basis[leaving_line - 1] = entering;

    for(int i = 0; i <= var; i ++){
        phase2_tableau[leaving_line][i] /= tmp[leaving_line][entering];
    }
    phase2_tableau[leaving_line][leaving] = 1 / tmp[leaving_line][entering];
    phase2_tableau[leaving_line][entering] = NAN;

    //残りの規定変数の行の操作
    for(int i = 0; i < con; i ++){
        if(basis[i] != entering){
            for(int j = 0; j <= var; j ++){
                phase2_tableau[1 + i][j] -= phase2_tableau[leaving_line][j] * tmp[1 + i][entering];
                if(j == leaving){
                    phase2_tableau[1 + i][j] = - phase2_tableau[leaving_line][j] * tmp[1 + i][entering];
                }
            }
        }
    }

    //zの行の操作
    for(int j = 0; j <= var; j ++){
        phase2_tableau[0][j] -= phase2_tableau[leaving_line][j] * tmp[0][entering];
        if(j == leaving){
            phase2_tableau[0][j] = - tmp[0][entering] * phase2_tableau[leaving_line][leaving];
        }
    }

    print_phase2_tableau();
}

void phase2_determine_entering_and_leaving_variable(){
    int entering_variable;  //交換する非基底変数
    int leaving_variable;  //交換する基底変数

    //entering_variableの決定
    entering_variable = -1;
    for(int i= 1; i <= var; i ++){
        for(int j = 0; j < con; j ++){
            if(i == basis[j]){
                break;
            }
            if(j == con -1){
                entering_variable = i;
            }
        }
        if(entering_variable != -1){
            break;
        }
    }

    for(int i = 1; i <= var; i ++){
        if(!isnan(phase2_tableau[0][i]) && phase2_tableau[0][entering_variable] > phase2_tableau[0][i]){
            entering_variable = i;
        }
    }

    //leaving_variableの決定
    double min;
    int min_num;
    for(int i = 0; i < con; i ++){  //minの初期値
        if(phase2_tableau[1 + i][entering_variable] > 0){
            min = phase2_tableau[1 + i][0] / phase2_tableau[1 + i][entering_variable];
            min_num = 1 + i;
            break;
        }
    }

    for(int i = 0; i < con; i ++){
        if(phase2_tableau[1 + i][entering_variable] > 0 && min > phase2_tableau[1 + i][0] / phase2_tableau[1 + i][entering_variable]){  //最小添え字規則
            min = phase2_tableau[1 + i][0] / phase2_tableau[1 + i][entering_variable];
            min_num = 1 + i;
        }
    }
    leaving_variable = basis[min_num  - 1];

    fprintf(stderr, "\nx_%dとx_%dの交換\n", entering_variable, leaving_variable);

    phase2_pivot(entering_variable, leaving_variable);
}

void print_answer(){
    int flag;
    fprintf(stderr, " 最適解は \n");
    for(int i = 1; i <= var; i ++){
        flag = 0;
        for(int j = 0; j < con; j ++){
            if(i == basis[j]){
                fprintf(stderr, " x_%d= %lf ", i, phase2_tableau[1 + j][0]);
                flag = 1;
                break;
            }
        }
        if(flag == 0){
            fprintf(stderr, " x_%d= %d ", i, 0);
        }
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "最適値は %lf \n", phase2_tableau[0][0]);
}

int has_negative_coefficient(){
    for(int i = 1; i <= var; i++){
        if(phase2_tableau[0][i] < 0){
            return 1;
        }
    }
    return 0;
}

void phase2(){
    fprintf(stderr, "\nphase2\n");

    for(int i = 0; i < MAX_CONSTRAINTS + 1; i ++){  //配列phase2_tableauの初期化
        for(int j = 0; j < MAX_TABLEAU_COLUMNS; j ++){
            phase2_tableau[i][j] = NAN;
        }
    }

    //phase2_tableauの設定
    for(int i = 0; i <= var; i ++){
        phase2_tableau[0][i] = phase1_tableau[0][i];
    }
    for(int i = 0; i < con; i ++){
        for(int j = 0; j <= var; j ++){
            phase2_tableau[1 + i][j] = phase1_tableau[2 + i][j];
        }
    }
    print_phase2_tableau();

    while(has_negative_coefficient()){
        phase2_determine_entering_and_leaving_variable();
    }

    print_answer();
}

int main(void){

    in_data();  //データの入力

    create_initial_tableau();  //初期タブロの作成

    if(phase1() == 0){
        phase2();  //フェーズ2
    }


    return 0;
}