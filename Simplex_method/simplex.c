#include <stdio.h>
#include <math.h>
#include <float.h>

#define MAX_VARIABLES 10  //変数の最大値
#define MAX_CONSTRAINTS 20  //制約条件の式の数の最大値
#define MAX_TABLEAU_COLUMNS 50

#define EPSILON 1e-9

typedef struct {
    double tableau[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS];  //シンプレックスタブロ
    int basis[MAX_CONSTRAINTS];  //基底変数を保持
    int variable_count;  //変数の数
    int constraint_count;  //制約条件の数
} SimplexTableau_t;

void check_data(SimplexTableau_t simplex, double objective_function[MAX_VARIABLES + 1], double constraints[MAX_CONSTRAINTS][MAX_VARIABLES + 1]){
    printf( "目的関数\nz = ");
    for(int i = 1; i <= simplex.variable_count; i ++){
        printf("%2fx_%d", objective_function[i], i);
        printf(" + ");
    }
    printf("%f", objective_function[0]);

    printf("  (最大化)\n");

    printf("\n制約条件\n");
    for(int i = 0; i < simplex.constraint_count; i ++){
        for(int j = 1; j <= simplex.variable_count; j ++){
            printf("%2fx_%d", constraints[i][j], j);
                if(j != simplex.variable_count){
                    printf(" + ");
            }
        }
        printf(" = %2f\n", constraints[i][0]);
    }

    printf("x_i>= 0  (i = 1,2,...,%d)\n", simplex.variable_count);
}

void in_data(SimplexTableau_t *simplex, double objective_function[MAX_VARIABLES + 1], double constraints[MAX_CONSTRAINTS][MAX_VARIABLES + 1]){
    fprintf(stderr, "変数の数を入力してください．-> ");
    scanf("%d", &simplex->variable_count);

    fprintf(stderr, "目的関数の係数の入力\n");
    for(int i = 1; i <= simplex->variable_count; i ++){
        fprintf(stderr, "x_%dの係数：", i);
        scanf("%lf", &objective_function[i]);
    }
    fprintf(stderr, "定数：");
    scanf("%lf", &objective_function[0]);

    fprintf(stderr, "制約条件の式の数の数を入力してください．-> ");
    scanf("%d", &simplex->constraint_count);

    fprintf(stderr, "制約条件の入力\n");
    for(int i = 0; i < simplex->constraint_count; i ++){
            fprintf(stderr, "制約条件 %d\n", i);
        for(int j = 1; j <= simplex->variable_count; j ++){
            fprintf(stderr, "x_%dの係数：", j);
            scanf("%lf", &constraints[i][j]);
        }
        fprintf(stderr, "右辺の定数：");
        scanf("%lf", &constraints[i][0]);
    }
    check_data(*simplex, objective_function, constraints);
}

void print_phase1_tableau(SimplexTableau_t simplex){  //phase1タブロの表示
    fprintf(stderr, "\n            定数    ");
    for(int i = 1; i <= simplex.variable_count + simplex.constraint_count; i ++){
        fprintf(stderr, "  -x_%d    ", i);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, " z   = ");
    for(int i = 0; i <= simplex.variable_count + simplex.constraint_count; i ++){
        fprintf(stderr," %8.2f ", simplex.tableau[0][i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, " -w  = ");
    for(int i = 0; i <= simplex.variable_count + simplex.constraint_count; i ++){
        fprintf(stderr," %8.2f ", simplex.tableau[1][i]);
    }
    fprintf(stderr, "\n");

    for(int i = 0; i < simplex.constraint_count; i ++){
        fprintf(stderr, " x_%d = ", simplex.basis[i]);
        for(int j = 0; j <= simplex.variable_count + simplex.constraint_count; j ++){
            fprintf(stderr," %8.2f ", simplex.tableau[2 + i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void print_phase2_tableau(SimplexTableau_t simplex){  //phase2タブロの表示
    fprintf(stderr, "\n            定数    ");
    for(int i = 1; i <= simplex.variable_count; i ++){
        fprintf(stderr, "  -x_%d    ", i);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, " z   = ");
    for(int i = 0; i <= simplex.variable_count; i ++){
        fprintf(stderr," %8.2f ", simplex.tableau[0][i]);
    }
    fprintf(stderr, "\n");

    for(int i = 0; i < simplex.constraint_count; i ++){
        fprintf(stderr, " x_%d = ", simplex.basis[i]);
        for(int j = 0; j <= simplex.variable_count; j ++){
            fprintf(stderr," %8.2f ", simplex.tableau[2 + i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void create_initial_tableau(SimplexTableau_t *simplex, double objective_function[MAX_VARIABLES + 1], double constraints[MAX_CONSTRAINTS][MAX_VARIABLES + 1]){
    //配列basisの初期化と決定
    for(int i = 0; i < MAX_CONSTRAINTS; i ++){
        simplex->basis[i] = -1;
    }
    for(int i = 0; i < simplex->constraint_count; i ++){
        simplex->basis[i] = simplex->variable_count + 1 + i;
    }

    for(int i = 0; i < MAX_CONSTRAINTS + 2; i ++){  //配列phase1_tableauの初期化
        for(int j = 0; j < MAX_TABLEAU_COLUMNS; j ++){
            simplex->tableau[i][j] = NAN;
        }
    }

    //元の目的関数z
    for(int i = 0; i <= simplex->variable_count; i ++){
        simplex->tableau[0][i] = -objective_function[i];
    }
    //人為変数x_variable_count+1からx_variable_count+constraint_countまで
    for(int i = 0; i < simplex->constraint_count; i ++){
        for(int j = 0; j <= simplex->variable_count; j ++){
            simplex->tableau[2 + i][j] = constraints[i][j];
        }
    }
    //人為変数導入後の目的関数-w
    for(int i = 0; i <= simplex->variable_count; i ++){
        simplex->tableau[1][i] = 0;
        for(int j = 0; j < simplex->constraint_count; j ++){
            simplex->tableau[1][i] += simplex->tableau[2 + j][i];
        }
        simplex->tableau[1][i] *= -1;
    }

    print_phase1_tableau(*simplex);
}

void pivot(SimplexTableau_t *simplex, int entering, int leaving, int phase){  //phase1のときphaseにはconstraint_countを格納，phase2のときphaseには0を格納
    double tmp[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS];

    for(int i = 0; i < MAX_CONSTRAINTS + 2; i ++){
        for(int j = 0; j < MAX_TABLEAU_COLUMNS; j ++){
            tmp[i][j] = simplex->tableau[i][j];
        }
    }

    //交換する基底変数の行の設定
    int leaving_line;
    for(int i = 0; i < simplex->constraint_count; i ++){
        if(simplex->basis[i] == leaving){
            leaving_line = 2 + i;
        }
    }

    simplex->basis[leaving_line - 2] = entering;

    for(int i = 0; i <= simplex->variable_count + phase; i ++){
        simplex->tableau[leaving_line][i] /= tmp[leaving_line][entering];
    }
    simplex->tableau[leaving_line][leaving] = 1 / tmp[leaving_line][entering];
    simplex->tableau[leaving_line][entering] = NAN;

    //残りの規定変数の行の操作
    for(int i = 0; i < simplex->constraint_count; i ++){
        if(simplex->basis[i] != entering){
            for(int j = 0; j <= simplex->variable_count + phase; j ++){
                simplex->tableau[2 + i][j] -= simplex->tableau[leaving_line][j] * tmp[2 + i][entering];
                if(j == leaving){
                    simplex->tableau[2 + i][j] = - simplex->tableau[leaving_line][j] * tmp[2 + i][entering];
                }
            }
        }
    }

    //zの行と-wの行の操作
    for(int j = 0; j <= simplex->variable_count + simplex->constraint_count; j ++){
        simplex->tableau[0][j] -= simplex->tableau[leaving_line][j] * tmp[0][entering];
        simplex->tableau[1][j] -= simplex->tableau[leaving_line][j] * tmp[1][entering];
        if(j == leaving){
            simplex->tableau[0][j] = - tmp[0][entering] * simplex->tableau[leaving_line][leaving];
            simplex->tableau[1][j] = - tmp[1][entering] * simplex->tableau[leaving_line][leaving];
        }
    }

    if(phase == simplex->constraint_count){
        print_phase1_tableau(*simplex);
    }else{
        print_phase2_tableau(*simplex);
    }
}

void determine_entering_and_leaving_variable(SimplexTableau_t *simplex, int phase){
    int entering_variable;  //交換する非基底変数
    int leaving_variable;  //交換する基底変数

    //entering_variableの決定
    entering_variable = -1;
    for(int i= 1; i <= simplex->variable_count + phase; i ++){
        for(int j = 0; j < simplex->constraint_count; j ++){
            if(i == simplex->basis[j]){
                break;
            }
            if(j == simplex->constraint_count -1){
                entering_variable = i;
            }
        }
        if(entering_variable != -1){
            break;
        }
    }

    int objective_value;
    if(phase == simplex->constraint_count){
        objective_value = 1;
    }else{
        objective_value = 0;
    }
    for(int i = 1; i <= simplex->variable_count + phase; i ++){
        if(!isnan(simplex->tableau[objective_value][i]) && simplex->tableau[objective_value][entering_variable] > simplex->tableau[objective_value][i]){
            entering_variable = i;
        }
    }

    //leaving_variableの決定
    double min_ratio;
    int leaving_row;

    for(int i = 0; i < simplex->constraint_count; i ++){  //min_ratioの初期値
        if(simplex->tableau[2 + i][entering_variable] > 0){
            min_ratio = simplex->tableau[2 + i][0] / simplex->tableau[2 + i][entering_variable];
            leaving_row = 2 + i;
            break;
        }
    }

    for(int i = 0; i < simplex->constraint_count; i ++){
        if(simplex->tableau[2 + i][entering_variable] > 0 && min_ratio > simplex->tableau[2 + i][0] / simplex->tableau[2 + i][entering_variable]){  //最小添え字規則
            min_ratio = simplex->tableau[2 + i][0] / simplex->tableau[2 + i][entering_variable];
            leaving_row = 2 + i;
        }
    }
    leaving_variable = simplex->basis[leaving_row  - 2];

    fprintf(stderr, "\nx_%dとx_%dの交換\n", entering_variable, leaving_variable);

    if(phase == simplex->constraint_count){
        pivot(simplex, entering_variable, leaving_variable, simplex->constraint_count);
    }else{
        pivot(simplex, entering_variable, leaving_variable, 0);
    }
}

int phase1(SimplexTableau_t *simplex){
    fprintf(stderr, "\nphase1\n");
    double previous_value_of_negative_w = simplex->tableau[1][0];  //前回の-wの目的関数値

    while(fabs(previous_value_of_negative_w) > EPSILON){
        determine_entering_and_leaving_variable(simplex, simplex->constraint_count);
        if(previous_value_of_negative_w > simplex->tableau[1][0]){
            fprintf(stderr, "\n実行不可能です\n");
            return 1;
        }
        previous_value_of_negative_w = simplex->tableau[1][0];  //更新
    }
    fprintf(stderr, "\n実行可能です\n");
    return 0;
}

void print_answer(SimplexTableau_t simplex){
    int flag;
    printf("\n最適解は\n");
    for(int i = 1; i <= simplex.variable_count; i ++){
        flag = 0;
        for(int j = 0; j < simplex.constraint_count; j ++){
            if(i == simplex.basis[j]){
                printf("x_%d= %lf", i, simplex.tableau[2 + j][0]);
                flag = 1;
                break;
            }
        }
        if(flag == 0){
            printf("x_%d= %d", i, 0);
        }
        if(i != simplex.variable_count){
            printf(" ");
        }
    }
    printf("\n");

    printf("最適値は %lf", simplex.tableau[0][0]);
}

int has_negative_coefficient(SimplexTableau_t simplex){
    for(int i = 1; i <= simplex.variable_count; i++){
        if(simplex.tableau[0][i] < 0){
            return 1;
        }
    }
    return 0;
}

void phase2(SimplexTableau_t *simplex){
    fprintf(stderr, "\nphase2\n");

    print_phase2_tableau(*simplex);

    while(has_negative_coefficient(*simplex)){
        determine_entering_and_leaving_variable(simplex, 0);
    }

    print_answer(*simplex);
}

int main(void){
    double objective_function[MAX_VARIABLES + 1];  //目的関数
    double constraints[MAX_CONSTRAINTS][MAX_VARIABLES + 1];  //制約条件

    SimplexTableau_t simplex;

    in_data(&simplex, objective_function, constraints);  //データの入力

    create_initial_tableau(&simplex, objective_function, constraints);  //初期タブロの作成

    if(phase1(&simplex) == 0){
        phase2(&simplex);  //フェーズ2
    }

    return 0;
}