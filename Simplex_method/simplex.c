#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

#define MAX_VARIABLES 10  //変数の最大値
#define MAX_CONSTRAINTS 20  //制約条件の式の数の最大値
#define MAX_TABLEAU_COLUMNS 50

#define Z_ROW 0
#define W_ROW 1
#define CONSTRAINT_ROW 2
#define RHS_COLUMN 0

#define PHASE1 1
#define PHASE2 2

#define EPSILON 1e-9

typedef struct {
    double objective_function[MAX_VARIABLES + 1];  //目的関数
    double constraints[MAX_CONSTRAINTS][MAX_VARIABLES + 1];  //制約条件
    double tableau[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS];  //シンプレックスタブロ
    int basis[MAX_CONSTRAINTS];  //基底変数を保持
    int variable_count;  //変数の数
    int constraint_count;  //制約条件の数
} SimplexTableau_t;

void print_problem(const SimplexTableau_t *simplex){
    printf( "目的関数\nz = ");
    for(int i = 1; i <= simplex->variable_count; i ++){
        printf("%2fx_%d", simplex->objective_function[i], i);
        printf(" + ");
    }
    printf("%f", simplex->objective_function[RHS_COLUMN]);

    printf("  (最大化)\n");

    printf("\n制約条件\n");
    for(int i = 0; i < simplex->constraint_count; i ++){
        for(int j = 1; j <= simplex->variable_count; j ++){
            printf("%2fx_%d", simplex->constraints[i][j], j);
                if(j != simplex->variable_count){
                    printf(" + ");
            }
        }
        printf(" = %2f\n", simplex->constraints[i][RHS_COLUMN]);
    }

    printf("x_i>= 0  (i = 1,2,...,%d)\n", simplex->variable_count);
}

int input_int(const char *message, int min, int max){
    int value;

    while(1){
        fprintf(stderr, "%s", message);

        if(scanf("%d", &value) == 1){
            if(value >= min && value <= max){
                return value;
            }
        }

        fprintf(stderr, "入力が不正です．%d〜%dの範囲で入力してください．\n", min, max);

        while(getchar() != '\n');
    }
}

double input_double(const char *message)
{
    double value;

    while(1){
        fprintf(stderr, "%s", message);

        if(scanf("%lf", &value) == 1){
            return value;
        }

        fprintf(stderr, "数値を入力してください．\n");

        while(getchar() != '\n');
    }
}

void in_data(SimplexTableau_t *simplex)
{
    simplex->variable_count = input_int("変数の数を入力してください．-> ", 1, MAX_VARIABLES);

    fprintf(stderr, "目的関数の係数の入力\n");

    char message[50];
    for(int i = 1; i <= simplex->variable_count; i++){

        sprintf(message, "x_%dの係数：", i);

        simplex->objective_function[i] = input_double(message);
    }

    simplex->objective_function[RHS_COLUMN] = input_double("定数：");


    simplex->constraint_count = input_int("制約条件の式の数を入力してください．-> ", 1, MAX_CONSTRAINTS);

    fprintf(stderr, "制約条件の入力\n");

    for(int i = 0; i < simplex->constraint_count; i++){

        fprintf(stderr,"制約条件 %d\n", i + 1);

        for(int j = 1; j <= simplex->variable_count; j++){
            char message[50];

            sprintf(message, "x_%dの係数：", j);

            simplex->constraints[i][j] = input_double(message);
        }

        simplex->constraints[i][RHS_COLUMN] = input_double("右辺の定数：");
    }

    print_problem(simplex);
}

int get_last_column(const SimplexTableau_t *simplex, int phase){
    return simplex->variable_count + (phase == PHASE1 ? simplex->constraint_count : 0);
}

void print_tableau(const SimplexTableau_t *simplex, int phase){  //phaseタブロの表示
    int last_column = get_last_column(simplex, phase);

    fprintf(stderr, "\n            定数    ");
    for(int i = 1; i <= last_column; i ++){
        fprintf(stderr, "  -x_%d    ", i);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, " z   = ");
    for(int i = 0; i <= last_column; i ++){
        fprintf(stderr," %8.2f ", simplex->tableau[Z_ROW][i]);
    }
    fprintf(stderr, "\n");

    if(phase==PHASE1){
        fprintf(stderr, " -w  = ");
        for(int i = 0; i <= last_column; i ++){
            fprintf(stderr," %8.2f ", simplex->tableau[W_ROW][i]);
        }
        fprintf(stderr, "\n");
    }

    for(int i = 0; i < simplex->constraint_count; i ++){
        fprintf(stderr, " x_%d = ", simplex->basis[i]);
        for(int j = 0; j <= last_column; j ++){
            fprintf(stderr," %8.2f ", simplex->tableau[CONSTRAINT_ROW + i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void create_initial_tableau(SimplexTableau_t *simplex){
    //配列basisの決定
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
        simplex->tableau[Z_ROW][i] = -simplex->objective_function[i];
    }
    //人為変数x_variable_count+1からx_variable_count+constraint_countまで
    for(int i = 0; i < simplex->constraint_count; i ++){
        for(int j = 0; j <= simplex->variable_count; j ++){
            simplex->tableau[CONSTRAINT_ROW + i][j] = simplex->constraints[i][j];
        }
    }
    //人為変数導入後の目的関数-w
    for(int i = 0; i <= simplex->variable_count; i ++){
        simplex->tableau[W_ROW][i] = 0;
        for(int j = 0; j < simplex->constraint_count; j ++){
            simplex->tableau[W_ROW][i] += simplex->tableau[CONSTRAINT_ROW + j][i];
        }
        simplex->tableau[W_ROW][i] *= -1;
    }

    print_tableau(simplex, PHASE1);
}

int find_leaving_line(const SimplexTableau_t *simplex, int leaving){
    for (int i = 0; i < simplex->constraint_count; i++) {
        if (simplex->basis[i] == leaving) {
            return CONSTRAINT_ROW + i;
        }
    }

    fprintf(stderr, "内部エラー: leaving_line が見つかりません\n");
    exit(EXIT_FAILURE);
}

void update_constraint_rows(SimplexTableau_t *simplex, int entering, int leaving, double previous_tableau[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS], int leaving_line, int last_column){
    simplex->basis[leaving_line - CONSTRAINT_ROW] = entering;

    for(int i = 0; i <= last_column; i ++){
        if(i == leaving){
            continue;
        }
        simplex->tableau[leaving_line][i] /= previous_tableau[leaving_line][entering];
    }
    simplex->tableau[leaving_line][leaving] = 1 / previous_tableau[leaving_line][entering];
    simplex->tableau[leaving_line][entering] = NAN;

    //残りの基底変数の行の操作
    for(int i = 0; i < simplex->constraint_count; i ++){
        if(simplex->basis[i] != entering){
            for(int j = 0; j <= last_column; j ++){
                if(j == leaving){
                    continue;
                }
                simplex->tableau[CONSTRAINT_ROW + i][j] -= simplex->tableau[leaving_line][j] * previous_tableau[CONSTRAINT_ROW + i][entering];
            }
            simplex->tableau[CONSTRAINT_ROW + i][leaving] = - simplex->tableau[leaving_line][leaving] * previous_tableau[CONSTRAINT_ROW + i][entering];
        }
    }
}

void update_objective_rows(SimplexTableau_t *simplex, int entering, int leaving, double previous_tableau[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS], int leaving_line, int last_column){

    for(int j = 0; j <= last_column; j ++){
        simplex->tableau[Z_ROW][j] -= simplex->tableau[leaving_line][j] * previous_tableau[Z_ROW][entering];
        simplex->tableau[W_ROW][j] -= simplex->tableau[leaving_line][j] * previous_tableau[W_ROW][entering];
        if(j == leaving){
            simplex->tableau[Z_ROW][j] = - previous_tableau[Z_ROW][entering] * simplex->tableau[leaving_line][leaving];
            simplex->tableau[W_ROW][j] = - previous_tableau[W_ROW][entering] * simplex->tableau[leaving_line][leaving];
        }
    }
}

void pivot(SimplexTableau_t *simplex, int entering, int leaving, int phase){
    double previous_tableau[MAX_CONSTRAINTS + 2][MAX_TABLEAU_COLUMNS];

    for(int i = 0; i < MAX_CONSTRAINTS + 2; i ++){
        for(int j = 0; j < MAX_TABLEAU_COLUMNS; j ++){
            previous_tableau[i][j] = simplex->tableau[i][j];
        }
    }

    int leaving_line;

    leaving_line = find_leaving_line(simplex, leaving);

    int last_column = get_last_column(simplex, phase);

    update_constraint_rows(simplex, entering, leaving, previous_tableau, leaving_line, last_column);

    update_objective_rows(simplex, entering, leaving, previous_tableau, leaving_line, last_column);

    if(phase == PHASE1){
        print_tableau(simplex, PHASE1);
    }else{
        print_tableau(simplex, PHASE2);
    }
}

int determine_leaving_variable(const SimplexTableau_t *simplex, int entering_variable, int phase){
    int leaving_variable;  //交換する基底変数

    double min_ratio = DBL_MAX;
    int leaving_row = -1;

    for(int i = 0; i < simplex->constraint_count; i ++){
        if(simplex->tableau[CONSTRAINT_ROW + i][entering_variable] > 0 && min_ratio > simplex->tableau[CONSTRAINT_ROW + i][RHS_COLUMN] / simplex->tableau[CONSTRAINT_ROW + i][entering_variable]){  //最小添え字規則
            min_ratio = simplex->tableau[CONSTRAINT_ROW + i][RHS_COLUMN] / simplex->tableau[CONSTRAINT_ROW + i][entering_variable];
            leaving_row = CONSTRAINT_ROW + i;
        }
    }

    if (leaving_row == -1) {
        fprintf(stderr, "内部エラー: leaving_row が見つかりません\n");
        exit(EXIT_FAILURE);
    }

    leaving_variable = simplex->basis[leaving_row  - CONSTRAINT_ROW];

    return leaving_variable;
}

int is_basis(const SimplexTableau_t *simplex, int variable){
    for(int i = 0; i < simplex->constraint_count; i++){
        if(simplex->basis[i] == variable){
            return 1;  // 基底変数
        }
    }
    return 0;  // 非基底変数
}

int determine_entering_variable(const SimplexTableau_t *simplex, int phase){
    int entering_variable = -1;  //交換する非基底変数

    int last_column = get_last_column(simplex, phase);

    int objective_value = (phase == PHASE1) ? W_ROW : Z_ROW;

    // 最初の非基底変数を探す
    for(int i = 1; i <= last_column; i++){
        if(!is_basis(simplex, i)){
            entering_variable = i;
            break;
        }
    }
    if(entering_variable == -1){
        fprintf(stderr, "内部エラー: entering_variable が見つかりません\n");
        exit(EXIT_FAILURE);
    }

    // 目的関数係数が最小の非基底変数を探す
    for(int i = 1; i <= last_column; i++){
        // 基底変数は候補から除外
        if(is_basis(simplex, i)){
            continue;
        }

        // NANは除外
        if(isnan(simplex->tableau[objective_value][i])){
            continue;
        }

        if(simplex->tableau[objective_value][entering_variable]
            > simplex->tableau[objective_value][i]){

            entering_variable = i;
        }
    }
    return entering_variable;
}

void determine_entering_and_leaving_variable(SimplexTableau_t *simplex, int phase){
    int entering_variable;  //交換する非基底変数
    int leaving_variable;  //交換する基底変数

    entering_variable = determine_entering_variable(simplex, phase);

    leaving_variable = determine_leaving_variable(simplex, entering_variable, phase);

    fprintf(stderr, "\nx_%dとx_%dの交換\n", entering_variable, leaving_variable);

    pivot(simplex, entering_variable, leaving_variable, phase);
}

int phase1(SimplexTableau_t *simplex){
    fprintf(stderr, "\nphase1\n");
    double previous_value_of_negative_w = simplex->tableau[W_ROW][RHS_COLUMN];  //前回の-wの目的関数値

    while(fabs(simplex->tableau[W_ROW][RHS_COLUMN]) > EPSILON){
        determine_entering_and_leaving_variable(simplex, PHASE1);
        if(previous_value_of_negative_w > simplex->tableau[W_ROW][RHS_COLUMN]){
            fprintf(stderr, "\n実行不可能です\n");
            return 1;
        }
        previous_value_of_negative_w = simplex->tableau[W_ROW][RHS_COLUMN];  //更新
    }
    fprintf(stderr, "\n実行可能です\n");
    return 0;
}

void print_answer(const SimplexTableau_t *simplex){
    int flag;
    printf("\n最適解は\n");
    for(int i = 1; i <= simplex->variable_count; i ++){
        flag = 0;
        for(int j = 0; j < simplex->constraint_count; j ++){
            if(i == simplex->basis[j]){
                printf("x_%d= %lf", i, simplex->tableau[CONSTRAINT_ROW + j][RHS_COLUMN]);
                flag = 1;
                break;
            }
        }
        if(flag == 0){
            printf("x_%d= %d", i, 0);
        }
        if(i != simplex->variable_count){
            printf(" ");
        }
    }
    printf("\n");

    printf("最適値は %lf", simplex->tableau[Z_ROW][RHS_COLUMN]);
}

int has_negative_coefficient(const SimplexTableau_t *simplex){
    for(int i = 1; i <= simplex->variable_count; i++){
        if(simplex->tableau[Z_ROW][i] < 0){
            return 1;
        }
    }
    return 0;
}

void phase2(SimplexTableau_t *simplex){
    fprintf(stderr, "\nphase2\n");

    print_tableau(simplex, PHASE2);

    while(has_negative_coefficient(simplex)){
        determine_entering_and_leaving_variable(simplex, PHASE2);
    }

    print_answer(simplex);
}

int main(void){
    SimplexTableau_t simplex;

    in_data(&simplex);  //データの入力

    create_initial_tableau(&simplex);  //初期タブロの作成

    if(phase1(&simplex) == 0){
        phase2(&simplex);  //フェーズ2
    }

    return 0;
}