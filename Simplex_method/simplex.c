#include <stdio.h>

#define MAX_VARIABLES 10  //変数の最大値
#define MAX_CONSTRAINTS 20  //制約条件の式の数の最大値

int objective_function[MAX_VARIABLES];  //目的関数
int constraints[MAX_CONSTRAINTS][MAX_VARIABLES];  //制約条件

int var;  //変数の数
int con;  //制約条件の数

void check_data(){
    fprintf(stderr, "目的関数  z = ");
    for(int i = 0; i < var; i ++){
        fprintf(stderr, "%2dx_%d", objective_function[i], i);
        fprintf(stderr, " + ");
    }
    fprintf(stderr, "  (最大化)\n");

    fprintf(stderr, "制約条件\n");
    for(int i = 0; i < con; i ++){
        for(int j = 0; j < var; j ++){
            fprintf(stderr, "%2dx_%d", constraints[i][j], j);
            fprintf(stderr, " + ");
        }
        fprintf(stderr, " = %2d\n", constraints[i][var + 1]);
    }
}

void in_data(){
    fprintf(stderr, "変数の数を入力してください．-> ");
    scanf("%d", &var);

    fprintf(stderr, "目的関数の係数の入力\n");
    for(int i = 0; i < var; i ++){
        fprintf(stderr, "x_%dの係数：", i);
        scanf("%d", &objective_function[i]);
    }

    fprintf(stderr, "制約条件の式の数の数を入力してください．-> ");
    scanf("%d", &con);

    fprintf(stderr, "制約条件の入力\n");
    for(int i = 0; i < con; i ++){
            fprintf(stderr, "制約条件 %d\n", i);
        for(int j = 0; j < var; j ++){
            fprintf(stderr, "x_%dの係数：", j);
            scanf("%d", &constraints[i][j]);
        }
        fprintf(stderr, "右辺の定数：");
        scanf("%d", &constraints[i][var + 1]);
    }
    check_data();
}

int main(void){

    in_data();

    return 0;
}