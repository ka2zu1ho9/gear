#include <stdio.h>
#include "../Inc/gear.h"
#include "../Inc/shaft.h"

int main(void){
    int is_keeping_while = 1;
    int mode = 0;
    int cnt = 0;

    while(1){
        if(cnt > 0){
            printf("計算を続けますか？\n");
            printf("yes:0,no:1\n");
            scanf("%d", &is_keeping_while);   
            printf("\n");

            if(is_keeping_while == 1){
                break;
            }
        }

        input_int("mode = (0:tooth_combination, 1:gear_strength, 2:shaft_strength)", &mode);

        if(mode == 0){
            pro_tooth_combination();
        }else if(mode == 1){
            pro_gear_strength();
        }else if(mode == 2){
            pro_shaft_strengh();
        }
        
        cnt = 1;
    }

    return 0;
}