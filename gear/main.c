#include <stdio.h>
#include "gear.h"

int main(void){
    int is_keeping_while = 1;
    int cnt = 0;

    while(1){
        if(cnt > 0){
            printf("Do you want to continue the calculation?\n");
            printf("yes:0,no:1\n");
            scanf("%d", &is_keeping_while);   
            printf("\n");

            if(is_keeping_while == 1){
                break;
            }
        }
        calculate_tooth_combination();
        cnt ++;
    }

    return 0;
}