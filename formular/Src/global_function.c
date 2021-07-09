#include "../Inc/global_function.h"

void input_int(char str[], int* data){
    printf("%s", str);
    scanf("%d", data);
}

void input_double(char str[], double* data){
    printf("%s", str);
    scanf("%lf", data);
}