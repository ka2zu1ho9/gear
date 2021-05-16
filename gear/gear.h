#include <stdio.h>
#include <math.h>

/* 計算モードの数 */
#define MODE_NUMBER 2
/* モータ出力(W[rpm])*/
#define OUTPUT_OF_MOTOR 7500

enum{
    known_sokuhi = 0,
    unknown_sokuhi = 1,
    calculate_error_of_sokuhi = 2
};

int gcd(int large_number, int small_number);

void processor_base_tooth_combination(double sokuhi, double error, int z_min, int z_max, int is_pinion);
void processor_other_tooth_combination(double target_sokuhi, double error);

void calculate_tooth_combination();