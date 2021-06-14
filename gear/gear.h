#include <stdio.h>
#include <math.h>

/* 計算モードの数 */
#define MODE_NUMBER 3
/* モータ出力(W[rpm])*/
#define OUTPUT_OF_MOTOR 7500

enum{
    known_sokuhi = 0,
    unknown_sokuhi = 1,
    calculate_error_of_sokuhi = 2
};

int gcd(int large_number, int small_number);

void cal_engagement_pitch_circle(double module, int z);

void cal_sokuhi_from_giving_condition(
    double weight_object, 
    double weight_one_unit, 
    double shita, 
    double u, 
    double radius_of_winding
);

void pro_base_tooth_combination(double sokuhi, double error, int z_min, int z_max, int is_pinion);
void pro_other_tooth_combination(double target_sokuhi, double error);

void cal_tooth_combination();