#include <stdio.h>
#include <math.h>

/* 計算モードの数 */
#define MODE_NUMBER 3
/* モータ出力(W[rpm])*/
#define OUTPUT_OF_MOTOR 7500

// const double elastic_module = 2 * 10*10*10*10*10*10;
// const double poason_rate = 0.33;
// const double shear_module = 0.8 * 10*10*10*10*10*10;

enum{
    known_sokuhi = 0,
    unknown_sokuhi = 1,
    calculate_error_of_sokuhi = 2
};

void input_int(char str[], int* data);
void input_double(char str[], double* data);

int gcd(int large_number, int small_number);

double inv(double angle);

void cal_base_pitch_circle(double module, int z);
void cal_engagement_pitch_circle(int mode, int z1, int z2, double a);
void cal_engagement_press_angle(int z1, int z2, double center_distance_increase);
void cal_backlush_on_pitch_circle(double engagement_press_angle, double module, int z1, int z2);
void cal_tolerance_unit(double base_pitch_circle, double module);
void cal_center_distance_increase(double module, double a, int z1, int z2);
void cal_tooth_tip_circle(int z, double module, double center_distance_increase, double transfer);
void cal_gear_base_circle(double base_pitch_circle);
void cal_engagement_length(double tooth_tip_radius[], double gear_base_radius[], double a, double press_angle);

void cal_sokuhi_from_giving_condition(
    double weight_object, 
    double weight_one_unit, 
    double shita, 
    double u, 
    double radius_of_winding
);

void pro_base_tooth_combination(double sokuhi, double error, int z_min, int z_max, int is_pinion);
void pro_other_tooth_combination(double target_sokuhi, double error);
void pro_tooth_combination();

void pro_gear_strength();