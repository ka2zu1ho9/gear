#ifndef _GEAR_H_
#define _GEAR_H_

#include <stdio.h>
#include <math.h>
#include "global_function.h"


/* ?v?Z???[?h??? */
#define MODE_NUMBER 3
/* ???[?^?o??(W[rpm])*/
#define OUTPUT_OF_MOTOR 7500

#define K1 1.0
#define K2 1.2
#define K3 1.5
#define E  (double)2 * 10*10*10*10*10*10

enum{
    known_sokuhi = 0,
    unknown_sokuhi = 1,
    calculate_error_of_sokuhi = 2
};

double l, te;
double y;
double engagement_press_angle;
double W[2];
double backlush;
double backlush_lmit[2];
double base_pitch_circle[2];
double tooth_tip_circle[2];
double gear_base_circle[2];
double shaft_torque[2];
double vartical_force[2];
double tanquent_force[2];
double bending_stress[2];
double curvanture_radius[2];
double compressive_stress[2];
double safe_rate[2];

int gcd(int large_number, int small_number);

double inv(double angle);

void cal_base_pitch_circle(int mode, double module, int z);
void cal_engagement_pitch_circle(int mode, int z1, int z2, double a);
void cal_engagement_press_angle(int z1, int z2, double center_distance_increase);
void cal_backlush_on_pitch_circle(double module, int z1, int z2);
void cal_tolerance_unit(double base_pitch_circle, double module);
void cal_max_and_min_backlush(int mode);
void cal_center_distance_increase(double module, double a, int z1, int z2);
void cal_tooth_tip_circle(int mode, int z, double module, double center_distance_increase, double transfer);
void cal_gear_base_circle(int mode);
void cal_engagement_length(double a);
void cal_normal_pitch(double module);
void cal_engagement_rate();
void cal_shaft_torque(int mode, double n_in, double H);
void cal_varitical_force(int mode);
void cal_tanquent_force(int mode);
void cal_bending_stress(int mode, double module, double b, double gear_fatigue);
void cal_safe_rate(int cal_mode, int mode,  double fatigue);
void cal_curvanture_radius(int mode);
void cal_compressive_stress(int mode, double b);

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

#endif