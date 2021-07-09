#ifndef _SHAFT_H_
#define _SHAFT_H_

#include <stdio.h>
#include <math.h>
#include "global_function.h"

typedef struct shaft_input{
    int mode;
    double moment_amplitude[2];
    double Tav;
    double sigma_y_p;
    double sigma_e_B;
    double tau_y_p;
    double d;

    // 
    double K;

    //  
    double sigma_B;
    double D;
    double rowe;


}input;

typedef struct shaft_output{
    double guzai[4];
    double notch_coefficient;
    double total_moment_amplitude;
    double using_shear_stress;
    double safe_rate;
}output;

input in;
output out;

extern double shaft_cofficient[5];

void cal_guzai();
void cal_notch_coefficient();
void cal_total_moment_amplitude();
void cal_using_shear_stress();
void cal_shaft_safe_rate();

void pro_shaft_strengh();

#endif