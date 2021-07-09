#include "shaft.h"

#include <stdio.h>
#include <math.h>

const double C[5] = { 0.71, 0.016, 0.1, 0.07, 5.75};

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

void input_int(char str[], int* data);
void input_double(char str[], double* data);

void cal_guzai();
void cal_notch_coefficient();
void cal_total_moment_amplitude();
void cal_using_shear_stress();
void cal_safe_rate();

void pro_shaft_strengh();