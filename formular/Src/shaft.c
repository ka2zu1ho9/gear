#include "../Inc/shaft.h"

double shaft_cofficient[5] = { 0.71, 0.016, 0.1, 0.07, 5.75};

void cal_guzai(){
    for (size_t i = 0; i < 4; i++){
        switch (i)
        {
        case 0:
            out.guzai[0] = shaft_cofficient[0] + shaft_cofficient[1]*in.sigma_B;
            printf("Ì1 = %f\n", out.guzai[0]);
            break;

        case 1:
            out.guzai[1] = 1 - exp(-shaft_cofficient[2]*in.d);
            printf("Ì2 = %f\n", out.guzai[1]);
            break;

        case 2:
            out.guzai[2] = 1 - exp( (-shaft_cofficient[3]*in.d) / in.rowe);
            printf("Ì3 = %f\n", out.guzai[2]);
            break;

        case 3:
            out.guzai[3] = 1 - exp( -shaft_cofficient[4]*(1 - in.d / in.D) );
            printf("Ì4 = %f\n", out.guzai[3]);
            break;
        default:
            break;
        }
    }
}

void cal_notch_coefficient(){
    double total_guzai;

    for (size_t i = 0; i < 4; i++)
    {
        total_guzai *= out.guzai[i];
    }

    out.notch_coefficient = 1 + total_guzai;

    printf("·Ý÷ÌØè«WK = %f\n", out.notch_coefficient);
    
}


void cal_total_moment_amplitude(){
    out.total_moment_amplitude = sqrt(in.moment_amplitude[0]*in.moment_amplitude[0] + in.moment_amplitude[1]*in.moment_amplitude[1]);
    printf("È°[gÌUMr = %f[kgfEmm]\n", out.total_moment_amplitude);
}


void cal_using_shear_stress(){
    out.using_shear_stress = (16*sqrt( ( (out.notch_coefficient*in.sigma_y_p*out.total_moment_amplitude) / in.sigma_e_B )*( (out.notch_coefficient*in.sigma_y_p*out.total_moment_amplitude) / in.sigma_e_B ) + in.Tav*in.Tav )) / (M_PI*in.d*in.d*in.d);
    printf("²Ìgp¹ñfÍÑ = %f[kg/mm^2]\n", out.using_shear_stress);
}

void cal_shaft_safe_rate(){
    out.safe_rate = in.tau_y_p / out.using_shear_stress;
    printf("%f\n", out.safe_rate);
}

void pro_shaft_strengh(){
    input_int("L[a:1,·Ý÷:2 = ", &in.mode);
    input_double("È°[gÌUx¬ª Mr_x[kgfEmm] = ", &in.moment_amplitude[0]);
    input_double("È°[gÌUy¬ª Mr_y[kgfEmm] = ", &in.moment_amplitude[1]);
    input_double("½ÏË¶è[g Tav[kgfEmm] = ", &in.Tav);
    input_double("ø£~_ Ðy_p[kgf/mm^2] = ", &in.sigma_y_p);
    input_double("ñ]È°æêÀx Ðe_B[kgf/mm^2] = ", &in.sigma_e_B);
    input_double("¹ñf~_ Ñy_p[kgf/mm^2] = ", &in.tau_y_p);

    if(in.mode == 1){
        input_double("¼a d[mm] = ", &in.d);
        input_double("Øè«W K = ", &in.K);
    }else{
        input_double("ø£­³ ÐB[kgf/mm^2] = ", &in.sigma_B);
        input_double("¼ai¾jD[mm] = ", &in.D);
        input_double("¼ai×jd[mm] = ", &in.d);
        input_double("·Ý÷¼a Ï[mm] = ", &in.rowe);
    }

    if(in.mode == 1){
        cal_total_moment_amplitude();
        cal_using_shear_stress();
        cal_shaft_safe_rate();
    }else{
        cal_guzai();
        cal_notch_coefficient();
        cal_total_moment_amplitude();
        cal_using_shear_stress();
        cal_shaft_safe_rate();
    }
}

