#include "../Inc/shaft.h"

double shaft_cofficient[5] = { 0.71, 0.016, 0.1, 0.07, 5.75};

void cal_guzai(){
    for (size_t i = 0; i < 4; i++){
        switch (i)
        {
        case 0:
            out.guzai[i] = shaft_cofficient[0] + shaft_cofficient[1]*in.sigma_B;
            printf("%f\n", out.guzai[i]);
            break;

        case 1:
            out.guzai[i] = 1 - exp(-shaft_cofficient[2]*in.d);
            printf("%f\n", out.guzai[i]);
            break;

        case 2:
            out.guzai[i] = 1 - exp( (-shaft_cofficient[3]*in.d) / in.rowe);
            printf("%f\n", out.guzai[i]);
            break;

        case 3:
            out.guzai[i] = 1 - exp( -shaft_cofficient[4]*(1 - in.d / in.D) );
            printf("%f\n", out.guzai[i]);
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

    printf("%f\n", out.notch_coefficient);
    
}


void cal_total_moment_amplitude(){
    out.total_moment_amplitude = sqrt(in.moment_amplitude[0]*in.moment_amplitude[0] + in.moment_amplitude[1]*in.moment_amplitude[1]);
    printf("%f\n", out.total_moment_amplitude);
}


void cal_using_shear_stress(){
    out.using_shear_stress = (16*sqrt( ( (out.notch_coefficient*in.sigma_y_p*out.total_moment_amplitude) / in.sigma_e_B )*( (out.notch_coefficient*in.sigma_y_p*out.total_moment_amplitude) / in.sigma_e_B ) + in.Tav*in.Tav )) / (M_PI*in.d*in.d*in.d);
    printf("%f[kg/mm^2]\n", out.using_shear_stress);
}

void cal_shaft_safe_rate(){
    out.safe_rate = in.tau_y_p / out.using_shear_stress;
    printf("%f\n", out.safe_rate);
}

void pro_shaft_strengh(){
    input_int("åãŽ¢åãŽ¼åã²å£òñáŽ¸åâÆå¢©åá¾å¢¢åáÊå¡äåà¤ñŽ¼°å¤¯åãŽ¼ì¼?: 1åà?åáÒå¢Áïâ?:2?Ž¼? = ", &in.mode);
    input_double("ëûŽ²åáÄå¦¤åãŽ¼åãŽ¡åãŽ³åã°è¸±ê»?Mr-x[kgfåãŽ»mm] = ", &in.moment_amplitude[0]);
    input_double("ëûŽ²åáÄå¦¤åãŽ¼åãŽ¡åãŽ³åã°è¸±ê»?Mr-x[kgfåãŽ»mm] = ", &in.moment_amplitude[1]);
    input_double("ê»Ž³éý?åáŽ­åáÐå£êåãŽ¢åãŽ¼åãŽ¡åãŽ³åã?[kgfåãŽ»mm] = ", &in.Tav);
    input_double("ê¾ÊçŽ¼ŽµñùºæŽ¼¾é¤»[kgf/mm^2] = ", &in.sigma_y_p);
    input_double("éûÜêŽ»Ž¢íöŽ²åâ¸ëÑðê¼Ž¦[kgf/mm^2] = ", &in.sigma_e_B);
    input_double("åáÖå£óëöŽ­ñùºæŽ¼¾é¤»[kgf/mm^2] = ", &in.tau_y_p);

    if(in.mode == 1){
        input_double("èºŽ¸ì¥ÄéÖ¶êÀ?d[mm] = ", &in.sigma_y_p);
        input_double("éè?åâ´èŽ¬?åáºæŽ¿¤èÊ² = ", &in.K);
    }else{
        input_double("ê¾ÊçŽ¼ŽµñùºçŽ¼Ž·åáÉ¼kgf/mm^2] = ", &in.sigma_B);
        input_double("èºŽ¸ì¥ÄéÖ¶êÀ??Ž¼°çŽ¤Žª?Ž¼±¥[mm] = ", &in.D);
        input_double("èºŽ¸ì¥ÄéÖ¶êÀ??Ž¼°éŽ´Ž°?Ž¼±Å[mm] = ", &in.d);
        input_double("åáÒå¢Áïâ²ç¹êêÀ§Ñ[mm] = ", &in.rowe);
    }

    cal_guzai();
    cal_notch_coefficient();
    cal_total_moment_amplitude();
    cal_using_shear_stress();
    cal_shaft_safe_rate();
}

