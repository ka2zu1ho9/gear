#include "../Inc/gear.h"

int gcd(int large_number, int small_number){
    int temp;
    int dividend,divisor;
    int r = 1;
    int cnt = 0;

    if(large_number < small_number){
        temp = large_number;
        large_number = small_number;
        small_number = temp;
    }

    dividend = large_number;
    divisor  = small_number;

    while(r != 0){
        if(cnt > 0){
            dividend = divisor;
            divisor  = r;
        }

        r = dividend % divisor;
        cnt ++;
    }

    return divisor;
}

double inv(double angle){
    return (tan(angle * (M_PI/180)) - angle*(M_PI/180));
}

void cal_base_pitch_circle(int mode, double module, int z){
    if(mode == 1){
        base_pitch_circle[0] = module*z;
        printf("%.4f[mm]\n", module*z);
    }else{
        base_pitch_circle[1] = module*z;
        printf("%.4f[mm]\n", module*z );
    }
}

void cal_engagement_pitch_circle(int mode, int z1, int z2, double a){
    if(mode == 1){
        printf("%.4f[mm]\n", (z1 * a) / (z1 + z2) );
    }else if(mode == 2){
        printf("%.4f[mm]\n", (z2 * a) / (z1 + z2) );
    }
}

void cal_backlush_on_pitch_circle(double module, int z1, int z2){

    backlush = ( cos(20*(M_PI/180.0)) / cos(engagement_press_angle*(M_PI/180.0)) )*module*(z1 + z2)*(inv(engagement_press_angle) - inv(20));
    printf("%.4f[mm]\n", backlush);
}

void cal_tolerance_unit(double mode, double module){
    if(mode == 1){
        W[0] = cbrt(base_pitch_circle[0]) + 0.65*module;
        printf("%.4f[??m]\n", W[0]);
    }else if(mode == 2){
        W[1] = cbrt(base_pitch_circle[1]) + 0.65*module;
        printf("%.4f[??m]\n", W[1]);
    }
}

void cal_max_and_min_backlush(int mode){
    double gear1_backlush[2];
    double gear2_backlush[2];
    if(mode == 1){
        //???
        gear1_backlush[0] = 35.5*(W[0] / pow(10,3));
    }else{
        //???
        gear1_backlush[1] = 10.0*(W[0] / pow(10,3));
    }
    
    if(mode == 1){
        // ???
        gear2_backlush[0] = 35.5*(W[1] / pow(10,3));
    }else{
        // ???
        gear2_backlush[1] = 10.0*(W[1] / pow(10,3));
    }

    if(mode == 1){
        backlush_lmit[0] = gear1_backlush[0] + gear2_backlush[0];
        printf("%f\n", backlush_lmit[0]);
    }else{
        backlush_lmit[1] = gear1_backlush[1] + gear2_backlush[1];
        printf("%f\n", backlush_lmit[1]);
    }


}

void cal_center_distance_increase(double module, double a, int z1, int z2){
    double result;

    y = (a/module) - ( (double)(z1 + z2)/2.0);

    printf("%.4f\n", y);
}

void cal_engagement_press_angle(int z1, int z2, double center_distance_increase){
    double result;

    result = acos( ( (double)(z1 + z2) / (double)(z1 + z2 + 2*y) ) * cos(20*M_PI/180.0)) * (180.0/M_PI);
    engagement_press_angle = result;
    printf("%.4f[deg]\n", result);
}

void cal_tooth_tip_circle(int mode, int z, double module, double center_distance_increase, double transfer){
    if(mode == 1){
        tooth_tip_circle[0] = (z + 2) * module + 2*module*(center_distance_increase - transfer);
        printf("%.4f[mm]\n", tooth_tip_circle[0]);
    }else{
        tooth_tip_circle[1] = (z + 2) * module + 2*module*(center_distance_increase - transfer);
        printf("%.4f[mm]\n", tooth_tip_circle[1]);
    }
}

void cal_gear_base_circle(int mode){
    if(mode == 1){
        gear_base_circle[0] = (base_pitch_circle[0]/2)*cos(20*(M_PI/180));
        printf("%.4f[mm]\n", gear_base_circle[0]);
    }else{
        gear_base_circle[1] = (base_pitch_circle[1]/2)*cos(20*(M_PI/180));
        printf("%.4f[mm]\n", gear_base_circle[1]);
    }
}

void cal_engagement_length(double a){
    double result;

    l = sqrt( (tooth_tip_circle[0]/2*tooth_tip_circle[0]/2) - gear_base_circle[0]*gear_base_circle[0]) + sqrt( (tooth_tip_circle[1]/2*tooth_tip_circle[1]/2) - gear_base_circle[1]*gear_base_circle[1] ) - a*sin(engagement_press_angle*(M_PI/180));

    printf("%.4f[mm]\n", l);
}

void cal_normal_pitch(double module){
    te = M_PI*module*cos(20*M_PI/180.0);
    printf("%.4f[mm]\n", te);
}

void cal_engagement_rate(){
    printf("%.4f\n", l/te);
}

void cal_shaft_torque(int mode, double n_in, double H){
    if(mode == 1){
        shaft_torque[0] = ((75 * 60 * H) / (2*M_PI*n_in));
        printf("%.4f[kgf?Emm]\n", shaft_torque[0]);
    }else{
        shaft_torque[1] = ((75 * 60 * H) / (2*M_PI*n_in));
        printf("%.4f[kgf?Emm]\n", shaft_torque[1]);
    }
}

void cal_varitical_force(int mode){
    if(mode == 1){
       vartical_force[0] = shaft_torque[0]*10*10*10/ (gear_base_circle[0]);
       printf("%.4f[kgf]\n", vartical_force[0]);
    }else{
       vartical_force[1] = shaft_torque[1]*10*10*10 / gear_base_circle[1];
       printf("%.4f[kgf]\n", vartical_force[1]);
    }
}

void cal_tanquent_force(int mode){
    if(mode == 1){
        tanquent_force[0] = vartical_force[0] * cos(20*M_PI/180.0);
        printf("%.4f[kgf]\n", tanquent_force[0]);
    }else{
        tanquent_force[1] = vartical_force[1] * cos(20*M_PI/180.0);
        printf("%.4f[kgf]\n", tanquent_force[1]);
    }
}

void cal_bending_stress(int mode, double module, double b, double gear_fatigure){
    if(mode == 1){
        bending_stress[0] = (tanquent_force[0]*gear_fatigure*K1*K2*K3/( module*b*cos(engagement_press_angle*(M_PI/180.0)) ) );
        printf("%.4f[kgf/mm2]\n", bending_stress[0]);
    }else{
        bending_stress[1] = (tanquent_force[1]*gear_fatigure*K1*K2*K3/( module*b*cos(engagement_press_angle*(M_PI/180.0)) ) );
        printf("%.4f[kgf/mm2]\n", bending_stress[1]);
    }
}

void cal_curvanture_radius(int mode){
    if(mode == 1){
        curvanture_radius[0] = gear_base_circle[0]*sin(engagement_press_angle*M_PI/180.0);
        printf("%.4f[mm]\n", curvanture_radius[0]);
    }else{
        curvanture_radius[1] = gear_base_circle[1]*sin(engagement_press_angle*M_PI/180.0);
        printf("%.4f[mm]\n", curvanture_radius[1]);
    }
}

void cal_compressive_stress(int mode, double b){
    double young_ratio = (2*pow(10,9)) / (9.8*pow(10, 6));
    if(mode == 1){
        compressive_stress[0] = sqrt( (0.175*vartical_force[0]*young_ratio*K1*K2*K3)/(b*curvanture_radius[0]) );
        printf("%.4f[kgf/mm2]\n", compressive_stress[0]);
    }else{
        compressive_stress[1] = sqrt( (0.175*vartical_force[1]*young_ratio*K1*K2*K3)/(b*curvanture_radius[1]) );
        printf("%.4f[kgf/mm2]\n", compressive_stress[1]);
    }
}

void cal_safe_rate(int cal_mode, int mode,  double fatigue){
    if(cal_mode == 1){
        // ???
        if(mode == 1){
            safe_rate[0] = ( ( (45*pow(10, 6)) / (9.8*pow(10, 6)) ) / bending_stress[0] );
            printf("%.4f\n",  safe_rate[0]);
        }else{
            safe_rate[0] = ( ( (45*pow(10, 6)) / (9.8*pow(10, 6)) ) / bending_stress[1]);
            printf("%.4f\n", safe_rate[0]);
        }
    }else{
        // ???k
        if(mode == 1){
            safe_rate[1] = ( ( (24.5*pow(10, 6)) / (9.8*pow(10, 6)) ) / compressive_stress[0]);
            printf("%.4f\n", safe_rate[1]);
        }else{
            safe_rate[1] = ( ( (24.5*pow(10, 6)) / (9.8*pow(10, 6)) ) / compressive_stress[1]);
            printf("%.4f\n", safe_rate[1]);
        }
    }
}

void pro_base_tooth_combination(
    double sokuhi, 
    double error, 
    int z_min, 
    int z_max, 
    int is_pinion)
{
    double sokuhi1, gosa, rate, repeat_frequency, gcd_result;

    if(is_pinion == 1){
        for(size_t i = z_min; i <= z_max - 1; i++){

            for (size_t j = z_min + 1; j <= z_max; j++){
                if(i < j){
                    sokuhi1 = ( (double)j / (double)i );
                    rate = sokuhi1 / sokuhi;
                    gosa = sokuhi1 - sokuhi;

                    repeat_frequency = (j*i)/(gcd(j,i)*gcd(j,i));
                    gcd_result = gcd(j,i);

                    if( fabs(gosa) < sokuhi * (error/100) ){
                        printf("%d\t%3d\t%.4f\t%.4f\t%.4f(%.4f)\n", i, j, sokuhi1, gosa, repeat_frequency, gcd_result);
                    }
                } 
            }

        }
    }else if(is_pinion == 0){
        for (size_t j = z_min + 1; j <= z_max; j++){
            sokuhi1 = ( (double)j / (double)z_min );
            rate = sokuhi1 / sokuhi;
            gosa = sokuhi1 - sokuhi;

            repeat_frequency = (j*z_min)/(gcd(j,z_min)*gcd(j,z_min));
            gcd_result = gcd(j,z_min);

            if( fabs(gosa) < sokuhi * (error/100) ){
                printf("%d\t%3d\t%.4f\t%.4f\t%.4f(%.4f)\n", z_min, j, sokuhi1, gosa, repeat_frequency, gcd_result);
            }
        }   
    }
}

void pro_other_tooth_combination(double target_sokuhi, double target_error){    
    int increase_first_stage_gear_ratio_flag = 0;
    double loop_v,r;
    double v[3];
    double real_sokuhi;
    double real_error;

    v[0] = pow(target_sokuhi, 1.0/3.0);

    printf("-----------------------------------------------------------------------------------------------\n");
    printf("1?i??????\t2?i??????\t3?i??????\t?S??????\n");
    printf("-----------------------------------------------------------------------------------------------\n");

    if(increase_first_stage_gear_ratio_flag == 0){
        for(double k = -3.0; k < 3.0; k += 0.0005){
            v[1] = v[0] + k;
            v[2] = v[1] + 2*k;
            real_sokuhi = v[0]*v[1]*v[2];

            if( ( (1-(target_error/100.0))*target_sokuhi < real_sokuhi ) && ( real_sokuhi < (1+(target_error/100.0))*target_sokuhi ) ){
                real_error = ( (real_sokuhi - target_sokuhi) / target_sokuhi )*100;
                printf("%.4f\t%.4f\t%.4f\t%.4f\n", v[0], v[1], v[2], real_sokuhi);
            }
        }
    }else if(increase_first_stage_gear_ratio_flag == 1){
        loop_v = pow(target_sokuhi, 1.0/3.0);
        r = loop_v / 3;
        v[0] = loop_v + (2 * r);

        for(double k = -3.0; k < 3.0; k += 0.0005){
            v[1] = loop_v + 2*k;
            v[2] = r + k;
            real_sokuhi = v[0]*v[1]*v[2];

            if( ( (1-(target_error/100.0))*target_sokuhi < real_sokuhi ) && ( real_sokuhi < (1+(target_error/100.0))*target_sokuhi ) 
            && v[0] > 0 && v[1] > 0 && v[2] > 0){

                real_error = ( (real_sokuhi - target_sokuhi) / target_sokuhi )*100;
                printf("%.4f\t%.4f\t%.4f\t%.4f\n", v[0], v[1], v[2], real_sokuhi);
            }
        }
    }
}

void pro_gear_strenrth(){
    int cnt = 0;
    int z[2] = {0};
    int module;

    double a,b,gear_coefficient;

    double n_in;
    double H;

    do
    {
        if(cnt > 0){
            printf("???????v???????B\n");
        }
        printf("????1 = \n");
        scanf("%d", &z[0]);

    } while (z[0] < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("???????v???????B\n");
        }
        printf("????2 = \n");
        scanf("%d", &z[1]);

    } while (z[1] < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("???????v???????B\n");
        }
        printf("???W???[?? = \n");
        scanf("%d", &module);

    } while (module < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("???????v???????B\n");
        }
        printf("???S????[mm] = \n");
        scanf("%d", &a);

    } while (a < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("???????v???????B\n");
        }
        printf("????[mm] = \n");
        scanf("%d", &a);

    } while (b < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("???????v???????B\n");
        }
        printf("?????W?? = \n");
        scanf("%d", &gear_coefficient);

    } while (gear_coefficient < 0);
    printf("\n");
    cnt = 0;            

    do
    {
        if(cnt > 0){
            printf("???????v???????B\n");
        }
        printf("?????W?? = \n");
        scanf("%d", &n_in);

    } while (n_in < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("???????v???????B\n");
        }
        printf("????`?B?W?? = \n");
        scanf("%d", &H);

    } while (H < 0);
    printf("\n");
    cnt = 0;
}

void pro_tooth_combination(){
    double sokuhi = 0;
    double error = 0;

    double v[3];
    double target_sokuhi = 0;
    double target_error = 0;

    int z[6] = {0};
    double module[3] = {0};

    int z_min = 0;
    int z_max = 0;
    int is_pinion = 2;
    int do_cal_engagement_pitch_circle = 2;
    int cnt = 0;
    int mode = 0;

    /* ???[?h????? */
    do
    {
        if(cnt > 0){
            printf("???????v???????B\n");
        }
        printf("???[?h??I??????????????B(0:1?i?????g????????????, 1:?S????????A?e?i?????g????????????, 2:?????????????????????)\n");
        scanf("%d", &mode);

    } while (mode < 0);
    printf("\n");
    cnt = 0;
    
    switch (mode)
    {
    case known_sokuhi:
        /* ???????? */
        do{
            printf("????????\n");
            printf("???? > 1\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("???? = ");
            scanf("%lf", &sokuhi);

            cnt ++;
        } while( sokuhi <= 1);
        printf("\n");
        cnt = 0;

        /* ???z?????????e??????? */
        do{
            printf("???z?????????e???????\n");
            printf("???z?????????e??? > 0\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("???z?????????e??? = ");
            scanf("%lf", &error);

            cnt ++;
        } while( error <= 0 );
        printf("\n");
        cnt = 0;

        do{
            printf("?s?j?I???????????? yes : 0 , no : 1\n");
            printf("???? = ");
            scanf("%d", &is_pinion);
        }while(is_pinion > 1);
        printf("\n");

        /* ???????????? */
        do{
            printf("????????????\n");
            printf("Minimum gear teeth >= 6\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("??????? = ");
            scanf("%d", &z_min);

            cnt ++;
        } while( z_min < 6 );
        printf("\n");
        cnt = 0;

        /* ?????????? */
        do{
            printf("??????????\n");
            printf("????? > ??????? + 1\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }
            printf("????? = ");
            scanf("%d", &z_max);

            cnt ++;
        } while( z_max < z_min + 1 );
        printf("\n");
        cnt = 0; 

        printf("----------------------------------------------------------------------------------\n");
        printf("z1\tz2\tz2/z1\t\t??\t\t?J?????p?x(??????)\n");
        printf("----------------------------------------------------------------------------------\n");
        pro_base_tooth_combination(sokuhi, error, z_min, z_max, is_pinion);   
        break;

    case unknown_sokuhi:
        /* ?S?????z???????? */
        do
        {
            printf("?S?????z????\n");
            printf("?S?????z???? > 1\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("?????S?????z???? = ");
            scanf("%lf", &target_sokuhi);

            cnt ++;
        } while (target_sokuhi <= 1);
        printf("\n");
        cnt = 0;

        do
        {
            printf("?S?????z?????????e???\n");
            printf("?S?????z?????????e??? > 0\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("?S?????z?????????e???= ");
            scanf("%lf", &target_error);

            cnt ++;            
        } while (target_error <= 0);
        printf("\n");        
        cnt = 0;

        pro_other_tooth_combination(target_sokuhi, target_error);
        printf("\n");

        for(int i = 0; i < 3; i++){
            do
            {
                printf("?g?p????%d?i??????????\n", i+1);
                printf("?g?p????%d?i?????? > 1\n", i+1);

                if(cnt > 0){
                    printf("????????v???????B\n");
                }

                printf("?g?p????%d?i?????? = ", i+1);
                scanf("%lf", &v[i]);

                cnt ++;
            } while ( v[i] <= 1);
            printf("\n");        
            cnt = 0;
        }

        /* ???z?????????e??????? */
        do{
            printf("???z?????????e???????\n");
            printf("???z?????????e??? > 0\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("???z?????????e??? = ");
            scanf("%lf", &error);

            cnt ++;
        } while( error <= 0 );
        printf("\n");
        cnt = 0;

        do{
            printf("?s?j?I???????????? yes : 0 , no : 1\n");
            printf("???? = ");
            scanf("%d", &is_pinion);
        }while(is_pinion > 1);
        printf("\n");

        /* ???????????? */
        do{
            printf("????????????\n");
            printf("Minimum gear teeth >= 6\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("??????? = ");
            scanf("%d", &z_min);

            cnt ++;
        } while( z_min < 6 );
        printf("\n");
        cnt = 0;

        /* ?????????? */
        do{
            printf("??????????\n");
            printf("????? > ??????? + 1\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("????? = ");
            scanf("%d", &z_max);

            cnt ++;
        } while( z_max < z_min + 1 );
        printf("\n");
        cnt = 0; 

        for(int i = 0; i < 3; i++){
            if(i > 0){
                is_pinion = 1;
            }
            printf("----------------------------------------------------------------------------------\n");
            printf("%d?i??????\n", i+1);
            printf("z1\tz2\tz2/z1\t\t??\t\t?J?????p?x(??????)\n");
            printf("----------------------------------------------------------------------------------\n");
            pro_base_tooth_combination(v[i], error, z_min, z_max, is_pinion);   
            printf("\n");
        }

        do{
            if(cnt > 0){
                printf("????????v???????B\n");
            }
            
            printf("????????s?b?`?~???v?Z????????Hy : 0, n : 1\n"); 
            scanf("%d", &do_cal_engagement_pitch_circle);
            cnt ++;
        } while(do_cal_engagement_pitch_circle > 1);
        printf("\n");
        cnt = 0; 

        if(do_cal_engagement_pitch_circle == 0){
            for(int i = 0; i < 3; i++){
                printf("?I??????module????i?????????????????????B\n");
                printf("module[%d] = ", i);
                scanf("%lf", &module[i]);
            }

            for(int i = 0; i < 6; i++){
                printf("?I??????????????i?????????????????????B\n");
                printf("z[%d] = ", i);
                scanf("%d", &z[i]);
            }

            for(int i = 0; i < 6; i++){
                printf("rk%d = ", i);
                if(i < 3){
                   cal_base_pitch_circle(1, module[i], z[i]);
                }else{
                   cal_base_pitch_circle(1, module[i-3], z[i]);
                }
            }
        }
    
        break;

    case calculate_error_of_sokuhi:
        do
        {
            printf("???z????????\n");
            printf("???z???? > 1\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("???z???? = ");
            scanf("%lf", &target_sokuhi);

            cnt ++;            
        } while (target_error > 1);
        printf("\n");        
        cnt = 0;
        
        do
        {
            printf("???z?????????e???????\n");
            printf("???z?????????e??? > 0\n");

            if(cnt > 0){
                printf("????????v???????B\n");
            }

            printf("???z?????????e??? = ");
            scanf("%lf", &target_error);

            cnt ++;            
        } while (target_error <= 0);
        printf("\n");        
        cnt = 0;

        printf("----------------------------------------------------------------------------------\n");
        printf("----------------------------------------------------------------------------------\n");
        printf("???????\t?????\n");
        printf("%.4f\t%.4f\n", (1-target_error)*target_sokuhi, (1+target_error)*target_sokuhi);
        printf("\n");

        break;
    default:
        break;
    }
}

void pro_gear_strength(){
    int z1, z2;
    double module, a, b, gear_fatigue, n_in, H;
    double bending_fatigue;
    double compression_fatigue;

    input_int("z1 = ", &z1);
    input_int("z2 = ", &z2);
    input_double("module = ", &module);
    input_double("a = ", &a);
    input_double("b = ", &b);
    input_double("???^?W?? = ", &gear_fatigue);
    input_double("n_in = ", &n_in);
    input_double("H = ", &H);
    input_double("??????????S???i?????j = ", &bending_fatigue);
    input_double("???k???????S???i?????j = ", &compression_fatigue);

    printf("\n");

    printf("\t\t   ????u = %.4f\n", (double)z2 / (double)z1);
    printf("     ???s?b?`?~???ad01 = ");
    cal_base_pitch_circle(1, module, z1);
    printf("     ???s?b?`?~???ad02 = ");
    cal_base_pitch_circle(2, module, z2);
    printf("   ????????s?b?`???arb1 = ");
    cal_engagement_pitch_circle(1, z1, z2, a);
    printf("   ????????s?b?`???arb2 = ");
    cal_engagement_pitch_circle(2, z1, z2, a);    
    printf("\t      ?????P??W1 = ");
    cal_tolerance_unit(1, module);
    printf("\t      ?????P??W2 = ");
    cal_tolerance_unit(2, module);    
    printf("       ???S?????????W??y = ");
    cal_center_distance_increase(module, a, z1, z2);
    printf("\t    ????????pab = ");
    cal_engagement_press_angle(z1, z2, a);
    printf("?s?b?`?~??o?b?N???b?VC0 = ");
    cal_backlush_on_pitch_circle(module, z1, z2);
    printf("  ???o?b?N???b?VC0_max = ");
    cal_max_and_min_backlush(1);
    printf("  ????o?b?N???b?VC0_min = ");
    cal_max_and_min_backlush(2);
    printf("\t   ?n??~???adk1 = ");
    cal_tooth_tip_circle(1, z1, module, y, 0);
    printf("\t   ?n??~???adk2 = ");
    cal_tooth_tip_circle(2, z2,module, y, 0);
    printf("\t   ??b?~???arg1 = ");
    cal_gear_base_circle(1);
    printf("\t   ??b?~???arg2 = ");
    cal_gear_base_circle(2);    
    printf("\t    ??????????? = ");
    cal_engagement_length(a);
    printf("\t    ?@???s?b?`te = ");
    cal_normal_pitch(module);
    printf("\t     ??????????? = ");
    cal_engagement_rate();
    printf("\t    ????g???NTa = ");
    cal_shaft_torque(1, n_in, H);
    printf("\t    ??????p??Fn = ");
    cal_varitical_force(1);
    printf("\t        ?????P0 = ");
    cal_tanquent_force(1);
    printf("       ????1?????????1 = ");
    cal_bending_stress(1, module, b, gear_fatigue);
    printf("\t ????1????S??Sb = ");
    cal_safe_rate(1, 1, bending_fatigue);
    printf("       ????1???????ap1 = ");
    cal_curvanture_radius(1);
    printf("       ????2???????ap1 = ");
    cal_curvanture_radius(2);
    printf("        ???????k?????c = ");
    cal_compressive_stress(1, b);
    printf("      ???k???????S??Sc = ");
    cal_safe_rate(2, 1, compression_fatigue);

    printf("???????\t");
    if(safe_rate[0] >= bending_fatigue){
        printf("OK\n");
    }else{
        printf("NG\n");
    }

    printf("???k????\t");
    if(safe_rate[1] >= compression_fatigue){
        printf("OK\n");
    }else{
        printf("NG\n");
    }

    printf("?o?b?N???b?V\t");
    if(backlush_lmit[1] <= backlush && backlush <= backlush_lmit[0]){
        printf("OK\n");
    }else{
        printf("NG\n");
    }
}