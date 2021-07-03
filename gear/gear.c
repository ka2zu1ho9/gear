#include "gear.h"

void input_int(char str[], int* data){
    printf("%s", str);
    scanf("%d", data);
}

void input_double(char str[], double* data){
    printf("%s", str);
    scanf("%lf", data);
}

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
    double result;

    result = ( cos(20*(M_PI/180.0)) / cos(engagement_press_angle*(M_PI/180.0)) )*module*(z1 + z2)*(inv(engagement_press_angle) - inv(20));
    printf("%.4f[mm]\n", result);
}

void cal_tolerance_unit(double mode, double module){
    if(mode == 1){
        printf("%.4f[μm]\n", cbrt(base_pitch_circle[0]) + 0.65*module);
    }else if(mode == 2){
        printf("%.4f[μm]\n", cbrt(base_pitch_circle[1]) + 0.65*module);
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
        printf("%.4f[kgf・mm]\n", shaft_torque[0]);
    }else{
        shaft_torque[1] = ((75 * 60 * H) / (2*M_PI*n_in));
        printf("%.4f[kgf・mm]\n", shaft_torque[1]);
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
        // 曲げ
        if(mode == 1){
            safe_rate[0] = ( ( (45*pow(10, 6)) / (9.8*pow(10, 6)) ) / bending_stress[0] );
            printf("%.4f\n",  safe_rate[0]);
        }else{
            safe_rate[0] = ( ( (45*pow(10, 6)) / (9.8*pow(10, 6)) ) / bending_stress[1]);
            printf("%.4f\n", safe_rate[0]);
        }
    }else{
        // 圧縮
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
    printf("1段目の速比\t2段目の速比\t3段目の速比\t全体の速比\n");
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
            printf("条件と一致しません。\n");
        }
        printf("歯数1 = \n");
        scanf("%d", &z[0]);

    } while (z[0] < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("条件と一致しません。\n");
        }
        printf("歯数2 = \n");
        scanf("%d", &z[1]);

    } while (z[1] < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("条件と一致しません。\n");
        }
        printf("モジュール = \n");
        scanf("%d", &module);

    } while (module < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("条件と一致しません。\n");
        }
        printf("中心距離[mm] = \n");
        scanf("%d", &a);

    } while (a < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("条件と一致しません。\n");
        }
        printf("歯幅[mm] = \n");
        scanf("%d", &a);

    } while (b < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("条件と一致しません。\n");
        }
        printf("歯数係数 = \n");
        scanf("%d", &gear_coefficient);

    } while (gear_coefficient < 0);
    printf("\n");
    cnt = 0;            

    do
    {
        if(cnt > 0){
            printf("条件と一致しません。\n");
        }
        printf("歯数係数 = \n");
        scanf("%d", &n_in);

    } while (n_in < 0);
    printf("\n");
    cnt = 0;

    do
    {
        if(cnt > 0){
            printf("条件と一致しません。\n");
        }
        printf("動力伝達係数 = \n");
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

    /* モードを入力 */
    do
    {
        if(cnt > 0){
            printf("条件と一致しません。\n");
        }
        printf("モードを選択してください。(0:1段の歯の組み合わせを求める, 1:全体の速比から、各段の歯の組み合わせを求める, 2:最大の速比と最小の速比を求める)\n");
        scanf("%d", &mode);

    } while (mode < 0);
    printf("\n");
    cnt = 0;
    
    switch (mode)
    {
    case known_sokuhi:
        /* 速比を入力 */
        do{
            printf("速比の条件\n");
            printf("速比 > 1\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("速比 = ");
            scanf("%lf", &sokuhi);

            cnt ++;
        } while( sokuhi <= 1);
        printf("\n");
        cnt = 0;

        /* 理想速比からの許容範囲を入力 */
        do{
            printf("理想速比からの許容範囲の条件\n");
            printf("理想速比からの許容範囲 > 0\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("理想速比からの許容範囲 = ");
            scanf("%lf", &error);

            cnt ++;
        } while( error <= 0 );
        printf("\n");
        cnt = 0;

        do{
            printf("ピニオンはありますか? yes : 0 , no : 1\n");
            printf("答え = ");
            scanf("%d", &is_pinion);
        }while(is_pinion > 1);
        printf("\n");

        /* 最小歯数を入力 */
        do{
            printf("最小歯数の条件\n");
            printf("Minimum gear teeth >= 6\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("最小歯数 = ");
            scanf("%d", &z_min);

            cnt ++;
        } while( z_min < 6 );
        printf("\n");
        cnt = 0;

        /* 最大歯数の入力 */
        do{
            printf("最大歯数の条件\n");
            printf("最大歯数 > 最小歯数 + 1\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("最大歯数 = ");
            scanf("%d", &z_max);

            cnt ++;
        } while( z_max < z_min + 1 );
        printf("\n");
        cnt = 0; 

        printf("----------------------------------------------------------------------------------\n");
        printf("z1\tz2\tz2/z1\t\t誤差\t\t繰り返し頻度(最大公約数)\n");
        printf("----------------------------------------------------------------------------------\n");
        pro_base_tooth_combination(sokuhi, error, z_min, z_max, is_pinion);   
        break;

    case unknown_sokuhi:
        /* 全体の理想速比を入力 */
        do
        {
            printf("全体の理想速比\n");
            printf("全体の理想速比 > 1\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("条件全体の理想速比 = ");
            scanf("%lf", &target_sokuhi);

            cnt ++;
        } while (target_sokuhi <= 1);
        printf("\n");
        cnt = 0;

        do
        {
            printf("全体の理想速比からの許容範囲\n");
            printf("全体の理想速比からの許容範囲 > 0\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("全体の理想速比からの許容範囲= ");
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
                printf("使用する%d段目の速比の条件\n", i+1);
                printf("使用する%d段目の速比 > 1\n", i+1);

                if(cnt > 0){
                    printf("条件が一致しません。\n");
                }

                printf("使用する%d段目の速比 = ", i+1);
                scanf("%lf", &v[i]);

                cnt ++;
            } while ( v[i] <= 1);
            printf("\n");        
            cnt = 0;
        }

        /* 理想速比からの許容範囲を入力 */
        do{
            printf("理想速比からの許容範囲の条件\n");
            printf("理想速比からの許容範囲 > 0\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("理想速比からの許容範囲 = ");
            scanf("%lf", &error);

            cnt ++;
        } while( error <= 0 );
        printf("\n");
        cnt = 0;

        do{
            printf("ピニオンはありますか? yes : 0 , no : 1\n");
            printf("答え = ");
            scanf("%d", &is_pinion);
        }while(is_pinion > 1);
        printf("\n");

        /* 最小歯数を入力 */
        do{
            printf("最小歯数の条件\n");
            printf("Minimum gear teeth >= 6\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("最小歯数 = ");
            scanf("%d", &z_min);

            cnt ++;
        } while( z_min < 6 );
        printf("\n");
        cnt = 0;

        /* 最大歯数の入力 */
        do{
            printf("最大歯数の条件\n");
            printf("最大歯数 > 最小歯数 + 1\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("最大歯数 = ");
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
            printf("%d段目の条件\n", i+1);
            printf("z1\tz2\tz2/z1\t\t誤差\t\t繰り返し頻度(最大公約数)\n");
            printf("----------------------------------------------------------------------------------\n");
            pro_base_tooth_combination(v[i], error, z_min, z_max, is_pinion);   
            printf("\n");
        }

        do{
            if(cnt > 0){
                printf("条件が一致しません。\n");
            }
            
            printf("かみ合いピッチ円を計算しますか？y : 0, n : 1\n"); 
            scanf("%d", &do_cal_engagement_pitch_circle);
            cnt ++;
        } while(do_cal_engagement_pitch_circle > 1);
        printf("\n");
        cnt = 0; 

        if(do_cal_engagement_pitch_circle == 0){
            for(int i = 0; i < 3; i++){
                printf("選択したmoduleを一段目から順に入力してください。\n");
                printf("module[%d] = ", i);
                scanf("%lf", &module[i]);
            }

            for(int i = 0; i < 6; i++){
                printf("選択した歯数を一段目から順に入力してください。\n");
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
            printf("理想速比の条件\n");
            printf("理想速比 > 1\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("理想速比 = ");
            scanf("%lf", &target_sokuhi);

            cnt ++;            
        } while (target_error > 1);
        printf("\n");        
        cnt = 0;
        
        do
        {
            printf("理想速比からの許容範囲の条件\n");
            printf("理想速比からの許容範囲 > 0\n");

            if(cnt > 0){
                printf("条件が一致しません。\n");
            }

            printf("理想速比からの許容範囲 = ");
            scanf("%lf", &target_error);

            cnt ++;            
        } while (target_error <= 0);
        printf("\n");        
        cnt = 0;

        printf("----------------------------------------------------------------------------------\n");
        printf("----------------------------------------------------------------------------------\n");
        printf("最小速比\t最大速比\n");
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
    input_double("歯型係数 = ", &gear_fatigue);
    input_double("n_in = ", &n_in);
    input_double("H = ", &H);
    input_double("曲げ応力の安全率（下限） = ", &bending_fatigue);
    input_double("圧縮応力の安全率（下限） = ", &compression_fatigue);

    printf("\n");

    printf("\t\t   速比u = %.4f\n", (double)z2 / (double)z1);
    printf("     基準ピッチ円直径d01 = ");
    cal_base_pitch_circle(1, module, z1);
    printf("     基準ピッチ円直径d02 = ");
    cal_base_pitch_circle(2, module, z2);
    printf("   かみ合いピッチ半径rb1 = ");
    cal_engagement_pitch_circle(1, z1, z2, a);
    printf("   かみ合いピッチ半径rb2 = ");
    cal_engagement_pitch_circle(2, z1, z2, a);    
    printf("\t      公差単位W1 = ");
    cal_tolerance_unit(1, module);
    printf("\t      公差単位W2 = ");
    cal_tolerance_unit(2, module);    
    printf("       中心距離増加係数y = ");
    cal_center_distance_increase(module, a, z1, z2);
    printf("\t    かみ合い角ab = ");
    cal_engagement_press_angle(z1, z2, a);
    printf("ピッチ円上バックラッシC0 = ");
    cal_backlush_on_pitch_circle(module, z1, z2);
    printf("\t   刃先円直径dk1 = ");
    cal_tooth_tip_circle(1, z1, module, y, 0);
    printf("\t   刃先円直径dk2 = ");
    cal_tooth_tip_circle(2, z2,module, y, 0);
    printf("\t   基礎円半径rg1 = ");
    cal_gear_base_circle(1);
    printf("\t   基礎円半径rg2 = ");
    cal_gear_base_circle(2);    
    printf("\t    かみ合い長さ = ");
    cal_engagement_length(a);
    printf("\t    法線ピッチte = ");
    cal_normal_pitch(module);
    printf("\t     かみ合い率ε = ");
    cal_engagement_rate();
    printf("\t    軸のトルクTa = ");
    cal_shaft_torque(1, n_in, H);
    printf("\t    垂直作用力Fn = ");
    cal_varitical_force(1);
    printf("\t        接線力P0 = ");
    cal_tanquent_force(1);
    printf("       歯車1の曲げ応力σ1 = ");
    cal_bending_stress(1, module, b, gear_fatigue);
    printf("\t 歯車1の安全率Sb = ");
    cal_safe_rate(1, 1, bending_fatigue);
    printf("       歯車1の曲率半径p1 = ");
    cal_curvanture_radius(1);
    printf("       歯車2の曲率半径p1 = ");
    cal_curvanture_radius(2);
    printf("        歯車の圧縮応力σc = ");
    cal_compressive_stress(1, b);
    printf("      圧縮応力の安全率Sc = ");
    cal_safe_rate(2, 1, compression_fatigue);

    printf("曲げ応力\t");
    if(safe_rate[0] >= bending_fatigue){
        printf("OK\n");
    }else{
        printf("NG\n");
    }

    printf("圧縮応力\t");
    if(safe_rate[1] >= compression_fatigue){
        printf("OK\n");
    }else{
        printf("NG\n");
    }
}