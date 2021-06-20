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
    return (tan(angle * (M_PI/180)) - angle);
}

void cal_base_pitch_circle(double module, int z){
    printf("%f\n", (module*z)/2.0 );
}

void cal_engagementpitch_circle(int mode, int z1, int z2, double a){
    if(mode == 1){
        printf("%f\n", (z1 * a) / (z1 + z2) );
    }else if(mode == 2){
        printf("%f\n", (z2 * a) / (z1 + z2) );
    }
}

void cal_engagement_press_angle(int z1, int z2, double center_distance_increase){
    double result;

    result = acos( ( (z1 + z2) / (z1 + z2 + 2*center_distance_increase) )*cos(20*(M_PI/180)) );
    printf("%f\n", result);
}

void cal_backlush_on_pitch_circle(double engagement_press_angle, double module, int z1, int z2){
    double result;

    result = ( cos(20*(M_PI/180)) / cos(engagement_press_angle*(M_PI/180)) )*module*(z1 + z2)*(inv(engagement_press_angle) - inv(20));
    printf("%f\n", result);
}

void cal_tolerance_unit(double base_pitch_circle, double module){
    printf("%f\n", cbrt(base_pitch_circle) + 0.65*module);
}

void cal_center_distance_increase(double module, double a, int z1, int z2){
    printf("%f\n", (a/module) - (z1 + z2)/2 );
}

void cal_tooth_tip_circle(int z, double module, double center_distance_increase, double transfer){
    printf("%f\n", (z + 2) * module + 2*module*(center_distance_increase - transfer));
}

void cal_gear_base_circle(double base_pitch_circle){
    printf("%f\n", (base_pitch_circle/2)*cos(20*(M_PI/180)));
}

void cal_engagement_length(double tooth_tip_radius[], double gear_base_radius[], double a, double press_angle){
    double result;

    result = sqrt( tooth_tip_radius[0]*tooth_tip_radius[0] - gear_base_radius[0]*gear_base_radius[0]) + sqrt( tooth_tip_radius[1]*tooth_tip_radius[1] - gear_base_radius[1]*gear_base_radius[1]);
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
                        printf("%d\t%3d\t%f\t%f\t%f(%f)\n", i, j, sokuhi1, gosa, repeat_frequency, gcd_result);
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
                printf("%d\t%3d\t%f\t%f\t%f(%f)\n", z_min, j, sokuhi1, gosa, repeat_frequency, gcd_result);
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
                printf("%f\t%f\t%f\t%f\n", v[0], v[1], v[2], real_sokuhi);
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
                printf("%f\t%f\t%f\t%f\n", v[0], v[1], v[2], real_sokuhi);
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
                   cal_base_pitch_circle(module[i], z[i]);
                }else{
                   cal_base_pitch_circle(module[i-3], z[i]);
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
        printf("%f\t%f\n", (1-target_error)*target_sokuhi, (1+target_error)*target_sokuhi);
        printf("\n");

        break;
    default:
        break;
    }
}

void pro_gear_strength(){
    int z1, z2;
    double module, a, b, y, n_in, H;
    double bending_fatigue;
    double compression_fatigue;

    input_int("z1 = ", &z1);
    input_int("z2 = ", &z2);
    input_double("module = ", &module);
    input_double("a = ", &a);
    input_double("b = ", &b);
    input_double("y = ", &y);
    input_double("n_in = ", &n_in);
    input_double("H = ", &H);
    input_double("曲げ疲れ限度 = ", &bending_fatigue);
    input_double("圧縮疲れ限度 = ", &compression_fatigue);
}