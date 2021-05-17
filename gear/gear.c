#include "gear.h"

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


void processor_base_tooth_combination(
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
                        printf("%d\t%3d\t%f\t%f\t%f(%f)\n", i, j, sokuhi1, rate, repeat_frequency, gcd_result);
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
                printf("%d\t%3d\t%f\t%f\t%f(%f)\n", z_min, j, sokuhi1, rate, repeat_frequency, gcd_result);
            }
        }   
    }
}

void processor_other_tooth_combination(double target_sokuhi, double target_error){
    double v[3];
    double real_sokuhi;
    double real_error;

    v[0] = pow(target_sokuhi, 1.0/3.0);

    printf("-----------------------------------------------------------------------------------------------\n");
    printf("1段目の速比\t2段目の速比\t3段目の速比\t全体の速比\n");
    printf("-----------------------------------------------------------------------------------------------\n");

    for(double k = -1.0; k < 1.0; k += 0.005){
        v[1] = v[0] + k;
        v[2] = v[1] + 2*k;
        real_sokuhi = v[0]*v[1]*v[2];

        if( ( (1-target_error)*target_sokuhi < real_sokuhi ) && ( real_sokuhi < (1+target_error)*target_sokuhi ) ){
            real_error = ( (real_sokuhi - target_sokuhi) / target_sokuhi )*100;
            printf("%f\t%f\t%f\t%f\n", v[0], v[1], v[2], real_sokuhi);
        }
    }
}

void calculate_tooth_combination(){
    double sokuhi = 0;
    double error = 0;

    double v[3];
    double target_sokuhi = 0;
    double target_error = 0;

    int z_min = 0;
    int z_max = 0;
    int is_pinion = 2;
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
        printf("z1\tz2\tz2/z1\t\trate\t\t繰り返し頻度(最大公約数)\n");
        printf("----------------------------------------------------------------------------------\n");
        processor_base_tooth_combination(sokuhi, error, z_min, z_max, is_pinion);   
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

        processor_other_tooth_combination(target_sokuhi, target_error);
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
            printf("z1\tz2\tz2/z1\t\trate\t\t繰り返し頻度(最大公約数)\n");
            printf("----------------------------------------------------------------------------------\n");
            processor_base_tooth_combination(v[i], error, z_min, z_max, is_pinion);   
            printf("\n");
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