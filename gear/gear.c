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
    printf("sokuhi1\t\tsokuhi2\t\tsokuhi3\t\toverall sokuhi\n");
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

    /* Enter mode */
    do
    {
        if(cnt > 0){
            printf("Don't meet the cnditions.\n");
        }
        printf("select mode(0:known_sokuhi, 1:unknown_sokuhi, 2:calculate_error_of_sokuhi)\n");
        scanf("%d", &mode);

    } while (mode < 0);
    printf("\n");
    cnt = 0;
    
    switch (mode)
    {
    case known_sokuhi:
        /* Enter velocity rate */
        do{
            printf("velocity rate cnditions\n");
            printf("velocity rate > 1\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("velocity rate = ");
            scanf("%lf", &sokuhi);

            cnt ++;
        } while( sokuhi <= 1);
        printf("\n");
        cnt = 0;

        /* Enter tolerance range */
        do{
            printf("error range conditions\n");
            printf("error range > 0\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("error range = ");
            scanf("%lf", &error);

            cnt ++;
        } while( error <= 0 );
        printf("\n");
        cnt = 0;

        do{
            printf("Is pinion? yes : 0 , no : 1\n");
            printf("Anser = ");
            scanf("%d", &is_pinion);
        }while(is_pinion > 1);
        printf("\n");

        /* Minimum gear teeth */
        do{
            printf("Minimum gear teeth conditions\n");
            printf("Minimum gear teeth > 6\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("z_min = ");
            scanf("%d", &z_min);

            cnt ++;
        } while( z_min < 6 );
        printf("\n");
        cnt = 0;

        /* Maximum gear teeth */
        do{
            printf("Maximum gear teeth conditions\n");
            printf("Minimum gear teeth > z_min + 1\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("z_max = ");
            scanf("%d", &z_max);

            cnt ++;
        } while( z_max < z_min + 1 );
        printf("\n");
        cnt = 0; 

        printf("----------------------------------------------------------------------------------\n");
        printf("z1\tz2\tz2/z1\t\trate\t\trepeat frequency(Greatest common divisor)\n");
        printf("----------------------------------------------------------------------------------\n");
        processor_base_tooth_combination(sokuhi, error, z_min, z_max, is_pinion);   
        break;

    case unknown_sokuhi:
        /* Enter target velocity rate */
        do
        {
            printf("overall velocity rate cnditions\n");
            printf("overall velocity rate > 1\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("overall velocity rate = ");
            scanf("%lf", &target_sokuhi);

            cnt ++;
        } while (target_sokuhi <= 1);
        printf("\n");
        cnt = 0;

        do
        {
            printf("overall velocity rate error cnditions\n");
            printf("overall velocity rate error > 0\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("overall velocity rate error = ");
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
                printf("using velocity rate %d cnditions\n", i+1);
                printf("using velocity rate %d > 1\n", i+1);

                if(cnt > 0){
                    printf("Don't meet the cnditions.\n");
                }

                printf("using velocity rate %d = ", i+1);
                scanf("%lf", &v[i]);

                cnt ++;
            } while ( v[i] <= 1);
            printf("\n");        
            cnt = 0;
        }

        /* Enter tolerance range */
        do{
            printf("error range conditions\n");
            printf("error range > 0\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("error range = ");
            scanf("%lf", &error);

            cnt ++;
        } while( error <= 0 );
        printf("\n");
        cnt = 0;

        do{
            printf("Is pinion? yes : 0 , no : 1\n");
            printf("Anser = ");
            scanf("%d", &is_pinion);
        }while(is_pinion > 1);
        printf("\n");

        /* Minimum gear teeth */
        do{
            printf("Minimum gear teeth conditions\n");
            printf("Minimum gear teeth > 6\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("z_min = ");
            scanf("%d", &z_min);

            cnt ++;
        } while( z_min < 6 );
        printf("\n");
        cnt = 0;

        /* Maximum gear teeth */
        do{
            printf("Maximum gear teeth conditions\n");
            printf("Minimum gear teeth > z_min + 1\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("z_max = ");
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
            printf("%d stage condition.\n", i+1);
            printf("z1\tz2\tz2/z1\t\trate\t\trepeat frequency(Greatest common divisor)\n");
            printf("----------------------------------------------------------------------------------\n");
            processor_base_tooth_combination(v[i], error, z_min, z_max, is_pinion);   
            printf("\n");
        }
    
        break;

    case calculate_error_of_sokuhi:
        do
        {
            printf("target velocity rate cnditions\n");
            printf("target velocity rate > 1\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("target velocity rate = ");
            scanf("%lf", &target_sokuhi);

            cnt ++;            
        } while (target_error > 1);
        printf("\n");        
        cnt = 0;
        
        do
        {
            printf("target velocity rate error cnditions\n");
            printf("target velocity rate error > 0\n");

            if(cnt > 0){
                printf("Don't meet the cnditions.\n");
            }

            printf("target velocity rate error = ");
            scanf("%lf", &target_error);

            cnt ++;            
        } while (target_error <= 0);
        printf("\n");        
        cnt = 0;

        printf("----------------------------------------------------------------------------------\n");
        printf("----------------------------------------------------------------------------------\n");
        printf("minimum velocity rate\tmaximum velocity rate\n");
        printf("%f\t\t%f\n", (1-target_error)*target_sokuhi, (1+target_error)*target_sokuhi);
        printf("\n");

        break;
    default:
        break;
    }
}