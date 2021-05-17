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
    printf("1�i�ڂ̑���\t2�i�ڂ̑���\t3�i�ڂ̑���\t�S�̂̑���\n");
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

    /* ���[�h����� */
    do
    {
        if(cnt > 0){
            printf("�����ƈ�v���܂���B\n");
        }
        printf("���[�h��I�����Ă��������B(0:1�i�̎��̑g�ݍ��킹�����߂�, 1:�S�̂̑��䂩��A�e�i�̎��̑g�ݍ��킹�����߂�, 2:�ő�̑���ƍŏ��̑�������߂�)\n");
        scanf("%d", &mode);

    } while (mode < 0);
    printf("\n");
    cnt = 0;
    
    switch (mode)
    {
    case known_sokuhi:
        /* �������� */
        do{
            printf("����̏���\n");
            printf("���� > 1\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("���� = ");
            scanf("%lf", &sokuhi);

            cnt ++;
        } while( sokuhi <= 1);
        printf("\n");
        cnt = 0;

        /* ���z���䂩��̋��e�͈͂���� */
        do{
            printf("���z���䂩��̋��e�͈͂̏���\n");
            printf("���z���䂩��̋��e�͈� > 0\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("���z���䂩��̋��e�͈� = ");
            scanf("%lf", &error);

            cnt ++;
        } while( error <= 0 );
        printf("\n");
        cnt = 0;

        do{
            printf("�s�j�I���͂���܂���? yes : 0 , no : 1\n");
            printf("���� = ");
            scanf("%d", &is_pinion);
        }while(is_pinion > 1);
        printf("\n");

        /* �ŏ���������� */
        do{
            printf("�ŏ������̏���\n");
            printf("Minimum gear teeth >= 6\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("�ŏ����� = ");
            scanf("%d", &z_min);

            cnt ++;
        } while( z_min < 6 );
        printf("\n");
        cnt = 0;

        /* �ő厕���̓��� */
        do{
            printf("�ő厕���̏���\n");
            printf("�ő厕�� > �ŏ����� + 1\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("�ő厕�� = ");
            scanf("%d", &z_max);

            cnt ++;
        } while( z_max < z_min + 1 );
        printf("\n");
        cnt = 0; 

        printf("----------------------------------------------------------------------------------\n");
        printf("z1\tz2\tz2/z1\t\trate\t\t�J��Ԃ��p�x(�ő����)\n");
        printf("----------------------------------------------------------------------------------\n");
        processor_base_tooth_combination(sokuhi, error, z_min, z_max, is_pinion);   
        break;

    case unknown_sokuhi:
        /* �S�̗̂��z�������� */
        do
        {
            printf("�S�̗̂��z����\n");
            printf("�S�̗̂��z���� > 1\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("�����S�̗̂��z���� = ");
            scanf("%lf", &target_sokuhi);

            cnt ++;
        } while (target_sokuhi <= 1);
        printf("\n");
        cnt = 0;

        do
        {
            printf("�S�̗̂��z���䂩��̋��e�͈�\n");
            printf("�S�̗̂��z���䂩��̋��e�͈� > 0\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("�S�̗̂��z���䂩��̋��e�͈�= ");
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
                printf("�g�p����%d�i�ڂ̑���̏���\n", i+1);
                printf("�g�p����%d�i�ڂ̑��� > 1\n", i+1);

                if(cnt > 0){
                    printf("��������v���܂���B\n");
                }

                printf("�g�p����%d�i�ڂ̑��� = ", i+1);
                scanf("%lf", &v[i]);

                cnt ++;
            } while ( v[i] <= 1);
            printf("\n");        
            cnt = 0;
        }

        /* ���z���䂩��̋��e�͈͂���� */
        do{
            printf("���z���䂩��̋��e�͈͂̏���\n");
            printf("���z���䂩��̋��e�͈� > 0\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("���z���䂩��̋��e�͈� = ");
            scanf("%lf", &error);

            cnt ++;
        } while( error <= 0 );
        printf("\n");
        cnt = 0;

        do{
            printf("�s�j�I���͂���܂���? yes : 0 , no : 1\n");
            printf("���� = ");
            scanf("%d", &is_pinion);
        }while(is_pinion > 1);
        printf("\n");

        /* �ŏ���������� */
        do{
            printf("�ŏ������̏���\n");
            printf("Minimum gear teeth >= 6\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("�ŏ����� = ");
            scanf("%d", &z_min);

            cnt ++;
        } while( z_min < 6 );
        printf("\n");
        cnt = 0;

        /* �ő厕���̓��� */
        do{
            printf("�ő厕���̏���\n");
            printf("�ő厕�� > �ŏ����� + 1\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("�ő厕�� = ");
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
            printf("%d�i�ڂ̏���\n", i+1);
            printf("z1\tz2\tz2/z1\t\trate\t\t�J��Ԃ��p�x(�ő����)\n");
            printf("----------------------------------------------------------------------------------\n");
            processor_base_tooth_combination(v[i], error, z_min, z_max, is_pinion);   
            printf("\n");
        }
    
        break;

    case calculate_error_of_sokuhi:
        do
        {
            printf("���z����̏���\n");
            printf("���z���� > 1\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("���z���� = ");
            scanf("%lf", &target_sokuhi);

            cnt ++;            
        } while (target_error > 1);
        printf("\n");        
        cnt = 0;
        
        do
        {
            printf("���z���䂩��̋��e�͈͂̏���\n");
            printf("���z���䂩��̋��e�͈� > 0\n");

            if(cnt > 0){
                printf("��������v���܂���B\n");
            }

            printf("���z���䂩��̋��e�͈� = ");
            scanf("%lf", &target_error);

            cnt ++;            
        } while (target_error <= 0);
        printf("\n");        
        cnt = 0;

        printf("----------------------------------------------------------------------------------\n");
        printf("----------------------------------------------------------------------------------\n");
        printf("�ŏ�����\t�ő呬��\n");
        printf("%f\t%f\n", (1-target_error)*target_sokuhi, (1+target_error)*target_sokuhi);
        printf("\n");

        break;
    default:
        break;
    }
}