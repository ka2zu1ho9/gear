void input_int(char str[], int* data){
    printf("%s", str);
    scanf("%d", data);
}

void input_double(char str[], double* data){
    printf("%s", str);
    scanf("%lf", data);
}

void cal_guzai(){
    for (size_t i = 0; i < 4; i++){
        switch ()
        {
        case 0:
            output.guzai[i] = C[0] + C[1]*input.sigma_B;
            printf("%f\n", output.guzai[i]);
            break;

        case 1:
            output.guzai[i] = 1 - exp(-C[2]*input.d);
            printf("%f\n", output.guzai[i]);
            break;

        case 2:
            output.guzai[i] = 1 - exp( (-C[3]*input.d) / input.rowe);
            printf("%f\n", output.guzai[i]);
            break;

        case 3:
            output.guzai[i] = 1 - exp( -C[4]*(1 - input.d / input.D) );
            printf("%f\n", output.guzai[i]);
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
        total_guzai *= output.guzai[i];
    }

    output.notch_coefficient = 1 + total_guzai;

    printf("%f\n", output.notch_coefficient);
    
}


void cal_total_moment_amplitude(){
    output.total_moment_amplitude = sqrt(input.moment_amplitude[0]*input.moment_amplitude[0] + input.moment_amplitude[1]*input.moment_amplitude[1]);
    printf("%f\n", output.total_moment_amplitude);
}


void cal_using_shear_stress(){
    output.using_shear_stress = (16*sqrt( ( (output.notch_coefficient*output.sigma_y_p*output.total_moment_amplitude) / output.sigma_e_B )*( (output.notch_coefficient*output.sigma_y_p*output.total_moment_amplitude) / output.sigma_e_B ) + input.Tav*input.Tav )) / (M_PI*input.d*input.d*input.d);
    printf("%f[kg/mm^2]\n", output.using_shear_stress);
}

void cal_safe_rate(){
    output.safe_rate = output.tau_y_p / output.using_shear_stress;
    printf("%f\n", output.safe_rate);
}

void pro_shaft_strengh(){
    input_int("モードを選んでください。（キー溝: 1　すみ肉:2） = ", &input.mode);
    input_double("曲げモーメント振幅Mr-x[kgf・mm] = ", &input.moment_amplitude[0]);
    input_double("曲げモーメント振幅Mr-x[kgf・mm] = ", &input.moment_amplitude[1]);
    input_double("平均ねじりモーメント[kgf・mm] = ", &input.Tav);
    input_double("引張降伏点[kgf/mm^2] = ", &input.sigma_y_p);
    input_double("回転疲れ限度[kgf/mm^2] = ", &input.sigma_e_B);
    input_double("せん断降伏点[kgf/mm^2] = ", &input.tau_y_p);

    if(input.mode == 1){
        input_double("丸棒直径d[mm] = ", &input.sigma_y_p);
        input_double("切り欠き係数 = ", &input.notch_coefficient);
    }else{
        input_double("引張降強さ[kgf/mm^2] = ", &input.sigma_B);
        input_double("丸棒直径（太）D[mm] = ", &input.D);
        input_double("丸棒直径（細）d[mm] = ", &input.d);
        input_double("すみ肉半径p[mm] = ", &input.rowe);
    }

    cal_guzai();
    cal_notch_coefficient();
    cal_total_moment_amplitude();
    cal_using_shear_stress();
    cal_safe_rate();
}

