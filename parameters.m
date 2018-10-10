function params = parameters(baseline)
    baseline.name = '';

    %% PARAMETERIZATION 1 - BASELINE (ANNUAL)
    params(1)       = baseline;
    params(1).name  = 'Baseline_A';

    %% PARAMETERIZATION 2 - BASELINE (QUARTERLY)
    params(2)       = baseline;
    params(2).name  = 'Baseline_Q';
    params(2).freq  = 4;
    params(2).sd_logyP = sqrt(0.0108);
    params(2).sd_logyT = sqrt(0.2087);
    params(2).rho_logyP = 0.9881;
    
    %% PARAMETERIZATION 3
    params(3) = baseline;
    params(3).name = 'i_NoTransShocks';
    params(3).sd_logyT = 0;
    
    %% PARAMETERIZATION 4
    params(4) = baseline;
    params(4).name = 'ii_MeasError';
    params(4).sd_logyT = sqrt(0.02);
    
    %% PARAMETERIZATION 5
    params(5) = baseline;
    params(5).name = 'iii_NoTrans_ReEst';
    params(5).rho_logyP = 0.8592;
    params(5).sd_logyP = sqrt(0.132);
    params(5).sd_logyT = 0;

    %% PARAMETERIZATION 6
    params(6) = baseline;
    params(6).name = 'iv_HighPers_Carrol';
    params(6).rho_logyP = 0.9995;
    params(6).sd_logyP = sqrt(0.015);
    params(6).sd_logyT = sqrt(0.01);
    
    %% PARAMETERIZATION 7
    params(7) = baseline;
    params(7).name = 'v_HighPers_NotReEst';
    params(7).rho_logyP = 0.99;
    params(7).sd_logyP = sqrt(0.0422);
    params(7).sd_logyT = sqrt(0.0497);

end