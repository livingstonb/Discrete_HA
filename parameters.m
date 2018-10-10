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
    params(3).name = 'a_i_NoTransShocks';
    params(3).sd_logyT = 0;
    
    %% PARAMETERIZATION 4
    params(4) = baseline;
    params(4).name = 'a_ii_MeasError';
    params(4).sd_logyT = sqrt(0.02);
    
    %% PARAMETERIZATION 5
    params(5) = baseline;
    params(5).name = 'a_iii_NoTrans_ReEst';
    params(5).rho_logyP = 0.8592;
    params(5).sd_logyP = sqrt(0.132);
    params(5).sd_logyT = 0;

    %% PARAMETERIZATION 6
    params(6) = baseline;
    params(6).name = 'a_iv_HighPers_Carrol';
    params(6).rho_logyP = 0.9995;
    params(6).sd_logyP = sqrt(0.015);
    params(6).sd_logyT = sqrt(0.01);
    
    %% PARAMETERIZATION 7
    params(7) = baseline;
    params(7).name = 'a_v_HighPers_NotReEst';
    params(7).rho_logyP = 0.99;
    
    %% PARAMATERIZATION 8
    params(8) = baseline;
    params(8).name = 'a_vi_LowPers_NotReEst';
    params(8).rho_logyP = 0.9;

    %% PARAMETERIZATION 9
    params(9) = baseline;
    params(9).name = 'a_vii_HighPers_ReEst';
    params(9).rho_logyP = 0.99;
    params(9).sd_logyP = sqrt(0.0088);
    params(9).sd_logyT = sqrt(0.0667);
    
    %% PARAMETERIZATION 10
    params(10) = baseline;
    params(10).name = 'a_viii_EvenHigherPers_ReEst';
    params(10).rho_logyP = 0.995;
    params(10).sd_logyP = 0.0043;
    params(10).sd_logyT = 0.0688;
    
    %% PARAMETERIZATION 11
    params(11) = baseline;
    params(11).name = 'a_ix_HighPersNoTrans_ReEst';
    params(11).rho_logyP = 0.99;
    params(11).sd_logyP = 0.0088;
    params(11).sd_logyT = 0;
    
end