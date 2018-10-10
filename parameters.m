function params = parameters(baseline)
    baseline.name = '';

    %----------------------------------------------------------------------
    %% BASELINES
    %----------------------------------------------------------------------
    % annual baseline
    params(1)       = baseline;
    params(1).name  = 'Baseline_A';

    % quarterly baseline
    params(end+1)       = baseline;
    params(end+1).name  = 'Baseline_Q';
    params(end+1).freq  = 4;
    params(end+1).sd_logyP = sqrt(0.0108);
    params(end+1).sd_logyT = sqrt(0.2087);
    params(end+1).rho_logyP = 0.9881;
    
    %----------------------------------------------------------------------
    %% PART 3a, ANNUAL MODEL
    %----------------------------------------------------------------------
    
    % i
    params(end+1) = baseline;
    params(end+1).name = 'a_i_NoTransShocks';
    params(end+1).sd_logyT = 0;
    
    % ii
    params(end+1) = baseline;
    params(end+1).name = 'a_ii_MeasError';
    params(end+1).sd_logyT = sqrt(0.02);
    
    % iii
    params(end+1) = baseline;
    params(end+1).name = 'a_iii_NoTrans_ReEst';
    params(end+1).rho_logyP = 0.8592;
    params(end+1).sd_logyP = sqrt(0.132);
    params(end+1).sd_logyT = 0;

    % iv
    params(end+1) = baseline;
    params(end+1).name = 'a_iv_HighPers_Carrol';
    params(end+1).rho_logyP = 0.9995;
    params(end+1).sd_logyP = sqrt(0.015);
    params(end+1).sd_logyT = sqrt(0.01);
    
    % v
    params(end+1) = baseline;
    params(end+1).name = 'a_v_HighPers_NotReEst';
    params(7).rho_logyP = 0.99;
    
    % vi
    params(end+1) = baseline;
    params(end+1).name = 'a_vi_LowPers_NotReEst';
    params(end+1).rho_logyP = 0.9;

    % vii
    params(end+1) = baseline;
    params(end+1).name = 'a_vii_HighPers_ReEst';
    params(end+1).rho_logyP = 0.99;
    params(end+1).sd_logyP = sqrt(0.0088);
    params(end+1).sd_logyT = sqrt(0.0667);
    
    % viii
    params(end+1) = baseline;
    params(end+1).name = 'a_viii_EvenHigherPers_ReEst';
    params(end+1).rho_logyP = 0.995;
    params(end+1).sd_logyP = 0.0043;
    params(end+1).sd_logyT = 0.0688;
    
    % ix
    params(end+1) = baseline;
    params(end+1).name = 'a_ix_HighPersNoTrans_ReEst';
    params(end+1).rho_logyP = 0.99;
    params(end+1).sd_logyP = 0.0088;
    params(end+1).sd_logyT = 0;
    
end