function params = parameters(baseline)
    baseline.name = '';

    %% PARAMETERIZATION 1 - BASELINE (ANNUAL)
    params(1)       = baseline;
    params(1).name  = 'Baseline_A';
    params(1).freq  = 1;


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
    params(3).freq = 1;
    params(3).sd_logyT = 0;


end