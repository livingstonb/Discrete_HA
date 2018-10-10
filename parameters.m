function params = parameters(baseline)
    params = ([]);
    baseline.name = 'baseline_A';

    %% PARAMETERIZATION 1 - BASELINE (ANNUAL)
    params(1)       = baseline;
    params(1).name  = 'baseline_A';

    %% PARAMETERIZATION 2 - BASELINE (QUARTERLY)
    params(2)       = baseline;
    params(2).name  = 'baseline_Q';
    params(2).freq  = 4;

    %% PARAMETERIZATION 3
    params(3) = baseline;


end