function params = parameters(baseline)
    baseline.name = '';

    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    % annual baseline
    params(1)       = baseline;
    params(1).name  = 'Baseline_A';

    % quarterly baseline
    params(end+1) = baseline;
    params(end).name = 'Baseline_Q';
    params(end).freq = 4;
    params(end).sd_logyP = sqrt(0.0108);
    params(end).sd_logyT = sqrt(0.2087);
    params(end).rho_logyP = 0.9881;
    
    %----------------------------------------------------------------------
    % PART 2, DIFFERENT ASSUMPTIONS
    %----------------------------------------------------------------------
    
    % different mean wealth targets
    for mw = [0.25, 0.5, 1, 2.5, 5]
        params(end+1) = baseline;
        params(end).name = ['AYtarget' num2str(mw)];
        params(end).targetAY = mw;
    end
    
    % different interest rates
    for ii = [0, 1, 3, 5]
        params(end+1) = baseline;
        params(end).name = ['IntRate' num2str(ii)];
        params(end).r = ii/100;
    end
    
    % different risk aversion coeffs
    for ira = [0.5, 1.5, 2, 4, 6]
        params(end+1) = baseline;
        params(end).name = ['RiskAver' num2str(ira)];
        params(end).risk_aver = ira;
    end
    
    % different tax rates
    for itax = [0.05, 0.1, 0.15, 0.25]
        params(end+1) = baseline;
        params(end).name = ['LabTax' num2str(itax)];
        params(end).labtaxlow = itax;
    end
    
    % no death
    params(end+1) = baseline;
    params(end).name = 'NoDeath';
    params(end).dieprob = 0;
    
    % no bequests
    params(end+1) = baseline;
    params(end).name = 'NoBequests';
    params(end).Bequests = 0;
    
    % perfect annuities
    params(end+1) = baseline;
    params(end).name = 'Annuities';
    params(end).r = r + params(end).dieprob;
    
    % luxury motive...
    
    %
    
    %----------------------------------------------------------------------
    % PART 3a, ANNUAL MODEL
    %----------------------------------------------------------------------
    
    % i
    params(end+1) = baseline;
    params(end).name = 'a_i_NoTransShocks';
    params(end).sd_logyT = 0;
    
    % ii
    params(end+1) = baseline;
    params(end).name = 'a_ii_MeasError';
    params(end).sd_logyT = sqrt(0.02);
    
    % iii
    params(end+1) = baseline;
    params(end).name = 'a_iii_NoTrans_ReEst';
    params(end).rho_logyP = 0.8592;
    params(end).sd_logyP = sqrt(0.132);
    params(end).sd_logyT = 0;

    % iv
    params(end+1) = baseline;
    params(end).name = 'a_iv_HighPers_Carrol';
    params(end).rho_logyP = 0.9995;
    params(end).sd_logyP = sqrt(0.015);
    params(end).sd_logyT = sqrt(0.01);
    
    % v
    params(end+1) = baseline;
    params(end).name = 'a_v_HighPers_NotReEst';
    params(end).rho_logyP = 0.99;
    
    % vi
    params(end+1) = baseline;
    params(end).name = 'a_vi_LowPers_NotReEst';
    params(end).rho_logyP = 0.9;

    % vii
    params(end+1) = baseline;
    params(end).name = 'a_vii_HighPers_ReEst';
    params(end).rho_logyP = 0.99;
    params(end).sd_logyP = sqrt(0.0088);
    params(end).sd_logyT = sqrt(0.0667);
    
    % viii
    params(end+1) = baseline;
    params(end).name = 'a_viii_EvenHigherPers_ReEst';
    params(end).rho_logyP = 0.995;
    params(end).sd_logyP = 0.0043;
    params(end).sd_logyT = 0.0688;
    
    % ix
    params(end+1) = baseline;
    params(end).name = 'a_ix_HighPersNoTrans_ReEst';
    params(end).rho_logyP = 0.99;
    params(end).sd_logyP = 0.0088;
    params(end).sd_logyT = 0;
    
end