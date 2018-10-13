function params = parameters()

    
    %% set baseline
    
    baseline.name = '';
    
    % data frequency 
    baseline.freq        = 1; % 1 yearly, 4 quarterly

    % returns
    baseline.r           = 0.02;

    % demographics
    baseline.dieprob     = 1/50;

    % preferences
    baseline.risk_aver   = 1;
    baseline.beta0       = 0.98; % annualized
    baseline.temptation  = 0;
    baseline.betaL       = 0.80;
    % betaH defined in main function file

    % warm glow bequests: bequest weight = 0 is accidental
    baseline.bequest_weight  = 0;
    baseline.bequest_luxury  = 0.01; %0.01, must be >0 to avoid NaN error;
    baseline.Bequests = 1; % 1 for wealth left as bequest, 0 for disappears
    baseline.Annuities = 0; % Automatically turns off bequests if set to 1

    % income risk: AR(1) + IID in logs
    baseline.LoadIncomeProcess = 0;
    baseline.nyT               = 11; %transitory component (not a state variable) (set to 1 for no Transitory Shocks)

    % yT,yP (only relevant if LoadIncomeProcess==0)
    baseline.NormalizeY   = 1; % 1 to normalize gross income, 0 otherwise
    baseline.yTContinuous = 0;
    baseline.sd_logyT     = sqrt(0.0497);  % 0.20, relevant if nyT>1
    baseline.lambdaT      = 1; % arrival rate of shocks;
    baseline.nyP          = 11; %11 persistent component
    baseline.sd_logyP     = sqrt(.0422); % 0.1950;
    baseline.rho_logyP    = 0.9525;
    baseline.nyF          = 1;
    baseline.sd_logyF     = 0;

    % cash on hand / savings grid
    baseline.nx          = 150;
    baseline.xmax        = 1000;  % need high if using high-variance income shocks
    baseline.xgrid_par   = 1/3; %1 for linear, 0 for L-shaped
    baseline.borrow_lim  = 0;

    %government
    baseline.labtaxlow       = 0; %proportional tax
    baseline.labtaxhigh      = 0; %additional tax on incomes above threshold
    baseline.labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in
    baseline.savtax          = 0; %0.0001;  %tax rate on savings
    baseline.savtaxthresh    = 0; %multiple of mean gross labor income

    %discount factor shocks
    baseline.nb          = 1; % higher numbers dramatically increase computing load
    baseline.betawidth   = 0.005; % too large and eigs hangs while finding stat distribution
    baseline.betaswitch  = 0; %0;

    % computation
    baseline.max_iter    = 1e5; % EGP
    baseline.tol_iter    = 1.0e-6; % EGP
    baseline.Nsim        = 100000; % 100000
    baseline.Tsim        = 200;
    baseline.nxlong      = 1500; % Grid size for final computations

    % beta iteration
    baseline.targetAY    = 3.5;
    baseline.maxiterAY   = 50;
    baseline.tolAY       = 1e-5;

    % mpc options
    baseline.Nmpcsim     = 1e6; % Number of draws to compute MPCs
    baseline.mpcfrac(1)  = -1e-5; %approximate thoeretical mpc
    baseline.mpcfrac(2)  = -0.01;
    baseline.mpcfrac(3)  = -0.1;
    baseline.mpcfrac(4)  = 1e-5; % approximate thoeretical mpc
    baseline.mpcfrac(5)  = 0.01; % used in decomposition, 1 percent of average gross labor income: approx $500
    baseline.mpcfrac(6)  = 0.1; % 10 percent of average gross labor income: approx $10000

    % wealth statistics options
    baseline.epsilon = [0, 0.005, 0.01, 0.02, 0.05, 0.1]; % fraction of mean labor income
    baseline.percentiles = [10, 25, 50, 75, 90, 95, 99, 99.9]; % in percent

    % decomposition
    baseline.abars = [0, 0.01, 0.05];

    % OPTIONS
    baseline.IterateBeta        = 1;
    baseline.Display            = 0;
    baseline.MakePlots          = 0;
    baseline.ComputeDirectMPC   = 1;
    baseline.Simulate           = 0;
    
    %% create params structure
    
    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    % annual baseline
    params(1)       = baseline;
    params(1).name  = 'Baseline_A';

    % quarterly baseline
    params(2) = baseline;
    params(2).name = 'Baseline_Q';
    params(2).freq = 4;
    params(2).sd_logyP = sqrt(0.0108);
    params(2).sd_logyT = sqrt(0.2087);
    params(2).rho_logyP = 0.9881;

    %----------------------------------------------------------------------
    % PART 2, DIFFERENT ASSUMPTIONS
    %----------------------------------------------------------------------
    
    % different mean wealth targets
    for mw = [0.25, 0.5, 1, 2.5, 5]
        params(end+1) = baseline;
        params(end).name = ['2 AYtarget' num2str(mw)];
        params(end).targetAY = mw;
    end
    
    % different interest rates
    for ii = [0, 1, 3, 5]
        params(end+1) = baseline;
        params(end).name = ['2 IntRate' num2str(ii)];
        params(end).r = ii/100;
    end
    
    % different risk aversion coeffs
    for ira = [0.5, 1.5, 2, 4, 6]
        params(end+1) = baseline;
        params(end).name = ['2 RiskAver' num2str(ira)];
        params(end).risk_aver = ira;
    end
    
    % different tax rates
    for itax = [0.05, 0.1, 0.15, 0.25]
        params(end+1) = baseline;
        params(end).name = ['2 LabTax' num2str(itax)];
        params(end).labtaxlow = itax;
    end
    
    % no death
    params(end+1) = baseline;
    params(end).name = '2 NoDeath';
    params(end).dieprob = 0;
    
    % no bequests
    params(end+1) = baseline;
    params(end).name = '2 NoBequests';
    params(end).Bequests = 0;
    
    % perfect annuities
    params(end+1) = baseline;
    params(end).name = '2 Annuities';
    params(end).Annuities = 1;
    
    % luxury motive...
    
    % fixed beta heterogeneity
    for ibw = [0.001, 0.005, 0.01]
        params(end+1) = baseline;
        params(end).name = ['2 FixedBetaHet_Width' num2str(ibw)];
        params(end).nb = 5;
        params(end).nxlong = 500;
        params(end).betawidth = ibw;
    end
    
    % random beta heterogeneity
    for ibw = [0.001, 0.005, 0.01]
        for bs = [1/10, 1/50]
            params(end+1) = baseline;
            params(end).name = ['2 RandomBetaHet Width' num2str(ibw) ' SwitchProb' num2str(bs)];
            params(end).nb = 5;
            params(end).nxlong = 500;
            params(end).betawidth = ibw;
            params(end).betaswitch = bs;
        end
    end
    
    %----------------------------------------------------------------------
    % PART 3a, ANNUAL MODEL
    %----------------------------------------------------------------------
    
    % i
    params(end+1) = baseline;
    params(end).name = '3a(i) NoTransShocks';
    params(end).nyT = 1;
    params(end).sd_logyT = 0;
    
    % ii
    params(end+1) = baseline;
    params(end).name = '3a(ii) MeasError';
    params(end).sd_logyT = sqrt(0.02);
    
    % iii
    params(end+1) = baseline;
    params(end).name = '3a(iii) NoTranReEst';
    params(end).rho_logyP = 0.8592;
    params(end).sd_logyP = sqrt(0.132);
    params(end).nyT = 1;
    params(end).sd_logyT = 0;

    % iv
    params(end+1) = baseline;
    params(end).name = '3a(iv) HighPersistCarrol';
    params(end).rho_logyP = 0.9995;
    params(end).sd_logyP = sqrt(0.015);
    params(end).sd_logyT = sqrt(0.01);
    
    % v
    params(end+1) = baseline;
    params(end).name = '3a(v) HighPersNotReEst';
    params(end).rho_logyP = 0.99;
    
    % vi
    params(end+1) = baseline;
    params(end).name = '3a(vi) LowPersNotReEst';
    params(end).rho_logyP = 0.9;

    % vii
    params(end+1) = baseline;
    params(end).name = '3a(vii) HighPersReEst';
    params(end).rho_logyP = 0.99;
    params(end).sd_logyP = sqrt(0.0088);
    params(end).sd_logyT = sqrt(0.0667);
    
    % viii
    params(end+1) = baseline;
    params(end).name = '3a(viii) EvenHigherPersReEst';
    params(end).rho_logyP = 0.995;
    params(end).sd_logyP = sqrt(0.0043);
    params(end).sd_logyT = sqrt(0.0688);
    
    % ix
    params(end+1) = baseline;
    params(end).name = '3a(ix) HighPersNoTransReEst';
    params(end).rho_logyP = 0.99;
    params(end).sd_logyP = sqrt(0.0088);
    params(end).nyT = 1;
    params(end).sd_logyT = sqrt(0);
    
    % xi
    params(end+1) = baseline;
    params(end).name = '3a(xi) MatchSS';
    params(end).rho_logyP = 0.9468;
    params(end).sd_logyP = sqrt(0.0641);
    params(end).sd_logyT = sqrt(0.0479);
    params(end).lambdaT  = 0.0821;
    
    % xii
    params(end+1) = baseline;
    params(end).name = '3a(xii) WithSCF m0';
    params(end).rho_logyP = 0.9787;
    params(end).sd_logyP = sqrt(0.0400);
    params(end).sd_logyT = sqrt(0.0508);
    
    % xiv
    params(end+1) = baseline;
    params(end).name = '3a(xiv) MassPointTrans';
    params(end).rho_logyP = sqrt(0.9516);
    params(end).sd_logyP = sqrt(0.0434);
    params(end).sd_logyT = sqrt(0.6431);
    params(end).lambdaT = 0.0760;
    
    %----------------------------------------------------------------------
    % PART 3b, QUARTERLY MODEL
    %----------------------------------------------------------------------
    
    % i
    params(end+1) = baseline;
    params(end).name = '3b(i) KMPTransf';
    params(end).freq = 4;
    params(end).rho_logyP = 0.9879;
    params(end).sd_logyP = sqrt(0.0109);
    params(end).sd_logyT = sqrt(0.0494);
    
    % iv
    params(end+1) = baseline;
    params(end).name = '3b(iv) PersEveryPeriod';
    params(end).freq = 4;
    params(end).rho_logyP = 0.9884;
    params(end).sd_logyP = sqrt(0.0105);
    params(end).sd_logyT = sqrt(1.5298);
    params(end).lambdaT = 0.0813;
    
    %----------------------------------------------------------------------
    % PART 4, Exotic Preferences
    %----------------------------------------------------------------------
    
    % temptation
    for itempt = [0.005 0.01 0.05]
        params(end+1) = baseline;
        params(end).name = ['4 Temptation' num2str(itempt)];
        params(end).temptation = itempt;
    end
end