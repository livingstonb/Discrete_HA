
%% HOUSEKEEPING

clear;
close all;

path = '/Users/Brian/Documents/GitHub/MPCrecode';
addpath([path '/Auxiliary Functions']);
cd(path);


%% SPECIFY BASELINE PARAMETERS

% data frequency 
baseline.freq        = 4; % 1 for yearly, 4 for quarterly

% returns
baseline.r           = 0.02;

% demographics
baseline.dieprob     = 1/50;

% preferences
baseline.risk_aver   = 1;
baseline.beta0       = 0.97; % annualized
baseline.temptation  = 0;
baseline.betaL       = 0.80;
% betaH defined in main function file

%warm glow bequests: bequest weight = 0 is accidental
baseline.bequest_weight  = 0; %0.07;
baseline.bequest_luxury  = 0.01; %0.01, must be >0 to avoid NaN error;
baseline.WealthInherited = 0; % 1 for wealth left as bequest, 0 for disappears

% income risk: AR(1) + IID in logs
baseline.LoadIncomeProcess = 0;
baseline.nyT               = 11; %transitory component (not a state variable) (set to 1 for no Transitory Shocks)

% yT,yP (only relevant if LoadIncomeProcess==0)
baseline.NormalizeY   = 1; % 1 to normalize gross income, 0 otherwise
baseline.yTContinuous = 0; % doesn't seem to work properly
baseline.sd_logyT.Q   = sqrt(0.2087);  % (Quarterly) 0.20, relevant if nyT>1
baseline.sd_logyT.A   = sqrt(0.0497);  % (Annual) 0.20, relevant if nyT>1
baseline.lambdaT      = 1; % arrival rate of shocks;
baseline.nyP          = 11; %11 persistent component
baseline.sd_logyP.Q   = sqrt(0.0108); % (Quarterly) 0.1950;
baseline.sd_logyP.A   = sqrt(0.0422); % (Annual)
baseline.rho_logyP.Q  = 0.9881;
baseline.rho_logyP.A  = 0.9525;
baseline.nyF          = 1;
baseline.sd_logyF     = 0;

% cash on hand / savings grid
baseline.nx          = 100;
baseline.xmax        = 1000;  % need high if using high-variance income shocks
baseline.xgrid_par   = 0.3; %1 for linear, 0 for L-shaped
baseline.borrow_lim  = 0;

%government
baseline.labtaxlow       = 0; %proportional tax
baseline.labtaxhigh      = 0; %additional tax on incomes above threshold
baseline.labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in
baseline.savtax          = 0; %0.0001;  %tax rate on savings
baseline.savtaxthresh    = 0; %multiple of mean gross labor income

%discount factor shocks
baseline.nb          = 1; % higher numbers dramatically increase computing load
baseline.betawidth   = 0.001; % too large and eigs hangs while finding stat distribution
baseline.betaswitch  = 1/50; %0;

% computation
baseline.max_iter    = 1e5; % EGP
baseline.tol_iter    = 1.0e-6; % EGP
baseline.Nsim        = 100000; % 100000
baseline.Tsim        = 200;
baseline.nxlong      = 500; % Grid size for final computations
 
% beta iteration
baseline.targetAY    = 3.5;
baseline.maxiterAY   = 40;
baseline.tolAY       = 1e-5;

% mpc options
baseline.Nmpcsim     = 1e6;
baseline.mpcfrac(1)  = -1e-5; %approximate thoeretical mpc
baseline.mpcfrac(2)  = -0.01;
baseline.mpcfrac(3)  = -0.05;
baseline.mpcfrac(4)  = 1e-5; % approximate thoeretical mpc
baseline.mpcfrac(5)  = 0.01; % used in decomposition, 1 percent of average gross labor income: approx $500
baseline.mpcfrac(6)  = 0.05; % 5 percent of average gross labor income: approx $5000

% wealth statistics options
baseline.epsilon = [0 0.005 0.01 0.02 0.05 0.1]; % fraction of mean labor income
baseline.percentiles = [10 25 50 75 90 95 99]; % in percent

% decomposition
baseline.abars = [0 0.01 0.05];

% OPTIONS
baseline.IterateBeta        = 0;
baseline.Display            = 1;
baseline.MakePlots          = 0;
baseline.ComputeDirectMPC   = 1;
baseline.Simulate           = 1;
Batch = 0; % Run alternate parameterizations

%% LOAD ALTERNATE PARAMETERIZATIONS, STRUCTURE ARRAY
if Batch == 0
    params = baseline;
else
    params = parameters(baseline);
end

%% CALL MAIN FUNCTION
Nparams = size(params,2);

direct_results = cell(1,Nparams);
norisk_results = cell(1,Nparams);
sim_results    = cell(1,Nparams);
exceptions     = cell(1,Nparams);
checks         = cell(1,Nparams);

for ip = 1:Nparams
    if Batch == 0
        [SR,DR,NR,checks{ip}] = egp_AR1_IID_tax_recode(params(ip));
        direct_results{ip}  = DR;
        norisk_results{ip}  = NR;
        sim_results{ip}     = SR;      
    else
        try
            % Main function
            [SR,DR,NR,checks{ip}] = egp_AR1_IID_tax_recode(params(ip));
            direct_results{ip}  = DR;
            norisk_results{ip}  = NR;
            sim_results{ip}     = SR;
        catch ME
            checks{ip} = 'EXCEPTION_THROWN';
            exceptions{ip} = ME;
        end
    end
end

%% SAVE
% save('Results','sim_results','direct_results','norisk_results',...
%                                                   'checks','exceptions')