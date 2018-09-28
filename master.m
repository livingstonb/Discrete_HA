
%% Housekeeping

clear;
close all;

path = '/Users/brianlivingston/Documents/GitHub/MPCrecode';
addpath([path '/Auxiliary Functions']);
cd(path);


%% Set parameters
prms = struct();

% data frequency
prms(1).freq        = 1; % 1 for yearly, 4 for quarterly

% returns
prms(1).r           = 0.02;
 
% demographics
prms(1).dieprob     = 1/50;

% preferences
prms(1).risk_aver   = 1;
prms(1).beta0       = 0.97;
prms(1).temptation  = 0;
prms(1).betaL       = 0.80;
% betaH defined in main function file

%warm glow bequests: bequessgrt_weight = 0 is accidental
prms(1).bequest_weight = 0; %0.07;
prms(1).bequest_luxury = 0.01; %0.01, must be >0 to avoid NaN error;
prms(1).WealthInherited = 1; % 1 for wealth left as bequest, 0 for disappears

% income risk: AR(1) + IID in logs
prms(1).LoadIncomeProcess   = 0;
prms(1).nyT                 = 11; %transitory component (not a state variable) (set to 1 for no Transitory Shocks)

% income normalization
prms(1).NormalizeY  = 1;

%only relevant if LoadIncomeProcess==0
prms(1).yTContinuous = 0;
prms(1).sd_logyT    = sqrt(0.2);  %0.20; %relevant if nyT>1
prms(1).lambdaT     = 1; % arrival rate of shocks;
prms(1).nyP         = 11; %11 persistent component
prms(1).sd_logyP    = sqrt(0.2); %0.1950;
prms(1).rho_logyP   = 0.9525;
prms(1).nyF         = 1;
prms(1).sd_logyF    = 0;

% cash on hand / savings grid
prms(1).nx          = 500;
prms(1).xmax        = 1000;  % need high if using high-variance income shocks
prms(1).xgrid_par   = 0.2; %1 for linear, 0 for L-shaped
prms(1).borrow_lim  = 0;

%government
prms(1).labtaxlow       = 0; %proportional tax
prms(1).labtaxhigh      = 0; %additional tax on incomes above threshold
prms(1).labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in

prms(1).savtax          = 0; %0.0001;  %tax rate on savings
prms(1).savtaxthresh    = 0; %multiple of mean gross labor income

%discount factor shocks;
prms(1).nb          = 2;  %1 or 2
prms(1).betawidth   = 0.02; % beta +/- beta width
prms(1).betaswitch  = 1/50; %0;

% computation
prms(1).max_iter    = 1e5; % for EGP
prms(1).tol_iter    = 1.0e-6;
prms(1).Nsim        = 100000;
prms(1).Tsim        = 200;
prms(1).nxlong_int  = 200; % Grid size for computing intermediate ergodic distributions
prms(1).nxlong      = 500; % Grid size for computing final ergodic distribution
prms(1).nxmpc       = 1000; % grid size for MPC. mpcamount = xmax/xmpc
 
prms(1).targetAY    = 3.5;
prms(1).maxiterAY   = 20;
prms(1).tolAY       = 1.0e-4;
prms(1).FastIter    = 1; % Use routine with low tolerance, small grid, until closer to target

%mpc options
prms(1).mpcfrac{1}  = -0.01; %approximate thoeretical mpc
prms(1).mpcfrac{2}  = -0.001; % 1 percent of average gross labor income: approx $500
prms(1).mpcfrac{3}  = -1e-5; % 5 percent of average gross labor income: approx $5000
prms(1).mpcfrac{4}  = 1e-5; % approximate thoeretical mpc
prms(1).mpcfrac{5}  = 0.01; % 1 percent of average gross labor income: approx $500
prms(1).mpcfrac{6}  = 0.05; % 5 percent of average gross labor income: approx $5000

% OPTIONS
prms(1).IterateBeta         = 0;
prms(1).Display             = 1;
prms(1).MakePlots           = 0;
prms(1).ComputeDistMPC      = 1;
prms(1).ComputeSimMPC       = 1;
prms(1).SolveDeterministic  = 0;
prms(1).Simulate            = 1;
prms(1).PrintStats          = 1;

%% Call model
Nprms = size(prms,2);
% Create structure arrays to store results
for ip = 1:Nprms
    [simulations(ip),results(ip)] = egp_AR1_IID_tax_recode(prms(ip));
end
