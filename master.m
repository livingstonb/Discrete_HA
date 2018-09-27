
%% Housekeeping

clear;
close all;

path = '/Users/Brian/Documents/GitHub/MPCrecode';
addpath([path '/Auxiliary Functions']);
cd(path);

% For MPC computations, half of grid points are in bottom 1% of asset
% space, so pick size of asset space according to formula:

% na = 0.9*xmax/(0.5*mpcamount)

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
prms(1).beta0       = 0.97365;
prms(1).temptation  = 0;
prms(1).betaL       = 0.90;
% betaH defined in main function file

%warm glow bequests: bequessgrt_weight = 0 is accidental
prms(1).bequest_weight = 0; %0.07;
prms(1).bequest_luxury = 0.01; %0.01, must be >0 to avoid NaN error;
prms(1).WealthInherited = 1; % 1 for wealth left as bequest, 0 for disappears

% income risk: AR(1) + IID in logs
prms(1).LoadIncomeProcess   = 1;
prms(1).nyT                 = 11; %transitory component (not a state variable) (set to 1 for no Transitory Shocks)

% income normalization
prms(1).NormalizeY  = 1;

%only relevant if LoadIncomeProcess==0
prms(1).yTContinuous = 0;
prms(1).sd_logyT    = sqrt(0.2);  %0.20; %relevant if nyT>1
prms(1).lambdaT     = 1; % arrival rate of shocks;
prms(1).nyP         = 11; %11 persistent component
prms(1).sd_logyP    = sqrt(0.1950); %0.1950;
prms(1).rho_logyP   = 0.9947;
prms(1).nyF         = 1;
prms(1).sd_logyF    = 0;

% cash on hand / savings grid
prms(1).nx          = 50;
prms(1).xmax        = 40;  %multiple of mean gross labor income
prms(1).xgrid_par   = 0.3; %1 for linear, 0 for L-shaped
prms(1).borrow_lim  = 0;

%government
prms(1).labtaxlow       = 0; %proportional tax
prms(1).labtaxhigh      = 0; %additional tax on incomes above threshold
prms(1).labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in

prms(1).savtax          = 0; %0.0001;  %tax rate on savings
prms(1).savtaxthresh    = 0; %multiple of mean gross labor income

%discount factor shocks;
prms(1).nb          = 1;  %1 or 2
prms(1).betawidth   = 0.02; % beta +/- beta width
prms(1).betaswitch  = 1/50; %0;

% computation
prms(1).max_evals   = 100; % for fzero to find beta
prms(1).max_iter    = 1e5; % for EGP
prms(1).tol_iter    = 1.0e-6;
prms(1).Nsim        = 100 %100000;
prms(1).Tsim        = 200;
prms(1).ExpandGridBetaIter  = 0; % Use larger grid for distributional computations (when iterating over betas)
prms(1).ExpandGridF = 1; % Use larger grid for final distributional computation
prms(1).nxlong      = 1000;  % larger grid size for ergodic distribution
prms(1).nxmpc       = 800; % larger grid size for MPC, must be divisible by 2
 
prms(1).targetAY    = 12.0;
prms(1).maxiterAY   = 20;
prms(1).tolAY       = 1.0e-4;

%mpc options
prms(1).mpcfrac{1}  = 1.0e-5; %approximate thoeretical mpc
prms(1).mpcfrac{2}  = 0.01; % 1 percent of average gross labor income: approx $500
prms(1).mpcfrac{3}  = 0.05; % 5 percent of average gross labor income: approx $5000

% OPTIONS
prms(1).IterateBeta         = 0;
prms(1).Display             = 0;
prms(1).MakePlots           = 1;
prms(1).ComputeDistMPC      = 1;
prms(1).ComputeSimMPC       = 1;
% prms(1).SolveDeterministic  = 0;
prms(1).Simulate            = 1;
prms(1).PrintStats          = 1;

%% Call model
Nprms = size(prms,2);
% Create structure array to store results
for ip = 1:Nprms
    results(ip) = egp_AR1_IID_tax_recode(prms(ip));
end
