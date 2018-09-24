
%% Housekeeping

clear;
close all;

path = '/Users/Brian/Documents/GitHub/MPCrecode';
addpath([path '/Auxiliary Functions']);
cd(path);

%% Set parameters
prms = struct();

% returns
prms(1).r           = 0.02;
prms(1).R           = 1 + prms(1).r;

% demographics
prms(1).dieprob     = 0;

% preferences
prms(1).risk_aver   = 1;
prms(1).beta0       = 0.97365;
prms(1).temptation  = 0;
prms(1).betaL       = 0.90;
prms(1).betaH       = 1/(prms(1).R*(1-prms.dieprob));

%warm glow bequests: bequessgrt_weight = 0 is accidental
prms(1).bequest_weight = 0; %0.07;
prms(1).bequest_luxury = 0.01; %0.01, must be >0 to avoid NaN error;

% income risk: AR(1) + IID in logs
prms(1).LoadIncomeProcess   = 1;
prms(1).nyT                 = 21; %transitory component (not a state variable) (set to 1 for no Transitory Shocks)

%only relevant if LoadIncomeProcess==0
prms(1).sd_logyT    = sqrt(0.2);  %0.20; %relevant if nyT>1
prms(1).lambdaT     = 1; % arrival rate of shocks;
prms(1).nyP         = 11; %11 persistent component
prms(1).sd_logyP    = sqrt(0.1950); %0.1950;
prms(1).rho_logyP   = 0.9947;
prms(1).nyF         = 1;
prms(1).sd_logyF    = 0;

% cash on hand / savings grid
prms(1).nx          = 1000;
prms(1).xmax        = 40;  %multiple of mean gross labor income
prms(1).xgrid_par   = 0.2; %1 for linear, 0 for L-shaped
prms(1).borrow_lim  = 0;

%government
prms(1).labtaxlow       = 0; %proportional tax
prms(1).labtaxhigh      = 0; %additional tax on incomes above threshold
prms(1).labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in

prms(1).savtax          = 0; %0.0001;  %tax rate on savings
prms(1).savtaxthresh    = 0; %multiple of mean gross labor income


%discount factor shocks;
prms(1).nb          = 1;  %1 or 2
prms(1).betawidth   = 0.065; % beta +/- beta width
prms(1).betaswitch  = 1/50; %0;

% computation
prms(1).NormalizeY  = 1;
prms(1).max_iter    = 1e5;
prms(1).tol_iter    = 1.0e-6;
prms(1).Nsim        = 100000;
prms(1).Tsim        = 200;
prms(1).ergodic_tol  = 1e-8;
 
prms(1).targetAY    = 12.0;
prms(1).maxiterAY   = 20;
prms(1).tolAY       = 1.0e-4;

%mpc options
prms(1).mpcfrac{1}  = 1.0e-10; %approximate thoeretical mpc
prms(1).mpcfrac{2}  = 0.01; % 1 percent of average gross labor income: approx $500
prms(1).mpcfrac{3}  = 0.05; % 5 percent of average gross labor income: approx $5000

% OPTIONS
prms(1).IterateBeta         = 0;
prms(1).Display             = 1;
prms(1).MakePlots           = 1;
prms(1).ComputeMPC          = 1;
prms(1).SolveDeterministic  = 0;
prms(1).Simulate            = 1;
prms(1).PrintStats          = 1;

%% Call model
Nprms = size(prms,2);
% Create structure array to store results
for ip = 1:Nprms
    results(ip) = egp_AR1_IID_tax_recode(prms(ip));
end
