
clear;
close all;

%% RUN OPTIONS
Batch = 1;
Server = 1;
SmallGrid = 0;

% empty string if not loading from file
IncomeProcess = '/home/livingstonb/GitHub/MPCrecode/IncomeVariables/quarterly_a.mat'; 

% select only a subset of experiments
names_to_run = {}; % empty cell array to run all names
suffix = '_Qa';
Frequencies = [4]; % [1 4], 1, or 4

%% PARAMETERS IF NOT RUNNING IN BATCH

% Do not change
params0.name = 'params0';
params0.index = 1;
params0.IncomeProcess = IncomeProcess;

% data frequency 
params0.freq        = 1; % 1 yearly, 4 quarterly

% returns
params0.r           = 0.02;

% demographics
params0.dieprob     = 1/50;

% preferences
params0.EpsteinZin  = 0; % no bequest utility or temptation allowed
params0.invies      = 2.5; % only relevant if doing Epstein-Zin
params0.risk_aver   = 1;
params0.beta0       = 0.95; % annualized
params0.temptation  = 0;
params0.betaL       = 0.80;
% betaH defined in main function file

% bequests/annuities. bequest weight = 0 is accidental
params0.bequest_weight  = 0; %0.07;
params0.bequest_curv        = 1;
params0.bequest_luxury  = 1; %0.01, must be >0 to avoid NaN error;
params0.Bequests = 1; % 1 for wealth left as bequest, 0 for disappears
params0.Annuities = 0; % Automatically turns off bequests if set to 1

% income risk: AR(1) + IID in logs
params0.nyT               = 11; %transitory component (not a state variable) (set to 1 for no Transitory Shocks)

% yT,yP (only relevant if LoadIncomeProcess==0)
params0.yTContinuous = 0;
params0.sd_logyT     = sqrt(0.0497);  % 0.20, relevant if nyT>1
params0.lambdaT      = 1; % arrival rate of shocks;
params0.nyP          = 11; %11 persistent component
params0.sd_logyP     = sqrt(.0422); % 0.1950;
params0.rho_logyP    = 0.9525;
params0.nyF          = 1;
params0.sd_logyF     = sqrt(0.1801);

% cash on hand / savings grid
params0.nx          = 150;
params0.xmax        = 1000;  % need high if using high-variance income shocks
params0.xgrid_par   = 1/3; %1 for linear, 0 for L-shaped
params0.borrow_lim  = 0;

%government
params0.labtaxlow       = 0; %proportional tax
params0.labtaxhigh      = 0; %additional tax on incomes above threshold
params0.labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in
params0.savtax          = 0; %0.0001;  %tax rate on savings
params0.savtaxthresh    = 0; %multiple of mean gross labor income

%discount factor shocks
params0.nb          = 1; % higher numbers dramatically increase computing load
params0.betawidth   = 0.001;
params0.betaswitch  = 0.1; %0;

% computation
params0.max_iter    = 2e5; % EGP
params0.tol_iter    = 1.0e-6; % EGP
params0.Nsim        = 100000; % 100000
params0.Tsim        = 100;
params0.nxlong      = 500; % Grid size for final computations

% beta iteration
params0.targetAY    = 3.5;
params0.maxiterAY   = 50;
params0.tolAY       = 1e-5;

% mpc options
params0.Nmpcsim     = 1e6;
params0.mpcfrac(1)  = -1e-5; %approximate thoeretical mpc
params0.mpcfrac(2)  = -0.01;
params0.mpcfrac(3)  = -0.1;
params0.mpcfrac(4)  = 1e-5; % approximate thoeretical mpc
params0.mpcfrac(5)  = 0.01; % used in decomposition, 1 percent of average gross labor income: approx $500
params0.mpcfrac(6)  = 0.1; % 5 percent of average gross labor income: approx $5000

% wealth statistics options
params0.epsilon = [0 0.005 0.01 0.02 0.05 0.1]; % fraction of mean labor income
params0.percentiles = [10 25 50 75 90 95 99 99.9]; % in percent

% decomposition
params0.abars = [0 0.01 0.05];

% OPTIONS
params0.IterateBeta        = 1;
params0.Display            = 1;
params0.MakePlots          = 0;
params0.Simulate           = 0;

% Add paths
if Server == 0
    path = '/Users/Brian/Documents/GitHub/MPCrecode';
    savetablepath_annual = '/Users/Brian/Documents/table_annual.xls';
    savetablepath_quarterly = '/Users/Brian/Documents/table_quarterly.xls';
    savematpath = '/Users/Brian/Documents/variables.mat';
else
    path = '/home/livingstonb/GitHub/MPCrecode';
    savetablepath_annual = ['/home/livingstonb/output/table_annual' suffix '.xls'];
    savetablepath_quarterly = ['/home/livingstonb/output/table_quarterly' suffix '.xls'];
    savematpath = ['/home/livingstonb/output/variables' suffix '.mat'];
end
addpath([path '/Auxiliary Functions']);
addpath([path '/MPC Functions']);
addpath([path '/Output Functions']);
addpath([path '/EGP']);
cd(path);

%% LOAD ALTERNATE PARAMETERIZATIONS
if Batch == 0
    params = params0;
else
    params = parameters(SmallGrid,IncomeProcess);
end

%% SELECT SUBSET OF SPECIFICATIONS (ONLY RELEVANT FOR BATCH)

if Batch == 1
    if isempty(names_to_run)
        params_to_run = 1:numel(params);
    else
        % Indices of selected names within params
        params_to_run = find(ismember({params.name},names_to_run));
    end
else
    params_to_run = 1;
end

params = params(params_to_run);

if all(ismember(Frequencies,[1 4]))
    params = params(ismember([params.freq],Frequencies));
else
    error('Frequencies must be [1 4], 1, or 4')
end

Nparams = numel(params);

%% CALL MAIN FUNCTION

direct_results = cell(1,Nparams); % Results from direct computations
norisk_results = cell(1,Nparams); % Results from norisk model
sim_results    = cell(1,Nparams); % Results from simulations
exceptions     = cell(1,Nparams); % ME objects on any exceptions thrown
checks         = cell(1,Nparams); % Information on failed sanity checks
decomps        = cell(1,Nparams); 
income         = cell(1,Nparams);

if Batch == 0
    tic
    [income{1},SR,DR,NR,checks{1},decomps{1}] = main(params(1));
    toc
    direct_results{1}  = DR;
    norisk_results{1}  = NR;
    sim_results{1}     = SR;      
else
    for ip = 1:Nparams
        tic
        disp(['Trying parameterization ' params(ip).name])
        try
            % Main function
            [income{ip},SR,DR,NR,checks{ip},decomps{ip}] = main(params(ip));
            direct_results{ip}  = DR;
            norisk_results{ip}  = NR;
            sim_results{ip}     = SR;
            exceptions{ip} = []; % main function completed
        catch ME
            checks{ip} = 'EXCEPTION_THROWN';
            exceptions{ip} = ME;
        end
        disp(['Finished parameterization ' params(ip).name])
            toc
            
        if Server == 1
        save(savematpath,'sim_results','direct_results','norisk_results',...
                            'checks','exceptions','params','decomps');
        end
    end
end

%% DECOMPOSITIONS - COMPARISONS WITH BASELINE
decomp2 = decomposition2(params,direct_results,exceptions);

%% SAVE VARIABLES AND CREATE TABLE
                                             
[T_annual,T_quarter] = create_table(params,direct_results,decomps,checks,exceptions,decomp2);
                            
if Server == 1
    if ~isempty(T_annual)
        writetable(T_annual,savetablepath_annual,'WriteRowNames',true);
    end
    if ~isempty(T_quarter)
        writetable(T_quarter,savetablepath_quarterly,'WriteRowNames',true);
    end
    exit
end
