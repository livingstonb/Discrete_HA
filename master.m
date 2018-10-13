
%% HOUSEKEEPING
clear;
close all;


%% PARAMETERS IF NOT RUNNING IN BATCH

params0.name = 'params0';

% data frequency 
params0.freq        = 1; % 1 yearly, 4 quarterly

% returns
params0.r           = 0.02;

% demographics
params0.dieprob     = 1/50;

% preferences
params0.risk_aver   = 1;
params0.beta0       = 0.98; % annualized
params0.temptation  = 0;
params0.betaL       = 0.80;
% betaH defined in main function file

% bequests/annuities. bequest weight = 0 is accidental
params0.bequest_weight  = 0; %0.07;
params0.bequest_luxury  = 0.01; %0.01, must be >0 to avoid NaN error;
params0.Bequests = 1; % 1 for wealth left as bequest, 0 for disappears
params0.Annuities = 0; % Automatically turns off bequests if set to 1

% income risk: AR(1) + IID in logs
params0.LoadIncomeProcess = 0;
params0.nyT               = 11; %transitory component (not a state variable) (set to 1 for no Transitory Shocks)

% yT,yP (only relevant if LoadIncomeProcess==0)
params0.NormalizeY   = 1; % 1 to normalize gross income, 0 otherwise
params0.yTContinuous = 0; % doesn't seem to work properly
params0.sd_logyT     = sqrt(0.0497);  % 0.20, relevant if nyT>1
params0.lambdaT      = 1; % arrival rate of shocks;
params0.nyP          = 11; %11 persistent component
params0.sd_logyP     = sqrt(.0422); % 0.1950;
params0.rho_logyP    = 0.9525;
params0.nyF          = 1;
params0.sd_logyF     = 0;

% cash on hand / savings grid
params0.nx          = 100;
params0.xmax        = 1000;  % need high if using high-variance income shocks
params0.xgrid_par   = 0.3; %1 for linear, 0 for L-shaped
params0.borrow_lim  = 0;

%government
params0.labtaxlow       = 0; %proportional tax
params0.labtaxhigh      = 0; %additional tax on incomes above threshold
params0.labtaxthreshpc  = 0.99; %percentile of earnings distribution where high tax rate kicks in
params0.savtax          = 0; %0.0001;  %tax rate on savings
params0.savtaxthresh    = 0; %multiple of mean gross labor income

%discount factor shocks
params0.nb          = 5; % higher numbers dramatically increase computing load
params0.betawidth   = 0.015; % too large and eigs hangs while finding stat distribution
params0.betaswitch  = 1/50; %0;

% computation
params0.max_iter    = 1e5; % EGP
params0.tol_iter    = 1.0e-6; % EGP
params0.Nsim        = 100000; % 100000
params0.Tsim        = 200;
params0.nxlong      = 300; % Grid size for final computations

% beta iteration
params0.targetAY    = 3.5;
params0.maxiterAY   = 40;
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
params0.IterateBeta        = 0;
params0.Display            = 1;
params0.MakePlots          = 0;
params0.ComputeDirectMPC   = 1;
params0.Simulate           = 0;
Batch = 0; % Run alternate parameterizations

% Add paths
if Batch == 0
    path = '/Users/Brian/Documents/GitHub/MPCrecode';
    savetablepath = '/Users/Brian/Documents/table.xls';
    savematpath = '/Users/Brian/Documents/variables.mat';
else
    path = '/home/livingstonb/GitHub/MPCrecode';
    savetablepath = '/home/livingstonb/output/table.xls';
    savematpath = '/home/livingstonb/output/variables.mat';
end
addpath([path '/Auxiliary Functions']);
addpath([path '/MPC Functions']);
addpath([path '/Output Functions']);
cd(path);

%% LOAD ALTERNATE PARAMETERIZATIONS, STRUCTURE ARRAY
if Batch == 0
    params = params0;
else
    % Load parameters as defined in function
    params = parameters();
end

%% CALL MAIN FUNCTION
Nparams = size(params,2);

direct_results = cell(1,Nparams); % Results from direct computations
norisk_results = cell(1,Nparams); % Results from norisk model
sim_results    = cell(1,Nparams); % Results from simulations
exceptions     = cell(1,Nparams); % ME objects on any exceptions thrown
checks         = cell(1,Nparams); % Information on failed sanity checks
decomps        = cell(1,Nparams);

if Batch == 0
    [SR,DR,NR,checks{1},decomps{1}] = main(params(1));
    direct_results{1}  = DR;
    norisk_results{1}  = NR;
    sim_results{1}     = SR;      
else
    for ip = 1:Nparams
            disp(['Trying parameterization ' params(ip).name])
            try
                % Main function
                [SR,DR,NR,checks{ip},decomps{ip}] = main(params(ip));
                direct_results{ip}  = DR;
                norisk_results{ip}  = NR;
                sim_results{ip}     = SR;
                exception{ip} = 0; % main function completed
            catch ME
                checks{ip} = 'EXCEPTION_THROWN';
                exceptions{ip} = ME;
            end
            disp(['Finished parameterization ' params(ip).name])
    end
end

%% DECOMPOSITIONS - COMPARISONS WITH BASELINE
% Parameterization
for ip = 1:Nparams
    decomp2{ip} = struct([]);
    % Decomposition around a <= abar
    for ia = 1:numel(params(ip).abars)
        if params(ip).freq == 1
            baseind = 1;
        elseif params(ia).freq == 4
            baseind = 2;
        end

        try
            m1 = direct_results{ip}.mpcs1_a_direct{5};
            m0 = direct_results{baseind}.mpcs1_a_direct{5};
            g1 = direct_results{ip}.agrid_dist;
            g0 = direct_results{baseind}.agrid_dist;

            decomp2{ip}(ia).Em1_less_Em0 = direct_results{ip}.avg_mpc1_agrid(5) ...
                            - direct_results{baseind}.avg_mpc1_agrid(5);
            decomp2{ip}(ia).Em1_less_Em0_check = g1'*m1 - g0'*m0;
            decomp2{ip}(ia).term1 = g0' * (m1 - m0);
            decomp2{ip}(ia).term2 = m0' * (g1 - g0);
            decomp2{ip}(ia).term3 = (m1 - m0)' * (g1 - g0);
            ME = struct();
        catch ME
        end
        
        if Batch==0 || numel(fieldnames(ME))>0
            % Either not running as batch, or error was thrown
            decomp2{ip}(ia).Em1_less_Em0 = NaN;
            decomp2{ip}(ia).Em1_less_Em0_check = NaN;
            decomp2{ip}(ia).term1 = NaN;
            decomp2{ip}(ia).term2 = NaN;
            decomp2{ip}(ia).term3 = NaN;
        end
    end
end

%% SAVE AND CREATE TABLE
if Batch == 1
    save(savematpath,'sim_results','direct_results','norisk_results',...
                                                 'checks','exceptions');
end
                                             
T = create_table(params,direct_results,decomps,checks,exceptions,decomp2);
                            
if Batch == 1
    writetable(T,savetablepath,'WriteRowNames',true);
    exit
end
