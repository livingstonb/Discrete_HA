%% ONE-ASSET HOUSEHOLD MODEL
% This is the main script for this code repository
% HA model

% Prior to running this script:

% (1) Set options in the section below. Directory only needs to be set for
% either the server or the local directory, depending on whether
% runopts.Server = 0 or 1

% (2) In code/Model_Setup, modify the parameters script 'parameters.m' and make sure that 
% runopts.mode is equal to 'parameters'. Alternatively, create a new
% parameters script using parameters.m as a guide. Note that the current
% 'parameters.m' script assumes that the main income process is in
% input/income_quarterly_b.mat. Also note that if frequency is set to 4
% (quarterly), then annual parameter values should be used and they will
% be automatically adjusted in Params.adjust_if_quarterly()
%
% Note that all parameter defaults
% are set in the class file code/+setup/Params.m, and parameters.m overrides
% these defaults. Any parameters not in parameters.m are set to their
% defaults. See the properties of code/Model_Setup/Params.m for a list of all
% parameters.

% (3) Set runopts.names_to_run equal to a cell array containing the name
% of the parameterization to run, or use an empty cell array to loop
% over all parameterizations.

% (4) If convergence fails, betaH0 and/or betaL may need to be adjusted.
% betaL is the lower bound picked for beta during iteration and
% betaH0 is the adjustment factor to the upper bound. The code will
% guess a theoretical upper bound, and then will add betaH0 to
% to that value.

% RUNNING ON THE SERVER: To run in batch on the server, use 
% code/batch/server.sbatch as a template. That script sends an array to SLURM 
% that runs all of the requested parameters in parameters.m. Make sure
% that the range of numbers in the slurm array match the number of 
% parameters in the parameters script. Output files
% are stored in the Output directory

% OUTPUT: Results are stored in the 'results' structure. Its 'direct' property
% contains results found from computing the stationary distribution
% using non-simulation numerical methods. The 'sim' property contains results
% found from simulation, if the option is turned on.

clear;
close all;

%% ------------------------------------------------------------------------
% SET OPTIONS
% -------------------------------------------------------------------------
% options
runopts.Server = 0; % use server paths
runopts.IterateBeta = 0;
runopts.fast = 0; % very small asset and income grids for testing
runopts.Simulate = 0; % also solve distribution via simulation
runopts.MakePlots = 1;
runopts.MPCs = 1;
runopts.MPCs_news = 0;
runopts.MPCs_loan_and_loss = 0;
runopts.DeterministicMPCs = 1; % must be on if decompositions are needed

% directories
runopts.localdir = '/home/brian/Documents/GitHub/Discrete_HA';
runopts.serverdir = '/home/livingstonb/GitHub/Discrete_HA';

% name of parameters script
runopts.mode = 'parameters'; % 'parameters', 'grid_tests1', etc...

% select only a subset of experiments (ignored when run on server)
% use empty cell array, {}, to run all
runopts.names_to_run = {'quarterly_b_nyT101'};

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE BELOW
% -------------------------------------------------------------------------
if runopts.Server == 0
    runopts.path = runopts.localdir;
    runopts.number = [];
    runopts.outdir = fullfile(runopts.localdir, 'output');
    runopts.savematpath = [runopts.localdir '/output/variables' num2str(runopts.number) '.mat'];
    if ~exist(runopts.localdir, 'dir')
        error('Directory not found')
    end
    
    if ~exist([runopts.localdir '/output'], 'dir')
        mkdir([runopts.localdir '/output']);
    end
else
    runopts.number = str2num(getenv('SLURM_ARRAY_TASK_ID'));
    runopts.path = runopts.serverdir;
    runopts.outdir = fullfile(runopts.serverdir, 'output');
    runopts.savematpath = [runopts.serverdir '/output/variables' num2str(runopts.number) '.mat'];
    if ~exist(runopts.serverdir, 'dir')
        error('Directory not found')
    end
    
    if ~exist([runopts.localdir '/output'], 'dir')
        mkdir([runopts.localdir '/output']);
    end
end

if exist(runopts.savematpath, 'file') == 2
    % Delete old results
    delete runopts.savematpath;
end

addpath(runopts.path);
addpath([runopts.path '/code']);
addpath([runopts.path '/code/aux_lib']);
cd(runopts.path);

% Load parameters
switch runopts.mode
    case 'parameters'
        params = setup.parameters(runopts);
    case 'EZtests'
        params = setup.parameters_EZtests(runopts);
    case 'other'
        params = setup.parameters_other(runopts);
    otherwise
        error('Parameters script not found')
end
Nparams = size(params,2);

%% ------------------------------------------------------------------------
% CALIBRATIONS
% -------------------------------------------------------------------------
n_calibrations = 0;

% Vary discount rate to match mean wealth = 3.5
param_name = 'beta0';
stat_name = 'mean_a';
stat_target = 3.5;
n_calibrations = n_calibrations + 1;

% % Vary r to match P(a < $1000) = 0.23
% param_name = 'r';
% stat_name = 'wealth_lt_1000';
% stat_target = 0.23;
% n_calibrations = n_calibrations + 1;

% % Vary discount rate to match median wealth = 1.6
% param_name = 'beta0';
% stat_name = 'median_a';
% stat_target = 1.6;
% n_calibrations = n_calibrations + 1;

%% ------------------------------------------------------------------------
% CALIBRATING WITH FSOLVE
% -------------------------------------------------------------------------
if n_calibrations == 1
    param_init = params.(param_name);
    calibrator = solver.Calibrator(params, param_name,...
        stat_name, stat_target);
    beta_final = fsolve(@(x) calibrator.fn_handle(x, params), param_init);
elseif n_calibrations > 1
    error("Ensure that a max of one calibration is selected")
end

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION
% -------------------------------------------------------------------------
decomp_meanmpc = cell(1, Nparams); 

% iterate through specifications (or run 1)
for ip = 1:Nparams
    if params(ip).freq == 1
        msgfreq = 'annual';
    else
        msgfreq = 'quarterly';
    end
    fprintf('\n Trying %s parameterization "%s"\n', msgfreq,params(ip).name)

    tic
    [results(ip), decomp_meanmpc{ip}] = main(params(ip));
    toc
    disp(['Finished parameterization ' params(ip).name])
    
end

%% ------------------------------------------------------------------------
% DECOMPOSITION 2 AND SAVING/TABLE CREATING
% -------------------------------------------------------------------------
disp('Check the results structure for detailed results')
% convert Params object to structure for saving
Sparams = aux.to_structure(params);
% Sincome = aux.to_structure(income);
save(runopts.savematpath, 'Sparams', 'results', 'decomp_meanmpc')

if runopts.Server == 1
    exit
end


%% ------------------------------------------------------------------------
% SOLVE AND CREATE TABLE OF RESULTS
% -------------------------------------------------------------------------
return_nans = false;
decomp_with_loose_borr_limit = false;
[~, repagent_decomps] = statistics.baseline_repagent_decomps(params, results, return_nans);

table_gen = statistics.TableGenerator();
table_gen.decomp_repagent = repagent_decomps;
table_gen.decomp_incrisk = decomp_meanmpc;

quarterly_results = table_gen.create(params, results, 4);
annual_results = table_gen.create(params, results, 1);