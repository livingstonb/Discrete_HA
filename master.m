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
runopts.Server = true; % use server paths
runopts.calibrate = true;
runopts.fast = false; % very small asset and income grids for testing
runopts.Simulate = false; % also solve distribution via simulation
runopts.MakePlots = false;
runopts.MPCs = true;
runopts.MPCs_news = false;
runopts.MPCs_loan_and_loss = false;
runopts.DeterministicMPCs = true; % must be on if decompositions are needed

% directories
runopts.serverdir = '/home/livingstonb/GitHub/Discrete_HA';

if ismac
    runopts.localdir = '/Users/brianlivingston/Documents/GitHub/Discrete_HA';
elseif isunix
    runopts.localdir = '/home/brian/Documents/GitHub/Discrete_HA';
elseif ispc
    runopts.localdir = '/Users/brian-laptop/Documents/GitHub/Discrete_HA';
end

% name of parameters script
runopts.mode = 'parameters'; % 'parameters', 'grid_tests1', etc...

% select only a subset of experiments (ignored when run on server)
% use empty cell array, {}, to run all
runopts.names_to_run = {'Q EZ w/riskaver betw exp(-2), exp(2)'};
runopts.number = [];

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE BELOW
% -------------------------------------------------------------------------
if runopts.Server == 0
    runopts.path = runopts.localdir;
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
all_names = table();
switch runopts.mode
    case 'parameters'
        [params, all_names] = setup.parameters(runopts);
    case 'EZtests'
        params = setup.parameters_EZtests(runopts);
    case 'other'
        params = setup.parameters_other(runopts);
    case 'loose_borr_constraint'
        params = setup.parameters_loose_borr_constraint(runopts);
    otherwise
        error('Parameters script not found')
end
Nparams = size(params,2);

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION
% -------------------------------------------------------------------------
% iterate through specifications (or run 1)
if params.freq == 1
    msgfreq = 'annual';
else
    msgfreq = 'quarterly';
end
fprintf('\n Trying %s parameterization "%s"\n', msgfreq, params.name)

if params.calibrate
    options = optimoptions(@lsqnonlin, 'MaxIterations', params.calibrate_maxiter,...
            'FunctionTolerance', params.calibrate_tol);
    solver_args = params.calibrator.get_args();
    [calibrated_params, resnorm] = lsqnonlin(params.calibrator.solver_handle,...
        solver_args{:}, options);
    
    if resnorm > 1e-5
        error('Could not match targets')
    end
end

results = main(params);
disp(['Finished parameterization ' params.name])

if runopts.Server
    exit
end

%% ------------------------------------------------------------------------
% CREATE TABLE OF RESULTS
% -------------------------------------------------------------------------
return_nans = false;
decomp_with_loose_borr_limit = false;

table_gen = statistics.TableGenerator();

quarterly_results = table_gen.create(params, results, 4);
annual_results = table_gen.create(params, results, 1);

if ~isempty(quarterly_results)
    xlxpath = fullfile(runopts.outdir, 'quarterly_results.xlsx');
    writetable(quarterly_results, xlxpath, 'WriteRowNames', true);
end

if ~isempty(annual_results)
    xlxpath = fullfile(runopts.outdir, 'annual_results.xlsx');
    writetable(annual_results, xlxpath, 'WriteRowNames', true);
end