%% ONE-ASSET HOUSEHOLD MODEL
% This is the main script for this code repository
% HA model

% Prior to running this script:

% (1) Set options in the section below.

% (2) In code/Model_Setup, modify the parameters script 'parameters.m' and make sure that 
% runopts.mode is equal to 'parameters'. Alternatively, create a new
% parameters script using parameters.m as a guide. Note that the current
% 'parameters.m' script assumes that the main income process is in
% input/income_quarterly_b.mat. Also note that if frequency is set to 4
% (quarterly), then annual parameter values should be used and they will
% be automatically adjusted in Params.adjust_if_quarterly()

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
runopts.calibrate = true;
runopts.fast = false; % very small asset and income grids for testing
runopts.Simulate = false; % also solve distribution via simulation
runopts.MakePlots = false;
runopts.MPCs = true;
runopts.MPCs_news = false;
runopts.MPCs_loan_and_loss = false;
runopts.DeterministicMPCs = true; % must be on if decompositions are needed
runopts.SaveOutput = true;

% name of parameters script
runopts.mode = 'parameters'; % 'parameters', 'grid_tests1', etc...

% select only a subset of experiments (ignored when run on server)
runopts.names_to_run = {};
runopts.number = [];

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE BELOW
% -------------------------------------------------------------------------
[~, currdir] = fileparts(pwd());
if ~strcmp(currdir, 'Discrete_HA')
    msg = 'The user must cd into the Discrete_HA directory';
    bad_dir = MException('Discrete_HA:master', msg);
    throw(bad_dir);
end

server_array_id = str2num(getenv('SLURM_ARRAY_TASK_ID'));
running_on_server = ~isempty(server_array_id);

if running_on_server
    runopts.number = server_array_id;
end

matname = sprintf('variables%d.mat', runopts.number);
runopts.savematpath = fullfile('output', matname);

warning('off', 'MATLAB:MKDIR:DirectoryExists')
mkdir('output')
mkdir('temp')

if exist(runopts.savematpath, 'file') == 2
    % Delete old results
    delete runopts.savematpath;
end

if ~exist('+EconTools', 'dir')
    error("EconTools not found")
end

addpath('code');
addpath(fullfile('code', 'aux_lib'));

% Load parameters
[params, all_names] = setup.(runopts.mode)(runopts);

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION
% -------------------------------------------------------------------------
fprintf('\nParameterization "%s" was chosen.\n', params.name)

if params.calibrate
    disp('Beginning model calibration...')
    options = optimoptions(@lsqnonlin, 'MaxIterations', params.calibrate_maxiter,...
            'FunctionTolerance', params.calibrate_tol);
    solver_args = params.calibrator.get_args();
    [calibrated_params, resnorm] = lsqnonlin(params.calibrator.solver_handle,...
        solver_args{:}, options);

    if resnorm > 1e-5
        error('Could not match targets')
    end
end

results = main(params, 'iterating', false);
fprintf('Finished parameterization %s\n', params.name)

if running_on_server
    exit
end

%% ------------------------------------------------------------------------
% CREATE TABLE OF RESULTS
% -------------------------------------------------------------------------

% table_gen_detailed = tables.TableGenDetailed(params, results, params.freq);
% % table_gen_final = tables.TableFinal_Main(params, results, 1);
% % table_gen_decomps = tables.TableFinal_BaselineDecomps(params, results, 2);

% table_detailed = table_gen_detailed.create(params, results, params.freq);
% % table_final = table_gen_final.create(params, results);
% % table_decomps = table_gen_decomps.create(params, results);

% table_gen_detailed.save_table();
% % table_gen_final.save_table();

% % save_tables = true;
% % tables_out = tables.create_final_tables(params, results,...
% %     save_tables);

table_gen = tables.StatsTable(params, {results.stats});
table_out = table_gen.create(params, {results.stats})
