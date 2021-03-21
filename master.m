%% ONE ASSET DISCRETE TIME HA MODEL
% This is the main script for this code repository.
% See the readme for details.

clear;
close all;

%% ------------------------------------------------------------------------
% SET OPTIONS
% -------------------------------------------------------------------------
% options
runopts.calibrate = false; % wrap code in nonlinear solver
runopts.fast = false; % very small asset and income grids for testing
runopts.Simulate = false; % also solve distribution via simulation
runopts.MakePlots = false; % not used
runopts.MPCs = true;
runopts.MPCs_news = false;
runopts.MPCs_loan_and_loss = false;
runopts.DeterministicMPCs = true; % must be on if decompositions are needed
runopts.SaveOutput = true;

% name of parameters script
runopts.mode = 'parameters'; % 'parameters'

% select only a subset of experiments (ignored when run on server)
runopts.names_to_run = {};
runopts.number = [2];

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE
% -------------------------------------------------------------------------
% Check current working directory
[~, currdir] = fileparts(pwd());
if ~strcmp(currdir, 'Discrete_HA')
    msg = 'The user must cd into the Discrete_HA directory';
    bad_dir = MException('Discrete_HA:master', msg);
    throw(bad_dir);
end

% Get task id if running on server
server_array_id = str2num(getenv('SLURM_ARRAY_TASK_ID'));
running_on_server = ~isempty(server_array_id);
if running_on_server
    runopts.number = server_array_id;
end

% Path for .mat output file
matname = sprintf('variables%d.mat', runopts.number);
runopts.savematpath = fullfile('output', matname);

% Directories
warning('off', 'MATLAB:MKDIR:DirectoryExists')
mkdir('output')
mkdir('temp')

if exist(runopts.savematpath, 'file') == 2
    % Delete old results
    delete runopts.savematpath;
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
    % Calibrate with nonlinear solver
    disp('Beginning model calibration...')
    calibrator = params.calibrator;
    options = optimoptions(@lsqnonlin, 'MaxIterations', params.calibrate_maxiter,...
            'FunctionTolerance', params.calibrate_tol);
    solver_args = params.calibrator.get_args();
    calibrated_params = lsqnonlin(params.calibrator.solver_handle,...
        solver_args{:}, options);
end

% Final run
results = main(params, 'iterating', false);
fprintf('Finished parameterization %s\n', params.name)

%% ------------------------------------------------------------------------
% CREATE TABLE OF RESULTS
% -------------------------------------------------------------------------
table_out = tables.SingleTable(params, results.stats)