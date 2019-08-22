%% ONE-ASSET HOUSEHOLD MODEL
% This is the main script for this code repository
% HA model

% Prior to running this script:

% (1) Set options in the section below. Directory only needs to be set for
% either the server or the local directory, depending on whether
% runopts.Server = 0 or 1

% (2) Modify the parameters script 'parameters.m' and make sure that 
% runopts.mode is equal to 'parameters'. Alternatively, create a new
% parameters script using parameters.m as a guide. Note that the current
% 'parameters.m' script assumes that the main income process is in
% IncomeGrids/quarterly_b.mat. Also note that if frequency is set to 4
% (quarterly), then annual parameter values should be used and they will
% be automatically adjusted in MPCParams.adjust_if_quarterly()
%
% Note that all parameter defaults
% are set in the class file Classes/MPCParams.m, and parameters.m overrides
% these defaults. Any parameters not in parameters.m are set to their
% defaults. See the properties of Classes/MPCParams.m for a list of all
% parameters.

% (3) Set runopts.names_to_run equal to a cell array containing the name
% of the parameterization to run, or use an empty cell array to loop
% over all parameterizations.

% (4) If convergence fails, betaH0 and/or betaL may need to be adjusted.

% RUNNING ON THE SERVER: To run in batch on the server, use 
% batch/server.sbatch as a template. That script sends an array to SLURM 
% that runs all of the requested parameters in parameters.m. Make sure
% that the range of numbers in the slurm array match the number of 
% parameters in the parameters script. Output files
% are stored in the Output directory

clear;
close all;

%% ------------------------------------------------------------------------
% SET OPTIONS
% -------------------------------------------------------------------------
% options
runopts.Server = 0; % use server paths
runopts.IterateBeta = 1;
runopts.fast = 0; % very small asset and income grids for speed
runopts.Simulate = 0; % also solve distribution via simulation
runopts.mpcshocks_after_period1 = 0; % compute mpcs for ishock > 1

% directories
runopts.localdir = '/Users/Brian-laptop/Documents/GitHub/Discrete_HA';
runopts.serverdir = '/home/livingstonb/GitHub/Discrete_HA';

% name of parameters script
runopts.mode = 'parameters'; % 'parameters', 'grid_tests1', etc...

% select only a subset of experiments (ignored when run on server)
% use empty cell array, {}, to run all
runopts.names_to_run = {'Q EZ ra1 ies2'}; % {'baseline_Q'}

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE BELOW
% -------------------------------------------------------------------------
if runopts.Server == 0
    runopts.path = runopts.localdir;
    runopts.number = [];
    runopts.savematpath = [runopts.localdir '/Output/variables' num2str(runopts.number) '.mat'];
else
    runopts.number = str2num(getenv('SLURM_ARRAY_TASK_ID'));
    runopts.path = runopts.serverdir;
    runopts.savematpath = [runopts.serverdir '/Output/variables' num2str(runopts.number) '.mat'];
end

if exist(runopts.savematpath, 'file') == 2
    % Delete old results
    delete runopts.savematpath;
end

addpath([runopts.path '/Auxiliary Functions']);
addpath([runopts.path '/Solution Functions']);
addpath([runopts.path '/Output Functions']);
addpath([runopts.path '/Parameters']);
addpath([runopts.path '/Classes']);
cd(runopts.path);

% Load parameters
switch runopts.mode
    case 'parameters'
        params = parameters(runopts);
    case 'grid_tests1'
        params = parameters_grid_tests1(runopts,'IncomeGrids/quarterly_b.mat');
    case 'grid_tests2'
        params = parameters_grid_tests2(runopts,'IncomeGrids/quarterly_b.mat');
    case 'grid_tests3'
        params = parameters_grid_tests3(runopts,'IncomeGrids/quarterly_b.mat');
    otherwise
        error('Parameters script not found')
end

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION
% -------------------------------------------------------------------------
Nparams = size(params,2);
decomps    = cell(1,Nparams); 

% iterate through specifications (or run 1)
for ip = 1:Nparams
    if params(ip).freq == 1
        msgfreq = 'annual';
    else
        msgfreq = 'quarterly';
    end
    fprintf('\n Trying %s parameterization "%s"\n',msgfreq,params(ip).name)

    tic
    [results(ip),decomps{ip}] = main(params(ip));
    toc
    disp(['Finished parameterization ' params(ip).name])
    
end

%% ------------------------------------------------------------------------
% DECOMPOSITION 2 AND SAVING/TABLE CREATING
% -------------------------------------------------------------------------

disp('Check the results structure for detailed results')
% convert MPCParams object to structure for saving
Sparams = MPCParams.to_struct(params);
save(runopts.savematpath,'Sparams','results','decomps')

if runopts.Server == 1
    exit
end
