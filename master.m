clear;
close all;

%% ------------------------------------------------------------------------
% SET OPTIONS 
% -------------------------------------------------------------------------
runopts.Batch = 0; % use parameters.m, not parameters_experiment.m
runopts.Display = 1;
runopts.Server = 0; % use server paths
runopts.fast = 0; % very small asset and income grids for speed
runopts.Simulate = 0;
runopts.localdir = '/Users/Brian/Documents/GitHub/Discrete_HA';

% local grid tests, 0 to turn off, 1 for transition probs, 2 for simulations
runopts.GRIDTEST = 0; % 

IncomeProcess = 'IncomeGrids/quarterly_b.mat';

% select only a subset of experiments (ignored when run on server)
% use empty cell array, {}, to run all
selection.names_to_run = {};

%% ------------------------------------------------------------------------
% APPLY OPTIONS AND LOAD PARAMETERS
% -------------------------------------------------------------------------
if runopts.Server == 0
    runopts.path = runopts.localdir;
    selection.number = [];
else
    selection.number = str2num(getenv('SLURM_ARRAY_TASK_ID'));
    runopts.path = '/home/livingstonb/GitHub/Discrete_HA';
    runopts.savematpath = ['/home/livingstonb/GitHub/Discrete_HA/Output/variables' num2str(selection.number) '.mat'];
    if exist(runopts.savematpath, 'file') == 2
        % Delete old results
        delete runopts.savematpath;
    end
end
addpath([runopts.path '/Auxiliary Functions']);
addpath([runopts.path '/Solution Functions']);
addpath([runopts.path '/Output Functions']);
addpath([runopts.path '/Parameters']);
cd(runopts.path);

% Load parameters
if runopts.GRIDTEST == 1
    % only setup to run locally
    params = parameters_grid_tests(runopts,IncomeProcess);
elseif runopts.GRIDTEST == 2
    params = parameters_grid_tests2(runopts,selection,IncomeProcess);
elseif runopts.Batch == 0
    params = parameters_experiment(runopts,IncomeProcess); 
else
    params = parameters(runopts,selection,IncomeProcess);
end

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION
% -------------------------------------------------------------------------
Nparams = size(params,2);
checks     = cell(1,Nparams); % Information on failed sanity checks
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
    [results(ip),checks{ip},decomps{ip}] = main(params(ip));
    toc
    disp(['Finished parameterization ' params(ip).name])
    
end

%% ------------------------------------------------------------------------
% DECOMPOSITION 2 AND SAVING/TABLE CREATING
% -------------------------------------------------------------------------
decomp2 = decomposition2(params,results);

if runopts.Server == 0
    % Create table
    [T_annual,T_quarter] = create_table(params,results,...
                                    decomps,checks,decomp2)
    disp('Check the results structure for detailed results')
else
    % convert MPCParams object to structure for saving
	Sparams = MPCParams.to_struct(params);
    save(runopts.savematpath,'Sparams','results','decomps','checks','decomp2')
    exit
end
