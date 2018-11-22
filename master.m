clear;
close all;

%% RUN OPTIONS
runopts.Batch = 1; % use parameters.m, not parameters_experiment.m
runopts.Display = 1;
runopts.Server = 0; % use server paths and limit display
runopts.fast = 0; % specify very small asset and income grids for speed
runopts.Simulate = 1;
runopts.localdir = '/Users/Brian/Documents/GitHub/Discrete_HA';
runopts.GRIDTEST = 2;

% empty string if not loading from file
IncomeProcess = 'IncomeGrids/quarterly_b.mat';

% select only a subset of experiments
% ignored when run on server
selection.names_to_run = {}; % cell array of strings or {} to run all

%% ADD PATHS
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

%% LOAD PARAMETERIZATIONS
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
% convert to structure for saving
Sparams = MPCParams.to_struct(params);

%% CALL MAIN FUNCTION
Nparams = size(params,2);
checks     = cell(1,Nparams); % Information on failed sanity checks
decomps    = cell(1,Nparams); 

% iterate through specifications
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

%% DECOMPOSITIONS - COMPARISONS WITH BASELINE
decomp2 = decomposition2(params,results);

%% CREATE TABLE/SAVE VARIABLES
if runopts.Server == 0
    % Create table
    [T_annual,T_quarter] = create_table(params,results,...
                                    decomps,checks,decomp2);
    if ~isempty(T_annual)
        T_annual
    end

    if ~isempty(T_quarter)
        T_quarter
    end

    disp('Check the results structure for detailed results')
else
    % Save this run
    save(runopts.savematpath,'Sparams','results','decomps','checks','decomp2')
    exit
end
