clear;
close all;

%% RUN OPTIONS
runopts.Batch = 1; % use parameters.m, not parameters_experiment.m
runopts.Display = 1;
runopts.Server = 1; % use server paths and limit display
runopts.fast = 0; % specify very small asset and income grids for speed
runopts.Simulate = 0;
runopts.localdir = '/Users/Brian/Documents/GitHub/Discrete_HA';

% empty string if not loading from file
IncomeProcess = 'IncomeGrids/quarterly_b.mat';

% select only a subset of experiments
% ignored when run on server
selection.names_to_run = {'baseline_Q'}; % cell array of strings or {} to run all

%% Add paths
if runopts.Server == 0
    runopts.path = runopts.localdir;
    selection.number = [];
else
    selection.number = str2num(getenv('SLURM_ARRAY_TASK_ID'));
    runopts.path = '/home/livingstonb/GitHub/Discrete_HA';
    runopts.savematpath = ['/home/livingstonb/GitHub/Discrete_HA/Output/variables' num2str(selection.number) '.mat'];
end
addpath([runopts.path '/Auxiliary Functions']);
addpath([runopts.path '/Solution Functions']);
addpath([runopts.path '/Output Functions']);
addpath([runopts.path '/Parameters']);
cd(runopts.path);

if exist(runopts.savematpath, 'file') == 2
    % Delete old results
    delete runopts.savematpath;
end

%% PARAMETERIZATIONS
if runopts.Batch == 0
    params = parameters_experiment(runopts,IncomeProcess); 
else
    params = parameters(runopts,selection,IncomeProcess);
end

Sparams = MPCParams.to_struct(params);

%% CALL MAIN FUNCTION
Nparams = size(params,2);
exceptions = cell(1,Nparams); % ME objects on any exceptions thrown
checks     = cell(1,Nparams); % Information on failed sanity checks
decomps    = cell(1,Nparams); 

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
    [T_annual,T_quarter] = create_table(params,results,...
                                    decomps,checks,exceptions,decomp2);
    if ~isempty(T_annual)
        T_annual
    end

    if ~isempty(T_quarter)
        T_quarter
    end

    disp('Check the results structure for detailed results')
else
    save(runopts.savematpath,'Sparams','results','decomps','checks','exceptions','decomp2')
    exit
end
