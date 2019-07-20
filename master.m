clear;
close all;

% master script for discrete HA model
% need to set runopts in "SET OPTIONS" prior to running

%% ------------------------------------------------------------------------
% SET OPTIONS 
% -------------------------------------------------------------------------
% options
runopts.Display = 1;
runopts.Server = 1; % use server paths
runopts.IterateBeta = 1;
runopts.fast = 0; % very small asset and income grids for speed
runopts.Simulate = 0; % also solve distribution via simulation
runopts.mpcshocks_after_period1 = 0; % compute mpcs for ishock > 1

% directories
runopts.localdir = '/home/brian/Documents/GitHub/Discrete_HA';
runopts.serverdir = '/home/livingstonb/GitHub/Discrete_HA';

% grid tests, 0 to turn off
% 1-3 to select grid test parameters
% 2 is for simulations
runopts.GRIDTEST = 0; % 

% location of quarterly income file
QIncome = 'IncomeGrids/quarterly_b.mat';

% select only a subset of experiments (ignored when run on server)
% use empty cell array, {}, to run all
runopts.names_to_run = {'Q FixedBetaHet5 Width0.01 Death'};
% runopts.names_to_run = {'baseline_Q'};

%% ------------------------------------------------------------------------
% APPLY OPTIONS AND LOAD PARAMETERS
% -------------------------------------------------------------------------
if runopts.Server == 0
    runopts.path = runopts.localdir;
    runopts.number = [];
else
    runopts.number = str2num(getenv('SLURM_ARRAY_TASK_ID'));
    runopts.path = runopts.serverdir;
    runopts.savematpath = [runopts.serverdir '/Output/variables' num2str(runopts.number) '.mat'];
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
    params = parameters_grid_tests(runopts,QIncome);
elseif runopts.GRIDTEST == 2 % simulations
    params = parameters_grid_tests2(runopts,QIncome);
elseif runopts.GRIDTEST == 3
    params = parameters_grid_tests3(runopts,QIncome);
else
    params = parameters(runopts,QIncome);
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

if runopts.Server == 0
%     [decomp2,decomp3] = decomposition2(params,results);
%     % Create table
%     [T_annual,T_quarter] = create_table(params,results,...
%                                     decomps,decomp2,decomp3)
	disp('Check the results structure for detailed results')
else
    % convert MPCParams object to structure for saving
	Sparams = MPCParams.to_struct(params);
    save(runopts.savematpath,'Sparams','results','decomps')
    exit
end
