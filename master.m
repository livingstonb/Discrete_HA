clear;
close all;

%% RUN OPTIONS
runopts.Batch = 1; % use parameters.m, not parameters_experiment.m
runopts.Display = 1;
runopts.Server = 0; % use server paths and limit display
runopts.TryCatch = 0; % use try-catch block in main loop (auto-on if Server=1)
runopts.fast = 0; % specify very small asset and income grids for speed
runopts.localdir = '/Users/Brian/Documents/GitHub/MPCrecode';

% empty string if not loading from file
% IncomeProcess = 'IncomeVariables/quarterly_a.mat';
IncomeProcess = '';

% select only a subset of experiments
selection.names_to_run = {}; % empty cell array to run all names
selection.suffix = '';
selection.frequencies = [1 4]; % [1 4], 1, or 4

%% Add paths
    if runopts.Server == 0
        runopts.path = runopts.localdir;
    else
        runopts.path = '/home/livingstonb/GitHub/MPCrecode';
        runopts.savetablepath_annual = ['/home/livingstonb/output/table_annual' selection.suffix '.xls'];
        runopts.savetablepath_quarterly = ['/home/livingstonb/output/table_quarterly' selection.suffix '.xls'];
        runopts.savematpath = ['/home/livingstonb/output/variables' selection.suffix '.mat'];
    end
    addpath([runopts.path '/Auxiliary Functions']);
    addpath([runopts.path '/MPC Functions']);
    addpath([runopts.path '/Output Functions']);
    addpath([runopts.path '/EGP']);
    cd(runopts.path);

%% FUNCTION CALL

[T_annual,T_quarter,results,checks,exceptions] = loop_through_main(runopts,IncomeProcess,selection);

if runopts.Server == 1
    exit
end