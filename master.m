clear;
close all;

%% RUN OPTIONS
runopts.Batch = 1; % use parameters.m, not parameters_experiment.m
runopts.Display = 1;
runopts.Server = 0; % use server paths and limit display
runopts.TryCatch = 1; % use try-catch block in main loop (auto-on if Server=1)
runopts.fast = 0; % specify very small asset and income grids for speed
runopts.localdir = '/Users/brianlivingston/Documents/GitHub/MPCrecode';

% empty string if not loading from file
% IncomeProcess = 'IncomeVariables/quarterly_a.mat';
IncomeProcess = '';

% select only a subset of experiments
selection.names_to_run = {'2 RandomBetaHet5 Width0.001 SwitchProb0.02 NoDeath'}; % empty cell array to run all names
selection.suffix = '';
selection.frequencies = [1]; % [1 4], 1, or 4

%% FUNCTION CALL

[T_annual,T_quarter,results,checks,exceptions] = loop_through_main(runopts,IncomeProcess,selection);

if runopts.Server == 1
    exit
end