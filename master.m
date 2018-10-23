clear;
close all;

%% RUN OPTIONS
runopts.Batch = 1;
runopts.Server = 0;
runopts.fast = 0;

% empty string if not loading from file
% IncomeProcess = 'IncomeVariables/quarterly_a.mat';
IncomeProcess = 'IncomeVariables/quarterly_a.mat';

% select only a subset of experiments
selection.names_to_run = {'EZ ra8 invies1'}; % empty cell array to run all names
selection.suffix = '';
selection.frequencies = [1]; % [1 4], 1, or 4

%% FUNCTION CALL

[T_annual,T_quarter,results,checks,exceptions] = loop_through_main(runopts,IncomeProcess,selection);

if runopts.Server == 1
    exit
end