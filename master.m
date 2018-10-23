clear;
close all;

%% RUN OPTIONS
runopts.Batch = 1;
runopts.Server = 0;
runopts.fast = 1;

% empty string if not loading from file
% IncomeProcess = 'IncomeVariables/quarterly_a.mat';
IncomeProcess = '';

% select only a subset of experiments
selection.names_to_run = {'baseline_A','baseline_Q'}; % empty cell array to run all names
selection.suffix = '';
selection.frequencies = [1 4]; % [1 4], 1, or 4


%% FUNCTION CALL

[T_annual,T_quarter,results] = loop_through_main(runopts,IncomeProcess,selection);

if runopts.Server == 1
    exit
end