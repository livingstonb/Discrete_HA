clear;
close all;

%% RUN OPTIONS
runopts.Batch = 0;
runopts.Server = 0;
runopts.fTast = 1;

% empty string if not loading from file
% IncomeProcess = 'IncomeVariables/quarterly_a.mat';
IncomeProcess = '';

% select only a subset of experiments
selection.names_to_run = {}; % empty cell array to run all names
selection.suffix = '';
selection.frequencies = [1 4]; % [1 4], 1, or 4

[T_annual,T_quarter] = loop_through_main(runopts,IncomeProcess,selection);

if runopts.Server == 1
    exit
end