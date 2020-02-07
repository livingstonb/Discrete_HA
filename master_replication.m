%% ONE-ASSET HOUSEHOLD MODEL
% This is the main script for this code repository
% HA model

clear;
close all;

%% ------------------------------------------------------------------------
% SET OPTIONS
% -------------------------------------------------------------------------
% Options
runopts.IterateBeta = 0; % irrelevant
runopts.fast = 0; % very small asset grids for testing
runopts.Simulate = 0; % also solve distribution via simulation
runopts.MPCs = 1;
runopts.MPCs_news = 1;
runopts.MPCs_loan_and_loss = 1;
runopts.DeterministicMPCs = 0; % irrelevant

% Directories
parent_dir = {'/home', 'brian', 'Documents', 'GitHub'};
runopts.localdir = fullfile(parent_dir{:}, 'Discrete_HA');

% Location of baseline income process
QIncome = fullfile('input', 'income_quarterly_b.mat');

% Shock sizes
% We assume that $500 is 0.81% of mean annual income
shocks = [-0.0081, -0.0405, -0.081, 0.0081, 0.0405, 0.081];

% Size of lump sum transfer (will be divided by 4)
lumptransfer = 0.0081 * 2.0 * 4.0;

%% ------------------------------------------------------------------------
% CHOOSE CALIBRATION
% -------------------------------------------------------------------------
% Which experiment to run
% 1 -- Beta heterogeneity
% 2 -- Baseline
% 3 -- Target fraction of hh with assets < $1000
runopts.number = 1;

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE BELOW
% -------------------------------------------------------------------------
runopts.path = runopts.localdir;
runopts.savematpath = fullfile(runopts.localdir, 'output', 'variables.mat');

if ~exist(runopts.localdir, 'dir')
    error('Directory not found')
end

outdir = fullfile(runopts.localdir, 'output');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

addpath(runopts.path);
addpath(fullfile(runopts.path, 'code'));
addpath(fullfile(runopts.path, 'code', 'aux_lib'));
cd(runopts.path);

%----------------------------------------------------------------------
% TARGET MEAN WEALTH = 3.2 * MEAN ANNUAL INCOME, W/DISC HETEROGENEITY
%----------------------------------------------------------------------
params = setup.Params(4, 'beta_heterogeneity', QIncome);
params.lumptransfer = lumptransfer;
params.shocks = shocks;
params.beta0 = 0.883375163;
params.betawidth = 0.0289;
params.nbeta = 3;

%----------------------------------------------------------------------
% TARGET MEAN WEALTH = 3.2 * MEAN ANNUAL INCOME
%----------------------------------------------------------------------
params(2) = setup.Params(4, 'baseline', QIncome);
params(2).lumptransfer = lumptransfer;
params(2).shocks = shocks;
params(2).beta0 = 0.984730508;

%----------------------------------------------------------------------
% TARGET FRACTION OF HOUSEHOLDS WITH LESS THAN $1000 IN ASSETS
%----------------------------------------------------------------------
params(3) = setup.Params(4, 'less_than_1000d', QIncome);
params(3).lumptransfer = lumptransfer;
params(3).shocks = shocks;
params(3).beta0 = 0.867871451;

%----------------------------------------------------------------------
% CALL METHODS/CHANGE SELECTED PARAMETERS, DO NOT CHANGE
%----------------------------------------------------------------------
params.set_index();
params = setup.Params.adjust_if_quarterly(params);
params.set_run_parameters(runopts);
params = setup.Params.select_by_number(params, runopts.number);

%% ------------------------------------------------------------------------
% SOLVE AND CREATE TABLE OF RESULTS
% -------------------------------------------------------------------------
results = main(params);

mpcs_on_table = true;
mpcs_news_on_table = true;
decomps_on_table = false;
table_gen = statistics.TableGenerator(...
	mpcs_on_table, mpcs_news_on_table, decomps_on_table);
results_table = table_gen.create(params, results, params.freq)