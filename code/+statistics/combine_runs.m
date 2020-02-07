
% this script combines .mat files named variablesX.mat into
% an excel spreadsheet
%
% 'codedir' is the location of the 'code' directory
%
% 'matdir' is the location of the .mat files, if used
%
% 'FROM_MATFILE' indicates whether this script should use
% .mat files or should be run immediately after master.m
%
% 'xlxdir' is the desired directory of the output spreadsheet

% basedir = '/home/livingstonb/GitHub/Discrete_HA';
% matdir = '/home/livingstonb/GitHub/Discrete_HA/output/';
% xlxdir = '/home/livingstonb/GitHub/Discrete_HA/output/';
basedir = '/home/brian/Documents/GitHub/Discrete_HA';
matdir = '/home/brian/Documents/GitHub/Discrete_HA/output/';
xlxdir = '/home/brian/Documents/GitHub/Discrete_HA/output/';
FROM_MATFILE = false;

addpath([basedir '/code']);

if ~FROM_MATFILE
    clearvars -except params results decomp_meanmpc xlxdir FROM_MATFILE
else
    clearvars -except basedir matdir xlxdir FROM_MATFILE
end

if FROM_MATFILE
    %% Read .mat files into a cell array
    decomp_meanmpc = cell(1,1);
    decomp_baseline = cell(1,1);
    decomp_repagent = cell(1,1);
    
    ind = 0;
    for run = 1:999
        runstr = num2str(run);
        fpath = [matdir,'variables',runstr,'.mat'];
        if exist(fpath,'file')
            ind = ind+1;
            
            S = load(fpath);
            params(ind) = S.Sparams;
            results(ind) = S.results;
            decomps{ind} = S.decomp_meanmpc{1};
        else
            continue
        end
    end
else
    decomps = decomp_meanmpc;
end

[decomps_baseline,decomps_repagent] ...
    	= statistics.baseline_repagent_decomps(params,results);
    
[T_annual,T_quarter] = statistics.create_table_old(...
    params,results,decomps,decomps_baseline,decomps_repagent);

try
    writetable(T_quarter, [xlxdir 'T_quarter.xlsx'], 'WriteRowNames', true);
    writetable(T_annual, [xlxdir 'T_annual.xlsx'], 'WriteRowNames', true);
end