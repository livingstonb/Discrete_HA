
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
decomp_with_loose_borr_limit = true;
index_loose_borr_limit_Q = 'baseline_Q_with_borrowing';
index_loose_borr_limit_A = 'baseline_A_with_borrowing';

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

[decomps_baseline, decomps_repagent] ...
    	= statistics.baseline_repagent_decomps(params, results);

if decomp_with_loose_borr_limit
    for ip = 1:ind
        if params(ip).freq == 1
            index_loose_borr_limit = find(...
                cellfun(@(z) strcmp(z,index_loose_borr_limit_A), {params.name}));
        else
            index_loose_borr_limit = find(...
                cellfun(@(z) strcmp(z,index_loose_borr_limit_Q), {params.name}));
        end
        p_no_bc = params(index_loose_borr_limit);
        results_no_bc = results{index_loose_borr_limit};

        decomp_alt(ip) = alternate_decomposition(params(ip), results{ip},...
            p_no_bc, results_no_bc);
    end
end
    
% [T_annual, T_quarter] = statistics.create_table_old(...
%     params, results, decomps, decomps_baseline, decomps_repagent);

mpcs_on_table = true;
mpcs_news_on_table = true;
MPCs_loan_and_loss_on_table = true;
decomps_on_table = true;
table_gen = statistics.TableGenerator(...
    mpcs_on_table, mpcs_news_on_table, MPCs_loan_and_loss_on_table, decomps_on_table);

for freq = [1 4]
    params_freq = struct();
    results_freq = struct();
    decomp_meanmpc_freq = struct();
    repagent_decomps_freq = struct();
    ifreq = 0;
    for ip = 1:ind
        if params(ip).freq == freq
            ifreq = ifreq + 1;
            params_freq(ifreq) = params(ip);
            results_freq(ifreq) = results(ip);
            decomp_meanmpc_freq(ifreq) = decomp_meanmpc(ip);
            repagent_decomps_freq(ifreq) = repagent_decomps(ip);
        end
    end

    if (freq == 1) && (ifreq > 0)
        T_annual = table_gen.create(...
            params_freq, results_freq, freq, decomp_meanmpc_freq, repagent_decomps_freq);
    elseif (freq == 1)
        T_annual = table();
    elseif (freq == 4) && (ifreq > 0)
        T_quarter = table_gen.create(...
            params_freq, results_freq, freq, decomp_meanmpc_freq, repagent_decomps_freq);
    elseif (freq == 4)
        T_quarter = table();
    end
end

if ~isempty(T_quarter)
    writetable(T_quarter, [xlxdir 'T_quarter.xlsx'], 'WriteRowNames', true);
end

if ~isempty(T_annual)
    writetable(T_annual, [xlxdir 'T_annual.xlsx'], 'WriteRowNames', true);
end
