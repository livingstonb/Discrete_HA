
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

FROM_MATFILE = true;

if ~FROM_MATFILE
    clearvars -except params results decomp_meanmpc xlxdir FROM_MATFILE
else
    clearvars -except FROM_MATFILE
end

basedir = '/home/brian/Documents/GitHub/Discrete_HA';
matdir = '/home/brian/Documents/GitHub/Discrete_HA/output/';
xlxdir = '/home/brian/Documents/GitHub/Discrete_HA/output/';

addpath([basedir '/code']);

decomp_with_loose_borr_limit = true;
index_loose_borr_limit_Q = 'baseline_Q_with_borrowing';
index_loose_borr_limit_A = 'baseline_A_with_borrowing';

mpcs_on_table = true;
mpcs_news_on_table = false;
MPCs_loan_and_loss_on_table = false;
decomps_on_table = true;


if FROM_MATFILE
    %% Read .mat files into a cell array
    decomp_meanmpc = cell(1,1);
    decomp_baseline = cell(1,1);
    decomp_repagent = cell(1,1);
    
    ind = 0;
    for run = 1:999
        runstr = num2str(run);
        fpath = [matdir,'variables', runstr, '.mat'];
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

if decomp_with_loose_borr_limit
    return_nans = true;
else
    return_nans = false;
end
[decomps_baseline, decomps_repagent] ...
    	= statistics.baseline_repagent_decomps(params, results, return_nans);

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
        results_no_bc = results(index_loose_borr_limit);

        return_nans = (ip == index_loose_borr_limit);
        decomp_alt(ip) = statistics.alternate_decomposition(...
            params(ip), results(ip),...
            p_no_bc, results_no_bc, return_nans);
    end
end
    
% [T_annual, T_quarter] = statistics.create_table_old(...
%     params, results, decomps, decomps_baseline, decomps_repagent);

table_gen = statistics.TableGenerator(...
    mpcs_on_table, mpcs_news_on_table,...
    MPCs_loan_and_loss_on_table,...
    decomps_on_table,...
    decomp_with_loose_borr_limit);

for freq = [1 4]
    ifreq = 0;
    for ip = 1:ind
        if params(ip).freq == freq
            ifreq = ifreq + 1;

            if ifreq == 1
                params_freq = params(ip);
                results_freq = results(ip);
                decomp_meanmpc_freq{1} = decomps{ip};
                decomps_repagent_freq = decomps_repagent(ip);
                decomp_alt_freq = decomp_alt(ip);
            else
                params_freq(ifreq) = params(ip);
                results_freq(ifreq) = results(ip);
                decomp_meanmpc_freq{ifreq} = decomps{ip};
                decomps_repagent_freq(ifreq) = decomps_repagent(ip);
                decomp_alt_freq(ifreq) = decomp_alt(ip);
            end
        end
    end

    if ifreq > 0
        T = table_gen.create(...
            params_freq, results_freq, freq,...
            decomp_meanmpc_freq,...
            decomps_repagent_freq,...
            decomp_alt_freq);

        if freq == 1
            T_annual = T;
            fname = 'T_annual.xlsx';
        elseif freq == 4
            T_quarter = T;
            fname = 'T_quarter.xlsx';
        end
        writetable(T, [xlxdir fname], 'WriteRowNames', true);
    end
end