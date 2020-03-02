
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

options.FROM_MATFILE = true;
options.server = true;
options.decomp_with_loose_borr_limit = false;
options.index_loose_borr_limit_Q = 'baseline_Q_with_borrowing';
options.index_loose_borr_limit_A = 'baseline_A_with_borrowing';


if options.server
    options.FROM_MATFILE = true;
end

if ~options.FROM_MATFILE
    clearvars -except params results options
else
    clearvars -except options
end

if ~options.server
    basedir = '/home/brian/Documents/GitHub/Discrete_HA';
    matdir = '/home/brian/Documents/GitHub/Discrete_HA/output/';
    xlxdir = '/home/brian/Documents/GitHub/Discrete_HA/output/';
else
    basedir = '/home/livingstonb/GitHub/Discrete_HA';
    matdir = '/home/livingstonb/GitHub/Discrete_HA/output/';
    xlxdir = '/home/livingstonb/GitHub/Discrete_HA/output/';
end

addpath([basedir '/code']);


mpcs_on_table = true;
mpcs_news_on_table = false;
MPCs_loan_and_loss_on_table = false;
decomps_on_table = true;


if options.FROM_MATFILE
    %% Read .mat files into a cell array
    ind = 0;
    for irun = 1:999
        runstr = num2str(irun);
        fpath = [matdir, 'variables', runstr, '.mat'];
        if exist(fpath,'file')
            ind = ind + 1;
            
            S = load(fpath);
            params(ind) = S.Sparams;
            results(ind) = S.results;
        end
    end
end

if options.decomp_with_loose_borr_limit
    return_nans = true;
else
    return_nans = false;
end
[decomps_baseline, ~] ...
    = statistics.baseline_repagent_decomps(params, results, return_nans);

if options.decomp_with_loose_borr_limit
    for ip = 1:ind
        if params(ip).freq == 1
            index_loose_borr_limit = find(...
                cellfun(@(z) strcmp(z,options.index_loose_borr_limit_A), {params.name}));
        else
            index_loose_borr_limit = find(...
                cellfun(@(z) strcmp(z,options.index_loose_borr_limit_Q), {params.name}));
        end
        p_no_bc = params(index_loose_borr_limit);
        results_no_bc = results(index_loose_borr_limit);

        return_nans = (ip == index_loose_borr_limit);
        decomp_alt(ip) = statistics.alternate_decomposition(...
            params(ip), results(ip),...
            p_no_bc, results_no_bc, return_nans);
    end
end

table_gen = statistics.TableGenerator();
table_gen.decomp_baseline = decomps_baseline;

if options.decomp_with_loose_borr_limit
    table_gen.decomp_incrisk_alt = decomp_alt;
end

quarterly_results = table_gen.create(params, results, 4);
annual_results = table_gen.create(params, results, 1);

if ~isempty(quarterly_results)
    xlxpath = fullfile(xlxdir, 'quarterly_results.xlsx');
    writetable(quarterly_results, xlxpath, 'WriteRowNames', true);
end

if ~isempty(annual_results)
    xlxpath = fullfile(xlxdir, 'annual_results.xlsx');
    writetable(annual_results, xlxpath, 'WriteRowNames', true);
end