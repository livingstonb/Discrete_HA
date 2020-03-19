
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

clear

options.server = false;
options.final_tables = false;
options.decomp_with_loose_borr_limit = false;
options.index_loose_borr_limit_Q = 'baseline_Q_with_borrowing';
options.index_loose_borr_limit_A = 'baseline_A_with_borrowing';

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
        stats{ind} = S.results.stats;
    end
end

for ip = 1:ind
    if params(ip).freq == 1
        baseind = find(ismember({params.name}, {'Annual'}));
    else
        baseind = find(ismember({params.name}, {'Quarterly'}));
    end

    if isempty(baseind)
        baseind = ip;
    end

    p0 = params(baseind);
    p1 = params(ip);
    stats0 = results(baseind).direct;
    stats1 = results(ip).direct;
    cdecomp = statistics.ComparisonDecomp(p0, p1, stats0, stats1);

    mpcs0 = reshape(stats0.mpcs(5).mpcs_1_t{1},...
        p0.nx_DST, []);
    mpcs1 = reshape(stats1.mpcs(5).mpcs_1_t{1},...
        p1.nx_DST, []);
    cdecomp.perform_decompositions(mpcs0, mpcs1);
    decomps_baseline(ip) = cdecomp.results;
end

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

table_gen = tables.StatsTable(params, stats);
table_gen.decomp_baseline = decomps_baseline;

if options.decomp_with_loose_borr_limit
    table_gen.decomp_incrisk_alt = decomp_alt;
end

table_out = table_gen.create(params, stats);

if ~isempty(table_out)
    xlxpath = fullfile(xlxdir, 'discrete_time_results.xlsx');
    writetable(table_out, xlxpath, 'WriteRowNames', true);
end

% quarterly_results = table_gen.create(params, results);
% annual_results = table_gen.create(params, results);

% if ~isempty(quarterly_results)
%     xlxpath = fullfile(xlxdir, 'quarterly_results.xlsx');
%     writetable(quarterly_results, xlxpath, 'WriteRowNames', true);
% end

% if ~isempty(annual_results)
%     xlxpath = fullfile(xlxdir, 'annual_results.xlsx');
%     writetable(annual_results, xlxpath, 'WriteRowNames', true);
% end

% save_tables = true;
% tables_out = tables.create_final_tables(params, results,...
%     decomps_baseline, save_tables);