classdef FinalTables
    
    properties (Constant)
        table_includes = {
            {'Quarterly', 'Annual'}
            {'Quarterly', 'Annual'}
            {'Q1a'}
            {'Q1b'}
            {'Q2'}
            {'Q3'}
            {'Q4'}
            {'Q7'}
            {'Q8'}
        };
    end

    methods (Static)
        function save_table1(params_in, results, dirpath, varargin)
            header = tables.FinalTables.table1header(params_in, results, varargin{:});
            panelA = tables.FinalTables.table1panelA(params_in, results, varargin{:});
            panelB = tables.FinalTables.table1panelB(params_in, results, varargin{:});
            panelC = tables.FinalTables.table1panelC(params_in, results, varargin{:});
            panelD = tables.FinalTables.table1panelD(params_in, results, varargin{:});

            headerpath = fullfile(dirpath, 'table1_header.xlsx');
            writetable(header, headerpath, 'WriteRowNames', true);
            
            panelApath = fullfile(dirpath, 'table1_panelA.xlsx');
            writetable(panelA, panelApath, 'WriteRowNames', true);
            
            panelBpath = fullfile(dirpath, 'table1_panelB.xlsx');
            writetable(panelB, panelBpath, 'WriteRowNames', true);

            panelCpath = fullfile(dirpath, 'table1_panelC.xlsx');
            writetable(panelC, panelCpath, 'WriteRowNames', true);

            panelDpath = fullfile(dirpath, 'table1_panelD.xlsx');
            writetable(panelD, panelDpath, 'WriteRowNames', true);
        end

        function save_table2(params_in, results, dirpath, varargin)
            header = tables.FinalTables.table2header(params_in, results, varargin{:});
            panelsABC = tables.FinalTables.table2panelABC(params_in, results, varargin{:});
            panelD = tables.FinalTables.table2panelD(params_in, results, varargin{:});

            headerpath = fullfile(dirpath, 'table2_header.xlsx');
            writetable(header, headerpath, 'WriteRowNames', true);

            panelApath = fullfile(dirpath, 'table2_panelA.xlsx');
            writetable(panelsABC{1}, panelApath, 'WriteRowNames', true);
            
            panelBpath = fullfile(dirpath, 'table2_panelB.xlsx');
            writetable(panelsABC{2}, panelBpath, 'WriteRowNames', true);

            panelCpath = fullfile(dirpath, 'table2_panelC.xlsx');
            writetable(panelsABC{3}, panelCpath, 'WriteRowNames', true);

            panelDpath = fullfile(dirpath, 'table2_panelD.xlsx');
            writetable(panelD, panelDpath, 'WriteRowNames', true);
        end

        function save_experiment_table(params_in, results, comparison_decomps, dirpath, tableno)
            header = tables.FinalTables.experiment_table_header(params_in, results, tableno);
            panels = {tables.FinalTables.experiment_table_panelA(params_in, comparison_decomps, tableno)};
            panels{2} = tables.FinalTables.experiment_table_panelA2(params_in, comparison_decomps, tableno);
            panels{3} = tables.FinalTables.experiment_table_panelB(params_in, results, tableno);
            panels{4} = tables.FinalTables.experiment_table_panelC(params_in, results, tableno);
            panels{5} = tables.FinalTables.experiment_table_panelD(params_in, results, tableno);

            panel_labels = {'A', 'A2', 'B', 'C', 'D'};

            headerpath = fullfile(dirpath, sprintf('table%d_header.xlsx', tableno));
            writetable(header, headerpath, 'WriteRowNames', true);

            for ii = 1:numel(panel_labels)
                label = panel_labels{ii};
                panel_path = fullfile(dirpath, sprintf('table%d_panel%s.xlsx', tableno, label));
                writetable(panels{ii}, panel_path, 'WriteRowNames', true);
            end
        end

        function table_out = table1header(params_in, results, ctimeresults)
            params = filter_param_names(params_in, tables.FinalTables.table_includes{1});
            statistics = cell(numel(params), 1);

            get_stats = @(x) {      x.stats.mpcs(5).quarterly
                                    x.stats.mpcs(5).annual
                                    x.stats.beta_A
                                  };

            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = get_stats(results(ip));
            end

            if (nargin == 3)
                statistics{end+1} = get_stats(ctimeresults);
            end
            table_out = make_table(statistics, params, 'ctime_header', 'Continuous Time');
        end

        function table_out = table1panelA(params_in, results, ctimeresults)
            params = filter_param_names(params_in, tables.FinalTables.table_includes{1});
            statistics = cell(numel(params), 1);

            get_stats = @(x) {      x.stats.mean_gross_y_annual
                                    x.stats.std_log_gross_y_annual
                                    x.stats.std_log_net_y_annual
                                  };
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = get_stats(results(ip));
            end

            if (nargin == 3)
                statistics{end+1} = get_stats(ctimeresults);
            end
            table_out = make_table(statistics, params, 'ctime_header', 'Continuous Time');
        end
        
        function table_out = table1panelB(params_in, results, ctimeresults)
            params = filter_param_names(params_in, tables.FinalTables.table_includes{1});
            statistics = cell(numel(params), 1);

            get_stats = @(x) {  x.stats.mean_a
                                x.stats.sav0
                                x.stats.constrained{1}
                                x.stats.constrained_dollars{1}
                                x.stats.constrained_dollars{2}
                                x.stats.constrained_dollars{3}
                                x.stats.constrained_dollars{4}
                                x.stats.a_lt_ysixth
                                x.stats.a_lt_ytwelfth
                                x.stats.wpercentiles{1}
                                x.stats.wpercentiles{2}
                                x.stats.wpercentiles{3}
                                x.stats.wpercentiles{5}
                                x.stats.wpercentiles{7}
                                x.stats.wpercentiles{8}
                                x.stats.w_top10share
                                x.stats.w_top1share
                                x.stats.wgini
                              };
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = get_stats(results(ip));
            end

            if (nargin == 3)
                statistics{end+1} = get_stats(ctimeresults);
            end
            table_out = make_table(statistics, params, 'ctime_header', 'Continuous Time');
        end

        function table_out = table1panelC(params_in, results, ctimeresults)
            params = filter_param_names(params_in, tables.FinalTables.table_includes{1});
            statistics = cell(numel(params), 1);
            get_stats = @(x) {  x.stats.mpcs(4).annual
                                x.stats.mpcs(6).annual
                                x.stats.mpcs(4).quarterly
                                x.stats.mpcs(6).quarterly
                              };
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = get_stats(results(ip));
            end

            if (nargin == 3)
                statistics{end+1} = get_stats(ctimeresults);
            end
            table_out = make_table(statistics, params, 'ctime_header', 'Continuous Time');
        end

        function table_out = table1panelD(params_in, results, ctimeresults)
            params = filter_param_names(params_in, tables.FinalTables.table_includes{1});
            statistics = cell(numel(params), 1);
            get_stats = @(x) {  x.stats.mpcs(1).annual
                                x.stats.mpcs(2).annual
                                x.stats.mpcs(3).annual
                                x.stats.mpcs(1).quarterly
                                x.stats.mpcs(2).quarterly
                                x.stats.mpcs(3).quarterly
                              };

            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = get_stats(results(ip));
            end

            if (nargin == 3)
                statistics{end+1} = get_stats(ctimeresults);
            end
            table_out = make_table(statistics, params, 'ctime_header', 'Continuous Time');
        end

        function table_out = table2header(params_in, results, ctimeresults)
            params = filter_param_names(params_in, tables.FinalTables.table_includes{2});
            import statistics.Statistics.sfill

            statistics = cell(numel(params), 1);

            for ii = 1:numel(params)
                ip = params(ii).index;
                if (params(ii).freq == 1)
                    mean_mpc = sfill(results(ip).stats.mpcs(5).annual.value,...
                        'MPC, quarterly or annual (\%)', 1, 'MPC, quarterly or annual (\%)');
                else
                    mean_mpc = sfill(results(ip).stats.mpcs(5).quarterly.value,...
                        'MPC, quarterly or annual (\%)', 1, 'MPC, quarterly or annual (\%)');
                end
                statistics{ii} = {  mean_mpc
                                    results(ip).stats.decomp_norisk.term1_pct
                                  };
            end

            if (nargin == 3)
                mean_mpc = sfill(ctimeresults.stats.mpcs(5).quarterly.value,...
                        'MPC, quarterly or annual (\%)', 1, 'MPC, quarterly or annual (\%)');
                statistics{end+1} = {   mean_mpc
                                        ctimeresults.stats.decomp_norisk.term1_pct
                                    };
            end
            table_out = make_table(statistics, params, 'ctime_header', 'Continuous Time');
        end

        function table_out = table2panelABC(params_in, results, ctimeresults)
            params = filter_param_names(params_in, tables.FinalTables.table_includes{2});
            table_out = cell(1, 3);

            for kk = 1:3
                statistics = cell(numel(params), 1);

                get_stats = @(x) {  x.stats.decomp_norisk.term2(kk)
                                    x.stats.decomp_norisk.term3(kk)
                                    x.stats.decomp_norisk.term4(kk)
                                  };

                for ii = 1:numel(params)
                    ip = params(ii).index;
                    statistics{ii} = get_stats(results(ip));
                end

                if (nargin == 3)
                    statistics{end+1} = get_stats(ctimeresults);
                end
                table_out{kk} = make_table(statistics, params, 'ctime_header', 'Continuous Time');
            end
        end

        function table_out = table2panelD(params_in, results, ctimeresults)
            params = filter_param_names(params_in, tables.FinalTables.table_includes{2});
            statistics = cell(numel(params), 1);
            
            get_stats = @(x) {  x.stats.decomp_RA.Em1_less_mRA
                                x.stats.decomp_RA.term1
                                x.stats.decomp_RA.term2
                                x.stats.decomp_RA.term3
                              };
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = get_stats(results(ip));
            end

            if (nargin == 3)
                statistics{end+1} = get_stats(ctimeresults);
            end
            table_out = make_table(statistics, params, 'ctime_header', 'Continuous Time');
        end

        function table_out = experiment_table_header(params_in, results, tableno)
            indices = filter_param_group(params_in, tables.FinalTables.table_includes{tableno});

            import statistics.Statistics.sfill

            for ii = 1:numel(indices)
                ip = indices(ii);
                if (tableno == 8) || (tableno == 9)
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {
                            sfill(string(tex_vals.description), 'Description')
                        };
                    else
                        variable_values = {
                            sfill(nan, 'Description')
                        };
                    end
                elseif (tableno == 6) || (tableno == 7)
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {
                            sfill(string(tex_vals.riskaver), 'Risk aversion')
                            sfill(tex_vals.ies, 'IES')
                        };
                    else
                        variable_values = {
                            sfill(nan, 'Risk aversion', 2)
                            sfill(nan, 'IES', 3)
                        };
                    end
                elseif tableno == 5
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {
                            sfill(tex_vals.pswitch, 'Switch probability', 2)
                            sfill(tex_vals.width, 'Spacing', 3)
                        };
                    else
                        variable_values = {
                            sfill(nan, 'Switch probability', 2)
                            sfill(nan, 'Spacing', 3)
                        };
                    end
                elseif tableno == 4
                    variable_values = {};
                elseif tableno == 3
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {sfill(tex_vals.value, 'Value', 2)};
                    else
                        variable_values = {sfill(nan, 'Value', 2)};
                    end
                end
                statistics{ii} = {  results(ip).stats.mpcs(5).quarterly
                                    results(ip).stats.mpcs(5).annual
                                    results(ip).stats.beta_A
                                  };
                statistics{ii} = [variable_values(:); statistics{ii}];
                params(ii) = params_in(ip);
            end

            table_out = make_table(statistics, params, 'experiment', true);
        end

        function table_out = experiment_table_panelA(params_in, comparison_decomps, tableno)
            indices = filter_param_group(params_in, tables.FinalTables.table_includes{tableno});

            for ii = 1:numel(indices)
                ip = indices(ii);
                statistics{ii} = {  comparison_decomps(ip).Em1_less_Em0
                                    comparison_decomps(ip).term1
                                    comparison_decomps(ip).term2
                                    comparison_decomps(ip).term2a(2)
                                    comparison_decomps(ip).term2b(2)
                                    comparison_decomps(ip).term3
                                  };
                params(ii) = params_in(ip);
            end
            table_out = make_table(statistics, params, 'experiment', true);
        end

        function table_out = experiment_table_panelA2(params_in, comparison_decomps, tableno)
            indices = filter_param_group(params_in, tables.FinalTables.table_includes{tableno});

            for ii = 1:numel(indices)
                ip = indices(ii);
                statistics{ii} = {  comparison_decomps(ip).term1_pct
                                    comparison_decomps(ip).term2_pct
                                    comparison_decomps(ip).term2a_pct(2)
                                    comparison_decomps(ip).term2b_pct(2)
                                    comparison_decomps(ip).term3_pct
                                  };

                params(ii) = params_in(ip);
            end
            table_out = make_table(statistics, params, 'experiment', true);
        end

        function table_out = experiment_table_panelB(params_in, results, tableno)
            indices = filter_param_group(params_in, tables.FinalTables.table_includes{tableno});

            for ii = 1:numel(indices)
                ip = indices(ii);
                statistics{ii} = {  results(ip).stats.mean_a
                                    results(ip).stats.sav0
                                    results(ip).stats.constrained{1}
                                    results(ip).stats.constrained_dollars{1}
                                    results(ip).stats.constrained_dollars{2}
                                    results(ip).stats.constrained_dollars{3}
                                    results(ip).stats.constrained_dollars{4}
                                    results(ip).stats.a_lt_ysixth
                                    results(ip).stats.a_lt_ytwelfth
                                    results(ip).stats.w_top10share
                                    results(ip).stats.w_top1share
                                    results(ip).stats.wgini
                                  };
                params(ii) = params_in(ip);
            end
            table_out = make_table(statistics, params, 'experiment', true);
        end

        function table_out = experiment_table_panelC(params_in, results, tableno)
            indices = filter_param_group(params_in, tables.FinalTables.table_includes{tableno});

            for ii = 1:numel(indices)
                ip = indices(ii);
                statistics{ii} = {  results(ip).stats.mpcs(4).quarterly
                                    results(ip).stats.mpcs(6).quarterly
                                  };
                params(ii) = params_in(ip);
            end
            table_out = make_table(statistics, params, 'experiment', true);
        end

        function table_out = experiment_table_panelD(params_in, results, tableno)
            indices = filter_param_group(params_in, tables.FinalTables.table_includes{tableno});

            for ii = 1:numel(indices)
                ip = indices(ii);
                statistics{ii} = {  results(ip).stats.mpcs(1).quarterly
                                    results(ip).stats.mpcs(2).quarterly
                                    results(ip).stats.mpcs(3).quarterly
                                  };
                params(ii) = params_in(ip);
            end
            table_out = make_table(statistics, params, 'experiment', true);
        end
    end
    
end

function params = filter_param_names(params_in, includes)
    jj = 1;
    for ii = 1:numel(params_in)
        if ismember(params_in(ii).name, includes)
            if jj == 1
                params = params_in(ii);
            else
                params(jj) = params_in(ii);
            end
            jj = jj + 1;
        end
    end
end

function indices = filter_param_group(params_in, includes)
    indices = [];
    jj = 1;
    for ii = 1:numel(params_in)
        keep = false;
        for kk = 1:numel(params_in(ii).group)
            if ismember(params_in(ii).group{kk}, includes)
                keep = true;
                break;
            end
        end

        if keep
            indices(jj) = ii;
            jj = jj + 1;
        end
    end
end

function table_out = make_table(statistics, params, varargin)
    parser = inputParser;
    addOptional(parser, 'experiment', false);
    addOptional(parser, 'ctime_header', []);
    parse(parser, varargin{:});
    experiment = parser.Results.experiment;
    ctime_header = parser.Results.ctime_header;

    n = numel(params);
    
    for ii = 1:n
        vars{ii} = get_values(statistics{ii});

        if experiment
            varnames{ii} = char(params(ii).tex_header + string(sprintf('__v%d__', ii)));
        else
            varnames{ii} = params(ii).name;
        end
    end
    m = n + 1;

    if ~isempty(ctime_header)
        vars{m} = get_values(statistics{m});
        varnames{m} = ctime_header;
        m = m + 1;
    end

    vars{m} = get_precision(statistics{1});
    varnames{m} = 'decimals';

    rows = get_names(statistics{1});
    
    table_out = table(vars{:}, 'RowNames', rows(:), 'VariableNames', varnames(:));
end

function values = get_values(entries)
    values = {};
    for ii = 1:numel(entries)
        values{ii} = entries{ii}.value;
    end
    values = values(:);
end

function decimals = get_precision(entries)
    decimals = {};
    for ii = 1:numel(entries)
        decimals{ii} = entries{ii}.decimals;
    end
    decimals = decimals(:);
end

function names = get_names(entries)
    names = {};
    for ii = 1:numel(entries)
        names{ii} = entries{ii}.tex_label;
    end
    names = names(:);
end