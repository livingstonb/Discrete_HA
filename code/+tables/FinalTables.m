classdef FinalTables
    
    properties (Constant)
        table1includes = {'Quarterly', 'Annual'};
        table2includes = {'Quarterly', 'Annual'};
        table3includes = {'Q1a'};
    end

    methods (Static)
        function save_table1(params_in, results, dirpath)
            header = tables.FinalTables.table1header(params_in, results);
            panelA = tables.FinalTables.table1panelA(params_in, results);
            panelB = tables.FinalTables.table1panelB(params_in, results);
            panelC = tables.FinalTables.table1panelC(params_in, results);
            panelD = tables.FinalTables.table1panelD(params_in, results);

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

        function save_table2(params_in, results, dirpath)
            header = tables.FinalTables.table2header(params_in, results);
            panelsABC = tables.FinalTables.table2panelABC(params_in, results);
            panelD = tables.FinalTables.table2panelD(params_in, results);

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

        function table_out = table1header(params_in, results)
            params = filter_param_names(params_in, tables.FinalTables.table1includes);
            statistics = cell(numel(params), 1);
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  results(ip).stats.mpcs(5).quarterly
                                    results(ip).stats.mpcs(5).annual
                                    results(ip).stats.beta_A
                                  };
            end
            table_out = make_table(statistics, params);
        end

        function table_out = table1panelA(params_in, results)
            params = filter_param_names(params_in, tables.FinalTables.table1includes);
            statistics = cell(numel(params), 1);
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  results(ip).stats.mean_gross_y_annual
                                    results(ip).stats.std_log_gross_y_annual
                                    results(ip).stats.std_log_net_y_annual
                                  };
            end
            table_out = make_table(statistics, params);
        end
        
        function table_out = table1panelB(params_in, results)
            params = filter_param_names(params_in, tables.FinalTables.table1includes);
            statistics = cell(numel(params), 1);
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  results(ip).stats.mean_a
                                    results(ip).stats.sav0
                                    results(ip).stats.constrained{1}
                                    results(ip).stats.constrained_dollars{1}
                                    results(ip).stats.constrained_dollars{2}
                                    results(ip).stats.constrained_dollars{3}
                                    results(ip).stats.constrained_dollars{4}
                                    results(ip).stats.a_lt_ysixth
                                    results(ip).stats.a_lt_ytwelfth
                                    results(ip).stats.wpercentiles{1}
                                    results(ip).stats.wpercentiles{2}
                                    results(ip).stats.wpercentiles{3}
                                    results(ip).stats.wpercentiles{5}
                                    results(ip).stats.wpercentiles{7}
                                    results(ip).stats.wpercentiles{8}
                                    results(ip).stats.w_top10share
                                    results(ip).stats.w_top1share
                                    results(ip).stats.wgini
                                  };
            end
            table_out = make_table(statistics, params);
        end

        function table_out = table1panelC(params_in, results)
            params = filter_param_names(params_in, tables.FinalTables.table1includes);
            statistics = cell(numel(params), 1);
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  results(ip).stats.mpcs(4).annual
                                    results(ip).stats.mpcs(6).annual
                                    results(ip).stats.mpcs(4).quarterly
                                    results(ip).stats.mpcs(6).quarterly
                                  };
            end
            table_out = make_table(statistics, params);
        end

        function table_out = table1panelD(params_in, results)
            params = filter_param_names(params_in, tables.FinalTables.table1includes);
            statistics = cell(numel(params), 1);
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  results(ip).stats.mpcs(1).annual
                                    results(ip).stats.mpcs(2).annual
                                    results(ip).stats.mpcs(3).annual
                                    results(ip).stats.mpcs(1).quarterly
                                    results(ip).stats.mpcs(2).quarterly
                                    results(ip).stats.mpcs(3).quarterly
                                  };
            end
            table_out = make_table(statistics, params);
        end

        function table_out = table2header(params_in, results)
            params = filter_param_names(params_in, tables.FinalTables.table2includes);
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
            table_out = make_table(statistics, params);
        end

        function table_out = table2panelABC(params_in, results)
            params = filter_param_names(params_in, tables.FinalTables.table2includes);
            table_out = cell(1, 3);

            for kk = 1:3
                statistics = cell(numel(params), 1);

                for ii = 1:numel(params)
                    ip = params(ii).index;
                    statistics{ii} = {  results(ip).stats.decomp_norisk.term2(kk)
                                        results(ip).stats.decomp_norisk.term3(kk)
                                        results(ip).stats.decomp_norisk.term4(kk)
                                      };
                end
                table_out{kk} = make_table(statistics, params);
            end
        end

        function table_out = table2panelD(params_in, results)
            params = filter_param_names(params_in, tables.FinalTables.table2includes);
            statistics = cell(numel(params), 1);
            
            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  results(ip).stats.decomp_RA.Em1_less_mRA
                                    results(ip).stats.decomp_RA.term1
                                    results(ip).stats.decomp_RA.term2
                                    results(ip).stats.decomp_RA.term3
                                  };
            end
            table_out = make_table(statistics, params);
        end

        function table_out = experiment_table_header(params_in, results, table)
            if (table == 3)
                params = filter_param_group(params_in, tables.FinalTables.table3includes);
            end

            import statistics.Statistics.sfill

            for ii = 1:numel(params)
                ip = params(ii).index;
                variable_value = sfill(params(ii).tex_header_value, 'Value', 2);
                statistics{ii} = {  variable_value
                                    results(ip).stats.mpcs(5).quarterly
                                    results(ip).stats.mpcs(5).annual
                                    results(ip).stats.beta_A
                                  };
            end
            table_out = make_table(statistics, params, true);
        end

        function table_out = experiment_table_panelA(params_in, comparison_decomps, tableno)
            if (tableno == 3)
                params = filter_param_group(params_in, tables.FinalTables.table3includes);
            end

            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  comparison_decomps(ip).Em1_less_Em0
                                    comparison_decomps(ip).term1
                                    comparison_decomps(ip).term2
                                    comparison_decomps(ip).term2a(2)
                                    comparison_decomps(ip).term2b(2)
                                    comparison_decomps(ip).term3
                                  };
            end
            table_out = make_table(statistics, params, true);
        end

        function table_out = experiment_table_panelA2(params_in, comparison_decomps, tableno)
            if (tableno == 3)
                params = filter_param_group(params_in, tables.FinalTables.table3includes);
            end

            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  comparison_decomps(ip).term1_pct
                                    comparison_decomps(ip).term2_pct
                                    comparison_decomps(ip).term2a_pct(2)
                                    comparison_decomps(ip).term2b_pct(2)
                                    comparison_decomps(ip).term3_pct
                                  };
            end
            table_out = make_table(statistics, params, true);
        end

        function table_out = experiment_table_panelB(params_in, results, tableno)
            if (tableno == 3)
                params = filter_param_group(params_in, tables.FinalTables.table3includes);
            end

            for ii = 1:numel(params)
                ip = params(ii).index;
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
            end
            table_out = make_table(statistics, params, true);
        end

        function table_out = experiment_table_panelC(params_in, results, tableno)
            if (tableno == 3)
                params = filter_param_group(params_in, tables.FinalTables.table3includes);
            end

            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  results(ip).stats.mpcs(4).quarterly
                                    results(ip).stats.mpcs(6).quarterly
                                  };
            end
            table_out = make_table(statistics, params, true);
        end

        function table_out = experiment_table_panelD(params_in, results, tableno)
            if (tableno == 3)
                params = filter_param_group(params_in, tables.FinalTables.table3includes);
            end

            for ii = 1:numel(params)
                ip = params(ii).index;
                statistics{ii} = {  results(ip).stats.mpcs(1).quarterly
                                    results(ip).stats.mpcs(2).quarterly
                                    results(ip).stats.mpcs(3).quarterly
                                  };
            end
            table_out = make_table(statistics, params, true);
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

function params = filter_param_group(params_in, includes)
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
            if jj == 1
                params = params_in(ii);
            else
                params(jj) = params_in(ii);
            end
            jj = jj + 1;
        end
    end
end

function table_out = make_table(statistics, params, experiment)
    if nargin < 3
        experiment = false;
    end

    n = numel(statistics);
    
    for ii = 1:n
        vars{ii} = get_values(statistics{ii});

        if experiment
            varnames{ii} = char(params(ii).tex_header + string(sprintf('__v%d__', ii)));
        else
            varnames{ii} = params(ii).name;
        end
    end
    vars{n+1} = get_precision(statistics{1});
    varnames{n+1} = 'decimals';

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