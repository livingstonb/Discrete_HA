classdef FinalTables
    
    properties (Constant)
        table1includes = {'Quarterly', 'Annual'};
        table2includes = {'Quarterly', 'Annual'};
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

            headerpath = fullfile(dirpath, 'table2_header.xlsx');
            writetable(header, headerpath, 'WriteRowNames', true);
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
            params = filter_param_names(params_in, tables.FinalTables.table1includes);
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

function table_out = make_table(statistics, params)
    n = numel(statistics);
    
    for ii = 1:n
        vars{ii} = get_values(statistics{ii});
        varnames{ii} = params(ii).name;
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