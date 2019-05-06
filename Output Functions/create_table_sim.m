function [T_annual,T_quarter] = create_table_sim(params,results,...
                                            decomps)
    %% Rownames
    rows = {'Specification'
            'Lookup Index'
            'Beta (Annualized)'
            'Mean gross annual income'
            'Stdev log annual gross income'
            'Stdev log annual net income'
            '    WEALTH STATISTICS'
            'Mean assets'
            'Fraction with a == 0'
            'Fraction with a <= 0.5% mean ann gross lab inc'
            'Fraction with a <= 1% mean ann gross lab inc'
            'Fraction with a <= 2% mean ann gross lab inc'
            'Fraction with a <= 5% mean ann gross lab inc'
            'Fraction with a <= 10% mean ann gross lab inc'
            'Wealth, 10th percentile'
            'Wealth, 25th percentile'
            'Wealth, 50th percentile'
            'Wealth, 75th percentile'
            'Wealth, 90th percentile'
            'Wealth, 95th percentile'
            'Wealth, 99th percentile'
            'Wealth, 99.9th percentile'
            'Wealth, top 10% share'
            'Wealth, top 1% share'
            'Gini coefficient'
            '    MPCS OUT OF FRACTION MEAN ANN INC'
            'IMPC(1,1) (size = -1e-5)'
            'IMPC(1,1) (size = -0.01)'
            'IMPC(1,1) (size = -0.1)'
            'IMPC(1,1) (size = 1e-5)'
            'IMPC(1,1) (size = 0.01)'
            'IMPC(1,1) (size = 0.1)'
            'IMPC(1,2) (size = -1e-5)'
            'IMPC(1,2) (size = -0.01)'
            'IMPC(1,2) (size = -0.1)'
            'IMPC(1,2) (size = 1e-5)'
            'IMPC(1,2) (size = 0.01)'
            'IMPC(1,2) (size = 0.1)'
            'IMPC(1,3) (size = -1e-5)'
            'IMPC(1,3) (size = -0.01)'
            'IMPC(1,3) (size = -0.1)'
            'IMPC(1,3) (size = 1e-5)'
            'IMPC(1,3) (size = 0.01)'
            'IMPC(1,3) (size = 0.1)'
            'IMPC(1,4) (size = -1e-5)'
            'IMPC(1,4) (size = -0.01)'
            'IMPC(1,4) (size = -0.1)'
            'IMPC(1,4) (size = 1e-5)'
            'IMPC(1,4) (size = 0.01)'
            'IMPC(1,4) (size = 0.1)'
            'IMPC(1,1-4) (size = -1e-5)'
            'IMPC(1,1-4) (size = -0.01)'
            'IMPC(1,1-4) (size = -0.1)'
            'IMPC(1,1-4) (size = 1e-5)'
            'IMPC(1,1-4) (size = 0.01)'
            'IMPC(1,1-4) (size = 0.1)'
            'IMPC(1,5-8) (size = -1e-5)'
            'IMPC(1,5-8) (size = -0.01)'
            'IMPC(1,5-8) (size = -0.1)'
            'IMPC(1,5-8) (size = 1e-5)'
            'IMPC(1,5-8) (size = 0.01)'
            'IMPC(1,5-8) (size = 0.1)'
            'IMPC(1,9-12) (size = -1e-5)'
            'IMPC(1,9-12) (size = -0.01)'
            'IMPC(1,9-12) (size = -0.1)'
            'IMPC(1,9-12) (size = 1e-5)'
            'IMPC(1,9-12) (size = 0.01)'
            'IMPC(1,9-12) (size = 0.1)'
            'IMPC(1,13-16) (size = -1e-5)'
            'IMPC(1,13-16) (size = -0.01)'
            'IMPC(1,13-16) (size = -0.1)'
            'IMPC(1,13-16) (size = 1e-5)'
            'IMPC(1,13-16) (size = 0.01)'
            'IMPC(1,13-16) (size = 0.1)'
            'Failed one or more checks'
            };
    Nrows = numel(rows) - 1;

    % Iterate over frequency
    for ifreq = [1 4]
        this_freq = find([params.freq]==ifreq);
        if isempty(this_freq)
            if ifreq == 1
                T_annual = [];
            else
                T_quarter = [];
            end
            continue
        end
        params_freq = params(this_freq);
        names = {params_freq.name};
        tablearray = zeros(Nrows,numel(params_freq));
        
        % Iterate over frequency ifreq
        ncolumn = 1;
        for ip = this_freq
            p = params(ip);
            
            % Check if column of NaNs must be used
            NaNcol = false;
            if numel(results(ip).checks) > 0
                if sum(ismember({'NoEGPConv','NoBetaConv'},results(ip).checks)) > 0
                    % Critical code failure
                    NaNcol = true;
                end
            end
            if numel(results(ip).checks) == 1
                if ismember('EXCEPTION_THROWN',results(ip).checks)
                    % Exception was thrown for this parameterization
                    NaNcol = true;
                end
            end

            if NaNcol == true
                column = [p.index;NaN(Nrows-1,1)];
            else
                
                column = [
                    p.index
                    results(ip).direct.beta_annualized      % Annualized beta
                    results(ip).sim.mean_grossy_A        % Mean annual gross labor income
                    sqrt(results(ip).sim.var_loggrossy_A)    % Stdev log annual gross income
                    sqrt(results(ip).sim.var_lognety_A)      % Stdev log annual net income
                    NaN
                    results(ip).sim.mean_a               % Mean assets
                    results(ip).sim.constrained(:)       % Fraction with a < eps * mean ann gross inc
                    results(ip).sim.wpercentiles(:)      % Wealth percentiles
                    results(ip).sim.top10share           % Top 10% wealth share
                    results(ip).sim.top1share            % Top 1% wealth share
                    results(ip).sim.wealthgini_A           % Gini coefficient
                    NaN
                    results(ip).sim.mpcs.avg_1_1(:)          % IMPC(1,1)
                    results(ip).sim.mpcs.avg_1_2(:)          % IMPC(1,2)
                    results(ip).sim.mpcs.avg_1_3(:)          % IMPC(1,3)
                    results(ip).sim.mpcs.avg_1_4(:)          % IMPC(1,4)
                    results(ip).sim.mpcs.avg_1_1to4(:)          % IMPC(1,1-4)
                    results(ip).sim.mpcs.avg_1_5to8(:)          % IMPC(1,5-8)
                    results(ip).sim.mpcs.avg_1_9to12(:)          % IMPC(1,9-12)
                    results(ip).sim.mpcs.avg_1_13to16(:)          % IMPC(1,13-16)
                    numel(results(ip).checks)>0];                
            end

            % Add this column to table
            tablearray(:,ncolumn) = column;
            ncolumn = ncolumn + 1;
        end
        
        tablearray = num2cell(tablearray);
        tablearray = [names;tablearray];
        T = array2table(tablearray);
        T.Properties.RowNames = rows;
        if ifreq == 1
            T_annual = T;
        else
            T_quarter = T;
        end
    end

end