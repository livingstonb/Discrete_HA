function [T_annual,T_quarter] = create_table(params,direct_results,decomps,checks,exceptions,decomp2)
    % Rownames
    rows = {'Specification'
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
            '    MPCS OUT OF 0.01 MEAN ANN INC)'
            'Mean annual MPC (size = -1e-5)'
            'Mean annual MPC (size = -0.01)'
            'Mean annual MPC (size = -0.1)'
            'Mean annual MPC (size = 1e-5)'
            'Mean annual MPC (size = 0.01)'
            'Mean annual MPC (size = 0.1)'
            'Mean quarterly MPC (size = -1e-5)'
            'Mean quarterly MPC (size = -0.01)'
            'Mean quarterly MPC (size = -0.1)'
            'Mean quarterly MPC (size = 1e-5)'
            'Mean quarterly MPC (size = 0.01)'
            'Mean quarterly MPC (size = 0.1)'
            '    DECOMPS (MPC OUT OF 0.01 MEAN ANN INC)'
            'Decomp of Em0 around 0, RA MPC'
            'Decomp of Em0 around 0, HtM Effect'
            'Decomp of Em0 around 0, Non-HtM, constraint'
            'Decomp of Em0 around 0, Non-HtM, inc risk'
            'Decomp of Em0 around 0.01, RA MPC'
            'Decomp of Em0 around 0.01, HtM Effect'
            'Decomp of Em0 around 0.01, Non-HtM, constraint'
            'Decomp of Em0 around 0.01, Non-HtM, inc risk'
            'Decomp of Em0 around 0.05, RA MPC'
            'Decomp of Em0 around 0.05, HtM Effect'
            'Decomp of Em0 around 0.05, Non-HtM, constraint'
            'Decomp of Em0 around 0.05, Non-HtM, inc risk'
            'Em1 - Em0'
            'Decomp of Em1-Em0, around 0, effect of MPC fcn'
            'Decomp of Em1-Em0, around 0, effect of distr'
            'Decomp of Em1-Em0, around 0, interaction'
            'Decomp of Em1-Em0, around 0.01, effect of MPC fcn'
            'Decomp of Em1-Em0, around 0.01, effect of distr'
            'Decomp of Em1-Em0, around 0.01, interaction'
            'Decomp of Em1-Em0, around 0.05, effect of MPC fcn'
            'Decomp of Em1-Em0, around 0.05, effect of distr'
            'Decomp of Em1-Em0, around 0.05, interaction'
            'Failed one or more checks'
            };
    Nrows = numel(rows) - 1;

    % Iterate over frequency
    for ifreq = [1 4]
        this_freq = find([params.freq]==ifreq & ~cellfun('isempty', direct_results));
        params_freq = params(this_freq);
        names = {params_freq.name};
        tablearray = zeros(Nrows,numel(params_freq));
        
        % Iterate over frequency ifreq
        ncolumn = 1;
        for ip = this_freq
            p = params(ip);

            if numel(exceptions{ip}) == 1
                % Exception was thrown for this parameterization
                column = NaN(Nrows,1);
            elseif numel(checks{ip}) > 0
                if sum(ismember({'NoEGPConv','NoBetaConv'},checks{ip})) > 0
                    % Critical code failure
                    column = NaN(Nrows,1);
                end
            else
                % Annual and quarterly MPCs
                if p.freq == 1
                    mpcs_A = direct_results{ip}.avg_mpc1_sim(:);
                    mpcs_Q = NaN(numel(p.mpcfrac),1);
                else
                    mpcs_A = direct_results{ip}.avg_mpc4_sim(:);
                    mpcs_Q = direct_results{ip}.avg_mpc1_sim(:);
                end
                
                % decomposition1
                dec1 = [[decomps{ip}.term1]
                        [decomps{ip}.term2]
                        [decomps{ip}.term3]
                        [decomps{ip}.term4]];
                
                % decomposition2
                dec2 = [[decomp2{ip}.term1]
                        [decomp2{ip}.term2]
                        [decomp2{ip}.term3]];
                
                column = [
                    direct_results{ip}.beta_annualized      % Annualized beta
                    direct_results{ip}.mean_grossy_A        % Mean annual gross labor income
                    direct_results{ip}.stdev_loggrossy_A    % Stdev log annual gross income
                    direct_results{ip}.stdev_lognety_A      % Stdev log annual net income
                    NaN
                    direct_results{ip}.mean_a               % Mean assets
                    direct_results{ip}.constrained(:)       % Fraction with a < eps * mean ann gross inc
                    direct_results{ip}.wpercentiles(:)      % Wealth percentiles
                    direct_results{ip}.top10share           % Top 10% wealth share
                    direct_results{ip}.top1share            % Top 1% wealth share
                    direct_results{ip}.wealthgini           % Gini coefficient
                    NaN
                    mpcs_A(:)                               % Annual MPCs
                    mpcs_Q(:)                               % Quarterly MPCs (if freq = 4)
                    NaN
                    dec1(:)                                 % Decomp1
                    decomp2{ip}(1).Em1_less_Em0             % Decomposition2
                    dec2(:)
                    numel(checks{ip})>0];                
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