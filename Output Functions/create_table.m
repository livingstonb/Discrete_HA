function [T_annual,T_quarter] = create_table(params,results,...
                                            decomps,checks,decomp2)
    %% Rownames
    rows = {'Specification'
            'Lookup Index'
            'Beta (Annualized)'
            'Mean gross annual income'
            'Stdev log annual gross income'
            'Stdev log annual net income'
            '____WEALTH STATISTICS'
            'Mean assets'
            'Fraction with a == 0'
            'Fraction with a <= 0.5% mean ann gross lab inc'
            'Fraction with a <= 1% mean ann gross lab inc'
            'Fraction with a <= 2% mean ann gross lab inc'
            'Fraction with a <= 5% mean ann gross lab inc'
            'Fraction with a <= 10% mean ann gross lab inc'
            'Fraction with a <= 1/6 own quarterly income'
            'Fraction with a <= 1/12 own quarterly income'
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
            '____MPCS OUT OF FRACTION MEAN ANN INC (SHOCK 0.01)'
            'IMPC(1,1), shock=0.01'
            'IMPC(1,2), shock=0.01'
            'IMPC(1,3), shock=0.01'
            'IMPC(1,4), shock=0.01'
            'IMPC(1,1-4), shock=0.01'
            'IMPC(1,5-9), shock=0.01'
            'IMPC(1,10-13), shock=0.01'
            'IMPC(1,13-16), shock=0.01'
            'IMPC(2,1), shock=0.01'
            'IMPC(2,2), shock=0.01'
            'IMPC(2,3), shock=0.01'
            'IMPC(2,4), shock=0.01'
            'IMPC(2,1-4), shock=0.01'
            'IMPC(2,5-9), shock=0.01'
            'IMPC(2,10-13), shock=0.01'
            'IMPC(2,13-16), shock=0.01'
            'IMPC(3,1), shock=0.01'
            'IMPC(3,2), shock=0.01'
            'IMPC(3,3), shock=0.01'
            'IMPC(3,4), shock=0.01'
            'IMPC(3,1-4), shock=0.01'
            'IMPC(3,5-9), shock=0.01'
            'IMPC(3,10-13), shock=0.01'
            'IMPC(3,13-16), shock=0.01'
            'IMPC(4,1), shock=0.01'
            'IMPC(4,2), shock=0.01'
            'IMPC(4,3), shock=0.01'
            'IMPC(4,4), shock=0.01'
            'IMPC(4,1-4), shock=0.01'
            'IMPC(4,5-9), shock=0.01'
            'IMPC(4,10-13), shock=0.01'
            'IMPC(4,13-16), shock=0.01'
            'IMPC(5,1-4), shock=0.01'
            'IMPC(5,5-9), shock=0.01'
            'IMPC(5,10-13), shock=0.01'
            'IMPC(5,13-16), shock=0.01'
            'IMPC(9,1-4), shock=0.01'
            'IMPC(9,5-9), shock=0.01'
            'IMPC(9,10-13), shock=0.01'
            'IMPC(9,13-16), shock=0.01'
            'IMPC(13,1-4), shock=0.01'
            'IMPC(13,5-9), shock=0.01'
            'IMPC(13,10-13), shock=0.01'
            'IMPC(13,13-16), shock=0.01'
            '____MPCS OUT OF FRACTION MEAN ANN INC (SHOCK 0.05)'
            'IMPC(1,1), shock=0.05'
            'IMPC(1,2), shock=0.05'
            'IMPC(1,3), shock=0.05'
            'IMPC(1,4), shock=0.05'
            'IMPC(1,1-4), shock=0.05'
            'IMPC(1,5-9), shock=0.05'
            'IMPC(1,10-13), shock=0.05'
            'IMPC(1,13-16), shock=0.05'
            'IMPC(2,1), shock=0.05'
            'IMPC(2,2), shock=0.05'
            'IMPC(2,3), shock=0.05'
            'IMPC(2,4), shock=0.05'
            'IMPC(2,1-4), shock=0.05'
            'IMPC(2,5-9), shock=0.05'
            'IMPC(2,10-13), shock=0.05'
            'IMPC(2,13-16), shock=0.05'
            'IMPC(3,1), shock=0.05'
            'IMPC(3,2), shock=0.05'
            'IMPC(3,3), shock=0.05'
            'IMPC(3,4), shock=0.05'
            'IMPC(3,1-4), shock=0.05'
            'IMPC(3,5-9), shock=0.05'
            'IMPC(3,10-13), shock=0.05'
            'IMPC(3,13-16), shock=0.05'
            'IMPC(4,1), shock=0.05'
            'IMPC(4,2), shock=0.05'
            'IMPC(4,3), shock=0.05'
            'IMPC(4,4), shock=0.05'
            'IMPC(4,1-4), shock=0.05'
            'IMPC(4,5-9), shock=0.05'
            'IMPC(4,10-13), shock=0.05'
            'IMPC(4,13-16), shock=0.05'
            'IMPC(5,1-4), shock=0.05'
            'IMPC(5,5-9), shock=0.05'
            'IMPC(5,10-13), shock=0.05'
            'IMPC(5,13-16), shock=0.05'
            'IMPC(9,1-4), shock=0.05'
            'IMPC(9,5-9), shock=0.05'
            'IMPC(9,10-13), shock=0.05'
            'IMPC(9,13-16), shock=0.05'
            'IMPC(13,1-4), shock=0.05'
            'IMPC(13,5-9), shock=0.05'
            'IMPC(13,10-13), shock=0.05'
            'IMPC(13,13-16), shock=0.05'
            '____DECOMP OF IMPC(1,1) (MPC OUT OF FRACTION MEAN ANN INC)'
            '(A/Q) Decomp of Em0 around 0, RA MPC'
            '(A/Q) Decomp of Em0 around 0, HtM Effect'
            '(A/Q) Decomp of Em0 around 0, Non-HtM, constraint'
            '(A/Q) Decomp of Em0 around 0, Non-HtM, inc risk'
            '(A/Q) Decomp of Em0 around 0.01, RA MPC'
            '(A/Q) Decomp of Em0 around 0.01, HtM Effect'
            '(A/Q) Decomp of Em0 around 0.01, Non-HtM, constraint'
            '(A/Q) Decomp of Em0 around 0.01, Non-HtM, inc risk'
            '(A/Q) Decomp of Em0 around 0.05, RA MPC'
            '(A/Q) Decomp of Em0 around 0.05, HtM Effect'
            '(A/Q) Decomp of Em0 around 0.05, Non-HtM, constraint'
            '(A/Q) Decomp of Em0 around 0.05, Non-HtM, inc risk'
            '____DECOMP OF ANNUAL EM1-EM0 (MPC OUT OF 0.01 MEAN ANN INC)'
            '(Annual) Em1 - Em0'
            '(Annual) Decomp of Em1-Em0, effect of MPC fcn'
            '(Annual) Decomp of Em1-Em0, effect of distr'
            '(Annual) Decomp of Em1-Em0, interaction'
            '(Annual) Decomp of the distr effect around 0, HtM households'
            '(Annual) Decomp of the distr effect around 0, non-HtM households'
            '(Annual) Decomp of the distr effect around 0.01, HtM households'
            '(Annual) Decomp of the distr effect around 0.01, non-HtM households'
            '(Annual) Decomp of the distr effect around 0.05, HtM households'
            '(Annual) Decomp of the distr effect around 0.05, non-HtM households'
            '____DECOMP OF QUARTERLY EM1-EM0 (MPC OUT OF 0.01 MEAN ANN INC)'
            '(Quarterly) Em1 - Em0'
            '(Quarterly) Decomp of Em1-Em0, effect of MPC fcn'
            '(Quarterly) Decomp of Em1-Em0, effect of distr'
            '(Quarterly) Decomp of Em1-Em0, interaction'
            '(Quarterly) Decomp of the distr effect around 0, HtM households'
            '(Quarterly) Decomp of the distr effect around 0, non-HtM households'
            '(Quarterly) Decomp of the distr effect around 0.01, HtM households'
            '(Quarterly) Decomp of the distr effect around 0.01, non-HtM households'
            '(Quarterly) Decomp of the distr effect around 0.05, HtM households'
            '(Quarterly) Decomp of the distr effect around 0.05, non-HtM households'
            'Failed one or more checks'
            };
    Nrows = numel(rows) - 1;

    %% Iterate over frequency
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
            if numel(checks{ip}) > 0
                if sum(ismember({'NoEGPConv','NoBetaConv'},checks{ip})) > 0
                    % Critical code failure
                    NaNcol = true;
                end
            end
            if numel(checks{ip}) == 1
                if ismember('EXCEPTION_THROWN',checks{ip})
                    % Exception was thrown for this parameterization
                    NaNcol = true;
                end
            end

            if NaNcol == true
                    column = [p.index;NaN(Nrows-1,1)];
            else
                
                % decomposition1
                dec1 = [[decomps{ip}.term1]
                        [decomps{ip}.term2]
                        [decomps{ip}.term3]
                        [decomps{ip}.term4]];
                    
                % decomposition2
                dec2_mpc1 = [decomp2(ip).mpc1_Em1_less_Em0
                                decomp2(ip).mpc1_term1                       
                                decomp2(ip).mpc1_term2
                                decomp2(ip).mpc1_term3
                                decomp2(ip).mpc1_term3a(1)   
                                decomp2(ip).mpc1_term3b(1)
                                decomp2(ip).mpc1_term3a(2)   
                                decomp2(ip).mpc1_term3b(2)
                                decomp2(ip).mpc1_term3a(3)   
                                decomp2(ip).mpc1_term3b(3)];
                dec2_mpc4 = [decomp2(ip).mpc4_Em1_less_Em0
                                decomp2(ip).mpc4_term1                       
                                decomp2(ip).mpc4_term2
                                decomp2(ip).mpc4_term3
                                decomp2(ip).mpc4_term3a(1)   
                                decomp2(ip).mpc4_term3b(1)
                                decomp2(ip).mpc4_term3a(2)   
                                decomp2(ip).mpc4_term3b(2)
                                decomp2(ip).mpc4_term3a(3)   
                                decomp2(ip).mpc4_term3b(3)];
                if p.freq == 1
                    dec2_A = dec2_mpc1;
                    dec2_Q = NaN(10,1);
                else
                    dec2_A = dec2_mpc4;
                    dec2_Q = dec2_mpc1;
                end

                column = [
                    p.index
                    results(ip).direct.beta_annualized      % Annualized beta
                    results(ip).direct.mean_grossy_A        % Mean annual gross labor income
                    results(ip).direct.stdev_loggrossy_A    % Stdev log annual gross income
                    results(ip).direct.stdev_lognety_A      % Stdev log annual net income
                    NaN
                    results(ip).direct.mean_a               % Mean assets
                    results(ip).direct.constrained(:)       % Fraction with a < eps * mean ann gross inc
                    results(ip).direct.HtM_one_sixth_Q          % Fraction with a < 1/6 quarterly income
                    results(ip).direct.HtM_one_twelfth_Q        % Fraction with a < 1/12 quarterly income
                    results(ip).direct.wpercentiles(:)      % Wealth percentiles
                    results(ip).direct.top10share           % Top 10% wealth share
                    results(ip).direct.top1share            % Top 1% wealth share
                    results(ip).direct.wealthgini           % Gini coefficient
                    NaN
                    results(ip).direct.mpcs01.avg_s_t{1,1}        % IMPC(1,1)
                    results(ip).direct.mpcs01.avg_s_t{1,2}        % IMPC(1,2)
                    results(ip).direct.mpcs01.avg_s_t{1,3}        % IMPC(1,3)
                    results(ip).direct.mpcs01.avg_s_t{1,4}        % IMPC(1,4)
                    results(ip).direct.mpcs01.avg_1_1to4          % IMPC(1,1-4)
                    results(ip).direct.mpcs01.avg_1_5to8          % IMPC(1,5-8)
                    results(ip).direct.mpcs01.avg_1_9to12         % IMPC(1,9-12)
                    results(ip).direct.mpcs01.avg_1_13to16        % IMPC(1,13-16)
                    results(ip).direct.mpcs01.avg_s_t{2,1}        % IMPC(2,1)
                    results(ip).direct.mpcs01.avg_s_t{2,2}        % IMPC(2,2)
                    results(ip).direct.mpcs01.avg_s_t{2,3}        % IMPC(2,3)
                    results(ip).direct.mpcs01.avg_s_t{2,4}        % IMPC(2,4)
                    results(ip).direct.mpcs01.avg_2_1to4          % IMPC(2,1-4)
                    results(ip).direct.mpcs01.avg_2_5to8          % IMPC(2,5-8)
                    results(ip).direct.mpcs01.avg_2_9to12         % IMPC(2,9-12)
                    results(ip).direct.mpcs01.avg_2_13to16        % IMPC(2,13-16)
                    results(ip).direct.mpcs01.avg_s_t{3,1}        % IMPC(3,1)
                    results(ip).direct.mpcs01.avg_s_t{3,2}        % IMPC(3,2)
                    results(ip).direct.mpcs01.avg_s_t{3,3}        % IMPC(3,3)
                    results(ip).direct.mpcs01.avg_s_t{3,4}        % IMPC(3,4)
                    results(ip).direct.mpcs01.avg_3_1to4          % IMPC(3,1-4)
                    results(ip).direct.mpcs01.avg_3_5to8          % IMPC(3,5-8)
                    results(ip).direct.mpcs01.avg_3_9to12         % IMPC(3,9-12)
                    results(ip).direct.mpcs01.avg_3_13to16        % IMPC(3,13-16)
                    results(ip).direct.mpcs01.avg_s_t{4,1}        % IMPC(4,1)
                    results(ip).direct.mpcs01.avg_s_t{4,2}        % IMPC(4,2)
                    results(ip).direct.mpcs01.avg_s_t{4,3}        % IMPC(4,3)
                    results(ip).direct.mpcs01.avg_s_t{4,4}        % IMPC(4,4)
                    results(ip).direct.mpcs01.avg_4_1to4          % IMPC(4,1-4)
                    results(ip).direct.mpcs01.avg_4_5to8          % IMPC(4,5-8)
                    results(ip).direct.mpcs01.avg_4_9to12         % IMPC(4,9-12)
                    results(ip).direct.mpcs01.avg_4_13to16        % IMPC(4,13-16)
                    results(ip).direct.mpcs01.avg_5_1to4          % IMPC(5,13-16)
                    results(ip).direct.mpcs01.avg_5_5to8          % IMPC(5,5-8)
                    results(ip).direct.mpcs01.avg_5_9to12         % IMPC(5,9-12)
                    results(ip).direct.mpcs01.avg_5_13to16        % IMPC(5,13-16)
                    results(ip).direct.mpcs01.avg_9_1to4          % IMPC(9,1-4)
                    results(ip).direct.mpcs01.avg_9_5to8          % IMPC(9,5-8)
                    results(ip).direct.mpcs01.avg_9_9to12         % IMPC(9,9-12)
                    results(ip).direct.mpcs01.avg_9_13to16        % IMPC(9,13-16)
                    results(ip).direct.mpcs01.avg_13_1to4         % IMPC(13,1-4)
                    results(ip).direct.mpcs01.avg_13_5to8         % IMPC(13,5-8)
                    results(ip).direct.mpcs01.avg_13_9to12        % IMPC(13,9-12)
                    results(ip).direct.mpcs01.avg_13_13to16       % IMPC(13,13-16)
                    NaN
                    results(ip).direct.mpcs05.avg_s_t{1,1}        % IMPC(1,1)
                    results(ip).direct.mpcs05.avg_s_t{1,2}        % IMPC(1,2)
                    results(ip).direct.mpcs05.avg_s_t{1,3}        % IMPC(1,3)
                    results(ip).direct.mpcs05.avg_s_t{1,4}        % IMPC(1,4)
                    results(ip).direct.mpcs05.avg_1_1to4          % IMPC(1,1-4)
                    results(ip).direct.mpcs05.avg_1_5to8          % IMPC(1,5-8)
                    results(ip).direct.mpcs05.avg_1_9to12         % IMPC(1,9-12)
                    results(ip).direct.mpcs05.avg_1_13to16        % IMPC(1,13-16)
                    results(ip).direct.mpcs05.avg_s_t{2,1}        % IMPC(2,1)
                    results(ip).direct.mpcs05.avg_s_t{2,2}        % IMPC(2,2)
                    results(ip).direct.mpcs05.avg_s_t{2,3}        % IMPC(2,3)
                    results(ip).direct.mpcs05.avg_s_t{2,4}        % IMPC(2,4)
                    results(ip).direct.mpcs05.avg_2_1to4          % IMPC(2,1-4)
                    results(ip).direct.mpcs05.avg_2_5to8          % IMPC(2,5-8)
                    results(ip).direct.mpcs05.avg_2_9to12         % IMPC(2,9-12)
                    results(ip).direct.mpcs05.avg_2_13to16        % IMPC(2,13-16)
                    results(ip).direct.mpcs05.avg_s_t{3,1}        % IMPC(3,1)
                    results(ip).direct.mpcs05.avg_s_t{3,2}        % IMPC(3,2)
                    results(ip).direct.mpcs05.avg_s_t{3,3}        % IMPC(3,3)
                    results(ip).direct.mpcs05.avg_s_t{3,4}        % IMPC(3,4)
                    results(ip).direct.mpcs05.avg_3_1to4          % IMPC(3,1-4)
                    results(ip).direct.mpcs05.avg_3_5to8          % IMPC(3,5-8)
                    results(ip).direct.mpcs05.avg_3_9to12         % IMPC(3,9-12)
                    results(ip).direct.mpcs05.avg_3_13to16        % IMPC(3,13-16)
                    results(ip).direct.mpcs05.avg_s_t{4,1}        % IMPC(4,1)
                    results(ip).direct.mpcs05.avg_s_t{4,2}        % IMPC(4,2)
                    results(ip).direct.mpcs05.avg_s_t{4,3}        % IMPC(4,3)
                    results(ip).direct.mpcs05.avg_s_t{4,4}        % IMPC(4,4)
                    results(ip).direct.mpcs05.avg_4_1to4          % IMPC(4,1-4)
                    results(ip).direct.mpcs05.avg_4_5to8          % IMPC(4,5-8)
                    results(ip).direct.mpcs05.avg_4_9to12         % IMPC(4,9-12)
                    results(ip).direct.mpcs05.avg_4_13to16        % IMPC(4,13-16)
                    results(ip).direct.mpcs05.avg_5_1to4          % IMPC(5,1-4)
                    results(ip).direct.mpcs05.avg_5_5to8          % IMPC(5,5-8)
                    results(ip).direct.mpcs05.avg_5_9to12         % IMPC(5,9-12)
                    results(ip).direct.mpcs05.avg_5_13to16        % IMPC(5,13-16)
                    results(ip).direct.mpcs05.avg_9_1to4          % IMPC(9,1-4)
                    results(ip).direct.mpcs05.avg_9_5to8          % IMPC(9,5-8)
                    results(ip).direct.mpcs05.avg_9_9to12         % IMPC(9,9-12)
                    results(ip).direct.mpcs05.avg_9_13to16        % IMPC(9,13-16)
                    results(ip).direct.mpcs05.avg_13_1to4         % IMPC(13,1-4)
                    results(ip).direct.mpcs05.avg_13_5to8         % IMPC(13,5-8)
                    results(ip).direct.mpcs05.avg_13_9to12        % IMPC(13,9-12)
                    results(ip).direct.mpcs05.avg_13_13to16       % IMPC(13,13-16)
                    NaN
                    dec1(:)                                 % Decomp1
                    NaN
                    dec2_A                                  % Decomposition2
                    NaN
                    dec2_Q
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