function [T_annual,T_quarter] = create_table(params,results,...
                                            decomps,decomp2,decomp3)
    %% Rownames
    rows = {'Specification'
            'Lookup Index'
            'Beta (Annualized)'
            'Mean gross annual income'
            'Stdev log annual gross income'
            'Stdev log annual net income'
            '____WEALTH CONSTRAINTS'
            'Mean assets'
            'Mean saving'
            'Fraction with s == 0'
            'Fraction with a == 0'
            'Fraction with a <= 0.5% mean ann gross lab inc'
            'Fraction with a <= 1% mean ann gross lab inc'
            'Fraction with a <= 2% mean ann gross lab inc'
            'Fraction with a <= 5% mean ann gross lab inc'
            'Fraction with a <= 10% mean ann gross lab inc'
            'Fraction with a <= 15% mean ann gross lab inc'
            'Fraction with a <= 1/6 own quarterly income'
            'Fraction with a <= 1/12 own quarterly income'
            'Fraction with x <= 1/6 own quarterly income'
            'Fraction with x <= 1/12 own quarterly income'
            '____WEALTH DISTRIBUTION'
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
            '____MPCS OUT OF FRACTION MEAN ANN INC (SHOCK -1e-5)'
            'PERIOD 1 MPC, shock = -1e-5'
            'PERIOD 2 MPC, shock = -1e-5'
            'PERIOD 3 MPC shock = -1e-5'
            'PERIOD 4 MPC shock = -1e-5'
            'FOUR PERIOD MPC, shock = -1e-5'
            '____MPCS OUT OF FRACTION MEAN ANN INC (SHOCK -0.01)'
            'PERIOD 1 MPC, shock = -0.01'
            'PERIOD 2 MPC, shock = -0.01'
            'PERIOD 3 MPC shock = -0.01'
            'PERIOD 4 MPC shock = -0.01'
            'FOUR PERIOD MPC, shock = -0.01'
            '____MPCS OUT OF FRACTION MEAN ANN INC (SHOCK -0.1)'
            'PERIOD 1 MPC, shock = -0.1'
            'PERIOD 2 MPC, shock = -0.1'
            'PERIOD 3 MPC shock = -0.1'
            'PERIOD 4 MPC shock = -0.1'
            'FOUR PERIOD MPC, shock = -0.1'
            '____MPCS OUT OF FRACTION MEAN ANN INC (SHOCK 1e-5)'
            'PERIOD 1 MPC, shock = 1e-5'
            'PERIOD 2 MPC, shock = 1e-5'
            'PERIOD 3 MPC shock = 1-5'
            'PERIOD 4 MPC shock = 1e-5'
            'FOUR PERIOD MPC, shock = 1e-5'
            '____MPCS OUT OF FRACTION MEAN ANN INC (SHOCK 0.01)'
            'PERIOD 1 MPC, shock = 0.01'
            'PERIOD 2 MPC, shock = 0.01'
            'PERIOD 3 MPC shock = 0.01'
            'PERIOD 4 MPC shock = 0.01'
            'FOUR PERIOD MPC, shock = 0.01'
            '____MPCS OUT OF FRACTION MEAN ANN INC (SHOCK 0.1)'
            'PERIOD 1 MPC, shock = 0.1'
            'PERIOD 2 MPC, shock = 0.1'
            'PERIOD 3 MPC shock = 0.1'
            'PERIOD 4 MPC shock = 0.1'
            'FOUR PERIOD MPC, shock = 0.1'
            '____MPC OUT OF NEWS'
            '(A/Q) MPC, shock = 0.01 next period'
            '(A/Q) MPC, shock = 0.1 next period'
            '(A/Q) MPC, shock = 0.01 in four periods'
            '(A/Q) MPC, shock = 0.1 in four periods'
            '(A/Q) E[MPC], negative shock in 8 periods'
            '(A/Q) E[MPC|MPC>0], negative shock in 8 periods'
            '(A/Q) P(MPC<0), negative shock in 8 periods'
            '(A/Q) P(MPC=0), negative shock in 8 periods'
            '(A/Q) P(MPC>0), negative shock in 8 periods'
            '(A/Q) E[MPC], loan for 4 periods'
            '(A/Q) E[MPC|MPC>0], loan for 4 periods'
            '(A/Q) P(MPC<0), loan for 4 periods'
            '(A/Q) P(MPC=0), loan for 4 periods'
            '(A/Q) P(MPC>0), loan for 4 periods'
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
            '____DECOMP OF ONE PERIOD EM1-EM_RA (MPC OUT OF 0.01 MEAN ANN INC)'
            'Em1 - EmRA'
            'Decomp of Em1-EmRA, effect of MPC fcn'
            'Decomp of Em1-EmRA, effect of distr'
            'Decomp of Em1-EmRA, interaction'
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

            if ~results(ip).Finished
                % code failed to run
                column = [p.index;NaN(Nrows-1,1)];
            else
                
%                 if ~NoDecomps
%                     % decomposition1
                    dec1 = [[decomps{ip}.term1]
                            [decomps{ip}.term2]
                            [decomps{ip}.term3]
                            [decomps{ip}.term4]];

                    % decomposition2
                    dec2_mpc1 = [decomp2(ip).mpc1_Em1_less_Em0
                                    decomp2(ip).mpc1_term1                       
                                    decomp2(ip).mpc1_term2
                                    decomp2(ip).mpc1_term3
                                    decomp2(ip).mpc1_term2a(1)   
                                    decomp2(ip).mpc1_term2b(1)
                                    decomp2(ip).mpc1_term2a(2)   
                                    decomp2(ip).mpc1_term2b(2)
                                    decomp2(ip).mpc1_term2a(3)   
                                    decomp2(ip).mpc1_term2b(3)];
                    dec2_mpc4 = [decomp2(ip).mpc4_Em1_less_Em0
                                    decomp2(ip).mpc4_term1                       
                                    decomp2(ip).mpc4_term2
                                    decomp2(ip).mpc4_term3
                                    decomp2(ip).mpc4_term2a(1)   
                                    decomp2(ip).mpc4_term2b(1)
                                    decomp2(ip).mpc4_term2a(2)   
                                    decomp2(ip).mpc4_term2b(2)
                                    decomp2(ip).mpc4_term2a(3)   
                                    decomp2(ip).mpc4_term2b(3)];
                    dec3_mpc1 = [decomp3(ip).mpc1_Em1_less_mRA
                                    decomp3(ip).mpc1_term1
                                    decomp3(ip).mpc1_term2
                                    decomp3(ip).mpc1_term3];
                    dec3_mpc4 = [decomp3(ip).mpc4_Em1_less_mRA
                                    decomp3(ip).mpc4_term1
                                    decomp3(ip).mpc4_term2
                                    decomp3(ip).mpc4_term3];

                    if p.freq == 1
                        dec2_A = dec2_mpc1;
                        dec2_Q = NaN(10,1);
                    else
                        dec2_A = dec2_mpc4;
                        dec2_Q = dec2_mpc1;
                    end

                s_t_mpcs = cell(1,6);

                for i = 1:3
%                     s_t_mpcs{i} = [results(ip).direct.mpcs_sim.avg_1_1(i)
%                                     results(ip).direct.mpcs_sim.avg_1_2(i)
%                                     results(ip).direct.mpcs_sim.avg_1_3(i)
%                                     results(ip).direct.mpcs_sim.avg_1_4(i)
%                                     results(ip).direct.mpcs_sim.avg_1_1to4(i)];
                end

                for i = 1:6
                    s_t_mpcs{i} = [results(ip).direct.mpcs(i).avg_s_t(1,1)          % IMPC(1,1)
                                    results(ip).direct.mpcs(i).avg_s_t(1,2)        % IMPC(1,2)
                                    results(ip).direct.mpcs(i).avg_s_t(1,3)        % IMPC(1,3)
                                    results(ip).direct.mpcs(i).avg_s_t(1,4)        % IMPC(1,4)
                                    results(ip).direct.mpcs(i).avg_1_1to4];          % IMPC(1,1-4)
                end

                column = [
                    p.index
                    results(ip).direct.beta_annualized      % Annualized beta
                    results(ip).direct.mean_grossy_A        % Mean annual gross labor income
                    results(ip).direct.stdev_loggrossy_A    % Stdev log annual gross income
                    results(ip).direct.stdev_lognety_A      % Stdev log annual net income
                    NaN
                    results(ip).direct.mean_a               % Mean assets
                    results(ip).direct.mean_s               % Mean saving
                    results(ip).direct.s0;
                    results(ip).direct.constrained(:)       % Fraction with a < eps * mean ann gross inc
                    results(ip).direct.a_sixth_sim
                    results(ip).direct.a_twelfth_sim       % Fraction with a < 1/12 quarterly income
                    results(ip).direct.x_sixth_sim
                    results(ip).direct.x_twelfth_sim
                    NaN
                    results(ip).direct.wpercentiles(:)      % Wealth percentiles
                    results(ip).direct.top10share           % Top 10% wealth share
                    results(ip).direct.top1share            % Top 1% wealth share
                    results(ip).direct.wealthgini           % Gini coefficient
                    NaN
                    s_t_mpcs{1}
                    NaN
                    s_t_mpcs{2}
                    NaN
                    s_t_mpcs{3}
                    NaN
                    s_t_mpcs{4}
                    NaN
                    s_t_mpcs{5}
                    NaN
                    s_t_mpcs{6}
                    NaN
                    results(ip).direct.mpcs(5).avg_s_t(2,1)
                    results(ip).direct.mpcs(6).avg_s_t(2,1)
                    results(ip).direct.mpcs(5).avg_s_t(5,1)
                    results(ip).direct.mpcs(6).avg_s_t(5,1)
                    results(ip).direct.loss_in_2_years.avg
                    results(ip).direct.loss_in_2_years.mpc_condl
                    results(ip).direct.loss_in_2_years.mpc_neg
                    results(ip).direct.loss_in_2_years.mpc0
                    results(ip).direct.loss_in_2_years.mpc_pos
                    results(ip).direct.loan.avg
                    results(ip).direct.loan.mpc_condl
                    results(ip).direct.loan.mpc_neg
                    results(ip).direct.loan.mpc0
                    results(ip).direct.loan.mpc_pos
                    NaN
                    dec1(:)                                 % Decomp1
                    NaN
                    dec2_A                                  % Decomposition2
                    NaN
                    dec2_Q
                    NaN
                    dec3_mpc1];                             % Decomposition3
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