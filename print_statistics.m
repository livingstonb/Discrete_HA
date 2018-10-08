function print_statistics(direct_results,sim_results,norisk_results,checks,p)
    % This function prints the main results associated with the
    % parameterization structure 'p'
    
    if p.Simulate == 0
        sim = [];
    end
    
    if p.freq == 1
        fprintf('\nANNUAL FREQUENCY\n')
    else
        fprintf('\nQUARTERLY FREQUENCY\n')
    end
    
    %% WEALTH DISTRIBUTION
    
    fprintf('\nWEALTH/SAVINGS: \n')
    
    % Mean wealth
    direct  = sprintf(' %2.3f (direct),',direct_results.mean_a);
    if p.Simulate == 1
        sim = sprintf(' %2.3f (Simulation)',sim_results.mean_a);
    end
    disp(['    Mean wealth:' direct sim]);
    
    % Borrowing constrained    
    for i = 1:numel(p.epsilon)
        label  = sprintf('    Fraction with s <= %3.1f%% mean ann gross inc:',p.epsilon(i)*100);
        direct = sprintf(' %3.3f (Direct)',direct_results.constrained(i));
        if p.Simulate == 1
            sim    = sprintf(', %3.3f (Simulation)',sim_results.constrained(i));
        end
        disp([label direct sim])
    end
    
    % Percentiles
    for i = 1:numel(p.percentiles)
        label = sprintf('    Wealth, %ith percentile:',p.percentiles(i));
        direct = sprintf(' %8.3f (Direct)',direct_results.wpercentiles(i));
        if p.Simulate == 1
            sim = sprintf(', %8.3f (Simulation)',sim_results.wpercentiles(i));
        end
        disp([label direct sim]);
    end
    
    % Top shares
    direct10 = sprintf(' %6.3f (Direct)',direct_results.top10share);
    direct1  = sprintf(' %6.3f (Direct)',direct_results.top1share);
    if p.Simulate == 1
        sim10 = sprintf(', %6.3f (Simulation)',sim_results.top10share);
        sim1 = sprintf(', %6.3f (Simulation)',sim_results.top1share);
    else
        sim10 = []; sim1 = [];
    end
    disp(['    Wealth, top 10% share:' direct10 sim10])
    disp(['    Wealth, top  1% share:' direct1 sim1])
    
    % Gini
    direct = sprintf(' %5.3f (Direct)',direct_results.wealthgini);
    if p.Simulate == 1
        sim = sprintf(' %5.3f (Direct)',sim_results.wealthgini);
    end
    disp(['    Gini coefficient:' direct sim])
    
    
    %% Income Distribution
    fprintf('\nINCOME DISTRIBUTION: \n')
    fprintf('  Cross-sectional variances\n')
    % Variance of log gross labor income
    direct = sprintf(' %2.3f (direct)',direct_results.var_loggrossy);
    if p.Simulate == 1
        sim = sprintf(', %2.3f (Simulation)',sim_results.var_loggrossy);
    end
    disp(['    Var Log Gross Earnings:' direct sim]);
    
    % Variance of log net labor income
    direct = sprintf('   %2.3f (direct)',direct_results.var_lognety);
    if p.Simulate == 1
        sim = sprintf(', %2.3f (Simulation)',sim_results.var_lognety);
    end
    disp(['    Var Log Net Earnings:' direct sim]);
    
    % Gini (gross earnings)
    direct = sprintf(' %5.3f (Direct)',direct_results.grossincgini);
    if p.Simulate == 1
        sim = sprintf(' %5.3f (Simulation)',sim_results.grossincgini);
    end
    disp(['    Gini coefficient, gross earnings:' direct sim])
    
    % Gini (net earnings)
    direct = sprintf(' %5.3f (Direct)',direct_results.netincgini);
    if p.Simulate == 1
        sim = sprintf(' %5.3f (Simulation)',sim_results.netincgini);
    end
    disp(['    Gini coefficient, net earnings:' direct sim])
    
    %% MPC
    fprintf('\nMPCs: \n')
    if p.freq == 1
        fprintf('  Annual (1-Period) MPCs\n')
    else
        fprintf('  Quarterly (1-Period) MPCs\n')
    end
    % Average MPC, 1 period
    for i = 1:numel(p.mpcfrac)
        if p.ComputeDirectMPC == 1
            direct = sprintf(' %2.3f (Direct)',direct_results.avg_mpc1{i});
        else
            direct = ' --- (Direct) ';
        end
        if p.Simulate == 1
            sim = sprintf(', %2.3f (Simulation)',sim_results.avg_mpc1{i});
        else
            sim = ', --- (Simulation) ';
        end
        
        label = sprintf('    MPC out of %6.2g of mean ann income:',p.mpcfrac(i));
        disp([label direct sim]);
    end
    
    % Average MPC, 4 periods
    if p.freq == 4
        fprintf('  Annual (4-Period) MPCs\n')
        for i = 1:numel(p.mpcfrac)
            if p.ComputeDirectMPC == 1
                direct = sprintf(', %2.3f (Direct)',direct_results.avg_mpc4{i});
            else
                direct = ' --- (Direct) ';
            end

            if p.Simulate == 1
                sim = sprintf(', %2.3f (Simulation)',sim_results.avg_mpc4{i});
            else
                sim = ', --- (Simulation) ';
            end

            label = sprintf('    MPC out of %6.2g of mean ann income:',p.mpcfrac(i));
            disp([label direct sim]);

        end
    end

    
    %% OTHER
    fprintf('\nOTHER: \n')
    
    % beta
    if p.freq == 4
        betamsg = sprintf(' %2.3f (Quarterly), %2.3f (Annualized),',direct_results.beta,direct_results.beta_annualized);
    elseif p.freq == 1
        betamsg = sprintf(' %2.3f (Annual),',direct_results.beta);
    end
    
    label = sprintf('    Beta =');
    if p.IterateBeta == 1
        disp([label betamsg ' found by iteration']);
    else
        disp([label betamsg ' set in parameters']);
    end
    
    % norisk MPC
    if p.freq == 1
        value = norisk_results.avg_mpc1;
        fprintf('    Mean annual MPC for model without risk (using dist of model with risk): %5.3f',value);
    else
        valueQ = norisk_results.avg_mpc1;
        valueY = norisk_results.avg_mpc4;
        fprintf('    Mean quarterly MPC for model without risk (using dist of model with risk): %5.3f\n',valueQ);
        fprintf('    Mean annual MPC for model without risk (using dist of model with risk): %5.3f\n',valueY);
    end
    
    %% ISSUES
    fprintf('\nERRORS: \n')
    for i = 1:numel(checks)
        fprintf('    %s\n',checks{i})
    end
end