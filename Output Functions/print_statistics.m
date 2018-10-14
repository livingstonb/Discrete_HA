function print_statistics(direct_results,sim_results,norisk_results,checks,p,decomp)
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
        label  = sprintf('    Fraction with a <= %3.1f%% mean ann gross inc:',p.epsilon(i)*100);
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
    % Stdev of log gross labor income
    direct = sprintf(' %2.3f (Direct)',direct_results.stdev_loggrossy_A);
    if p.Simulate == 1
        sim = sprintf(', %2.3f (Simulation)',sim_results.var_loggrossy);
    end
    disp(['    StDev Log Gross Earnings:' direct sim]);
    
    % Stdev of log net labor income
    direct = sprintf('   %2.3f (direct)',direct_results.stdev_lognety_A);
    if p.Simulate == 1
        sim = sprintf(', %2.3f (Simulation)',sim_results.var_lognety);
    end
    disp(['    StDev Log Net Earnings:' direct sim]);
    
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
            direct = sprintf(' %2.3f (Direct)',direct_results.avg_mpc1_sim(i));
        else
            direct = ' --- (Direct) ';
        end
        if p.Simulate == 1
            sim = sprintf(', %2.3f (Simulation)',sim_results.avg_mpc1(i));
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
                direct = sprintf(', %2.3f (Direct)',direct_results.avg_mpc4_sim(i));
            else
                direct = ' --- (Direct) ';
            end

            if p.Simulate == 1
                sim = sprintf(', %2.3f (Simulation)',sim_results.avg_mpc4(i));
            else
                sim = ', --- (Simulation) ';
            end

            label = sprintf('    MPC out of %6.2g of mean ann income:',p.mpcfrac(i));
            disp([label direct sim]);

        end
    end
    
    %% DECOMPOSITION
    fprintf('\nDECOMPOSITION: \n')
    disp(' For a == 0')
    fprintf('    RA MPC:                                               %5.3f\n',decomp(1).term1)
    fprintf('    Effect of HtM households:                             %5.3f\n',decomp(1).term2)
    fprintf('    Effect of borrowing constraint on non-HtM households: %5.3f\n',decomp(1).term3)
    fprintf('    Effect of idiosyncratic risk on non-HtM households:   %5.3f\n',decomp(1).term4)
    
    disp(' For a <=0.01')
    fprintf('    RA MPC:                                               %5.3f\n',decomp(2).term1)
    fprintf('    Effect of HtM households:                             %5.3f\n',decomp(2).term2)
    fprintf('    Effect of borrowing constraint on non-HtM households: %5.3f\n',decomp(2).term3)
    fprintf('    Effect of idiosyncratic risk on non-HtM households:   %5.3f\n',decomp(2).term4)
    
    disp(' For a <= 0.05')
    fprintf('    RA MPC:                                               %5.3f\n',decomp(3).term1)
    fprintf('    Effect of HtM households:                             %5.3f\n',decomp(3).term2)
    fprintf('    Effect of borrowing constraint on non-HtM households: %5.3f\n',decomp(3).term3)
    fprintf('    Effect of idiosyncratic risk on non-HtM households:   %5.3f\n',decomp(3).term4)
    
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

    
    %% ISSUES
    fprintf('\nERRORS: \n')
    for i = 1:numel(checks)
        fprintf('    %s\n',checks{i})
    end
end