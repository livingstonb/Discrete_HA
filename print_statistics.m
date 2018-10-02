function print_statistics(direct_results,sim_results,p)
    % This function prints the main results associated with the
    % parameterization structure 'p'
    
    if p.Simulate == 0
        sim = [];
    end
    
    %% WEALTH DISTRIBUTION
    
    fprintf('\nWEALTH/SAVINGS: \n')
    
    % Mean wealth
    direct  = sprintf(' %2.3f (Direct Comp),',direct_results.mean_a);
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
    
    
    %% Income Distribution
    fprintf('\nINCOME DISTRIBUTION: \n')
    fprintf('  Cross-sectional variances\n')
    % Variance of log gross labor income
    direct = sprintf(' %2.3f (Direct Comp)',direct_results.var_loggrossy);
    if p.Simulate == 1
        sim = sprintf(', %2.3f (Simulation)',sim_results.var_loggrossy);
    end
    disp(['    Var Log Gross Earnings:' direct sim]);
    
    % Variance of log net labor income
    direct = sprintf('   %2.3f (Direct Comp)',direct_results.var_lognety);
    if p.Simulate == 1
        sim = sprintf(', %2.3f (Simulation)',sim_results.var_lognety);
    end
    disp(['    Var Log Net Earnings:' direct sim]);
    
    %% MPC
    fprintf('\nMPCs: \n')
    fprintf('  1-Period MPCs\n')
    % Average MPC, 1 period
    for i = 1:size(p.mpcfrac,2)
        if p.ComputeDirectMPC == 1
            direct = sprintf(' %2.3f (Direct Comp)',direct_results.avg_mpc1{i});
        else
            direct = ' --- (Direct Comp) ';
        end
        if p.Simulate == 1
            sim = sprintf(', %2.3f (Simulation)',sim_results.avg_mpc1{i});
        else
            sim = ', --- (Simulation) ';
        end
        
        label = sprintf('    MPC out of %6.2g of mean income:',p.mpcfrac{i});
        disp([label direct sim]);
    end
    
    % Average MPC, 4 periods (simulation)
    fprintf('  4-Period MPCs\n')
    for i = 1:numel(p.mpcfrac)
        direct = ' --- (Direct Comp) ';
        if p.Simulate == 1
            sim = sprintf(', %2.3f (Simulation)',sim_results.avg_mpc4{i});
        else
            sim = ', --- (Simulation) ';
        end
        
        label = sprintf('    MPC out of %6.2g of mean income:',p.mpcfrac{i});
        disp([label direct sim]);
        
    end
    
    % Average MPC, 4 periods (direct)
    if p.ComputeDirectMPC == 1
        msg = sprintf('    MPC out of %6.2g of mean income:',direct_results.avg_mpc4.mpcfrac);
        direct = sprintf(' %4.3f (Direct Comp)',direct_results.avg_mpc4.value);
        disp([msg direct]);
    end
    
    %% OTHER
    fprintf('\nOTHER: \n')
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
    fprintf('\nISSUES: \n')
    for i = 1:numel(direct_results.issues)
        fprintf('    %s\n',direct_results.issues{i})
    end
end