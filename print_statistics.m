function print_statistics(results,simulations,p)
    % This function prints the main results associated with the
    % parameterization structure 'p'
    
    if p.Simulate == 0
        sim = [];
    end
    
    %% WEALTH DISTRIBUTION
    
    fprintf('\nWEALTH/SAVINGS: \n')
    
    % Mean wealth
    direct  = sprintf(' %2.3f (Direct Comp),',results.mean_a);
    if p.Simulate == 1
        sim = sprintf(' %2.3f (Simulation)',simulations.mean_a);
    end
    disp(['    Mean wealth:' direct sim]);
    
    % Borrowing constrained    
    for i = 1:numel(p.epsilon)
        label  = sprintf('    Fraction with s <= 43.1f%% mean ann gross inc:',p.epsilon(i)*100);
        direct = sprintf(' %3.3f (Direct)',results.constrained(i));
        if p.Simulate == 1
            sim    = sprintf(', %3.3f (Simulation)',simulations.constrained(i));
        end
        disp([label direct sim])
    end
    
    % Percentiles
    for i = 1:numel(p.percentiles)
        label = sprintf('    Wealth, %ith percentile:',p.percentiles(i));
        direct = sprintf(' %8.3f (Direct)',results.wpercentiles(i));
        if p.Simulate == 1
            sim = sprintf(', %8.3f (Simulation)',simulations.wpercentiles(i));
        end
        disp([label direct sim]);
    end
    
    % Top shares
    direct10 = sprintf(' %6.3f (Direct)',results.top10share);
    direct1  = sprintf(' %6.3f (Direct)',results.top1share);
    if p.Simulate == 1
        sim10 = sprintf(', %6.3f (Simulation)',simulations.top10share);
        sim1 = sprintf(', %6.3f (Simulation)',simulations.top1share);
    else
        sim10 = []; sim1 = [];
    end
    disp(['    Wealth, top 10% share:' direct10 sim10])
    disp(['    Wealth, top  1% share:' direct1 sim1])
    
    
    %% Income Distribution
    fprintf('\nINCOME DISTRIBUTION: \n')
    fprintf('  Cross-sectional variances\n')
    % Variance of log gross labor income
    direct = sprintf(' %2.3f (Direct Comp)',results.var_loggrossy);
    if p.Simulate == 1
        sim = sprintf(', %2.3f (Simulation)',simulations.var_loggrossy);
    end
    disp(['    Var Log Gross Earnings:' direct sim]);
    
    % Variance of log net labor income
    direct = sprintf('   %2.3f (Direct Comp)',results.var_lognety);
    if p.Simulate == 1
        sim = sprintf(', %2.3f (Simulation)',simulations.var_lognety);
    end
    disp(['    Var Log Net Earnings:' direct sim]);
    
    %% MPC
    fprintf('\nMPCs: \n')
    fprintf('  1-Period MPCs\n')
    % Average MPC, 1 period
    for i = 1:size(results.mpcamount,2)
        if p.ComputeDirectMPC == 1
            direct = sprintf(' %2.3f (Direct Comp)',results.avg_mpc1{i});
        else
            direct = ' --- (Direct Comp) ';
        end
        if p.Simulate == 1
            sim = sprintf(', %2.3f (Simulation)',simulations.avg_mpc1{i});
        else
            sim = ', --- (Simulation) ';
        end
        
        label = sprintf('    MPC out of %6.2g of mean income:',results.mpcamount{i});
        disp([label direct sim]);
    end
    
    % Average MPC, 4 periods (simulation)
    fprintf('  4-Period MPCs\n')
    for i = 1:size(results.mpcamount,2)
        direct = ' --- (Direct Comp) ';
        if p.Simulate == 1
            sim = sprintf(', %2.3f (Simulation)',simulations.avg_mpc4{i});
        else
            sim = ', --- (Simulation) ';
        end
        
        label = sprintf('    MPC out of %6.2g of mean income:',results.mpcamount{i});
        disp([label direct sim]);
        
    end
    
    % Average MPC, 4 periods (direct)
    if p.ComputeDirectMPC == 1
        msg = sprintf('    MPC out of %6.2g:',results.distMPCamount');
        direct = sprintf(' %4.3f (Direct Comp)',results.avg_mpc4);
        disp([msg direct]);
    end
    
    %% OTHER
    fprintf('\nOTHER: \n')
    if p.freq == 4
        betafreq = sprintf(' %2.3f (Quarterly),',results.betaquarterly);
    elseif p.freq == 1
        betamsg = sprintf(' %2.3f (Annual),',results.betaannual);
    end
    label = sprintf('    Beta =');
    if p.IterateBeta == 1
        disp([label betamsg ' found by iteration']);
    else
        disp([label betamsg ' set in parameters']);
    end
    
    %% ISSUES
    fprintf('\nISSUES: \n')
    for i = 1:numel(results.issues)
        fprintf('    %s\n',results.issues{i})
    end
end