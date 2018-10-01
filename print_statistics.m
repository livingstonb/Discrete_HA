function print_statistics(results,simulations,p)
    % This function prints the main results associated with the
    % parameterization structure 'p'
    
    if p.Simulate == 0
        sim = [];
    end
    
    %% WEALTH DISTRIBUTION
    
    fprintf('\nWEALTH DISTRIBUTION: \n')
    
    % Mean wealth
    direct  = sprintf(' %2.3f (Direct Comp),',results.mean_s);
    if p.Simulate == 1
        sim = sprintf(' %2.3f (Simulation)',simulations.mean_s);
    end
    disp(['    Mean wealth:' direct sim]);
    
    % Borrowing constrained    
    for i = 1:numel(p.epsilon)
        label  = sprintf('    Fraction with a <= %3.1f%% mean ann gross inc:',p.epsilon(i)*100);
        direct = sprintf(' %3.3f (Direct)',results.constrained(i));
        if p.Simulate == 1
            sim    = sprintf(', %3.3f (Simulation)',simulations.constrained(i));
        end
        disp([label direct sim])
    end
    
    % Percentiles
    direct10 = sprintf(' %2.3f (Direct Comp)',results.p10wealth);
    direct25 = sprintf(' %2.3f (Direct Comp)',results.p25wealth);
    direct50 = sprintf(' %2.3f (Direct Comp)',results.p50wealth);
    direct90 = sprintf(' %2.3f (Direct Comp)',results.p90wealth);
    direct99 = sprintf(' %2.3f (Direct Comp)',results.p99wealth);
    if p.Simulate == 1
        sim10 = sprintf(', %2.3f (Simulation)',simulations.p10wealth);
        sim25 = sprintf(', %2.3f (Simulation)',simulations.p25wealth);
        sim50 = sprintf(', %2.3f (Simulation)',simulations.p50wealth);
        sim90 = sprintf(', %2.3f (Simulation)',simulations.p90wealth);
        sim99 = sprintf(', %2.3f (Simulation)',simulations.p99wealth);
    else
        sim10 =[]; sim25 =[]; sim50 =[]; sim90 =[]; sim99 =[];
    end
    disp(['    10th percentile:' direct10 sim10]);
    disp(['    25th percentile:' direct25 sim25]);
    disp(['    50th percentile:' direct50 sim50]);    
    disp(['    90th percentile:' direct90 sim90]);
    disp(['    99th percentile:' direct99 sim99]);

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
    direct = sprintf(' %2.3f (Direct Comp)',results.var_lognety);
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
        if p.ComputeSimMPC == 1
            sim = sprintf(', %2.3f (Simulation)',simulations.avg_mpc1{i});
        else
            sim = ', --- (Simulation) ';
        end
        
        label = sprintf('    MPC out of %4.2g of mean income:',results.mpcamount{i});
        disp([label direct sim]);
    end
    
    % Average MPC, 4 periods (simulation)
    fprintf('  4-Period MPCs\n')
    for i = 1:size(results.mpcamount,2)
        direct = ' --- (Direct Comp) ';
        if p.ComputeSimMPC == 1
            sim = sprintf(', %2.3f (Simulation)',simulations.avg_mpc4{i});
        else
            sim = ', --- (Simulation) ';
        end
        
        label = sprintf('    MPC out of %4.2g of mean income:',results.mpcamount{i});
        disp([label direct sim]);
        
    end
    
    % Average MPC, 4 periods (direct)
    if p.ComputeDirectMPC == 1
        msg = sprintf('    MPC out of %4.2g:',results.distMPCamount');
        direct = sprintf(' %4.3f (Direct Comp)',results.avg_mpc4);
        disp([msg direct]);
    end
    
    %% OTHER
    fprintf('\nOTHER: \n')
    quarterb = sprintf(' %2.3f (Quarterly),',results.betaquarterly);
    annualb = sprintf(' %2.3f (Annual),',results.betaannual);
    label = sprintf('    Beta =');
    if p.IterateBeta
        disp([label quarterb annualb ' found by iteration']);
    else
        disp([label quarterb annualb ' set in parameters']);
    end
    
    %% ISSUES
    fprintf('\nISSUES: \n')
    for i = 1:numel(results.issues)
        fprintf('    %s\n',results.issues{i})
    end
end