function print_statistics(results,simulations,p)
    if p.Simulate == 0
        sim = [];
    end
    
    %% WEALTH DISTRIBUTION
    
    fprintf('\nWEALTH DISTRIBUTION: \n')
    
    % Mean wealth
    direct  = sprintf(' %2.3f (Direct Comp)',results.mean_s);
    if p.Simulate == 1
        sim = sprintf(' %2.3f (Simulation)',simulations.mean_s);
    end
    disp(['    Mean wealth:' direct sim]);
    
    % Borrowing constrained
    direct  = sprintf(' %2.3f (Direct Comp)',results.frac_constrained);
    if p.Simulate == 1
        sim = sprintf(' %2.3f (Simulation)',simulations.frac_constrained);
    end
    disp(['    Fraction constrained:' direct sim]);
    
    % Percent whose wealth is less than 5% labor income
    direct = sprintf(' %2.3f (Direct Comp)',results.frac_less5perc_labincome);
    if p.Simulate == 1
        sim = sprintf(' %2.3f (Simulation)',simulations.frac_less5perc_labincome);
    end
    disp(['    Percent with wealth < 5 percent labinc:' direct sim]);
    
    % Percentiles
    direct10 = sprintf(' %2.3f (Direct Comp)',results.p10wealth);
    direct25 = sprintf(' %2.3f (Direct Comp)',results.p25wealth);
    direct50 = sprintf(' %2.3f (Direct Comp)',results.p50wealth);
    direct90 = sprintf(' %2.3f (Direct Comp)',results.p90wealth);
    direct99 = sprintf(' %2.3f (Direct Comp)',results.p99wealth);
    if p.Simulate == 1
        sim10 = sprintf(' %2.3f (Simulation)',simulations.p10wealth);
        sim25 = sprintf(' %2.3f (Simulation)',simulations.p25wealth);
        sim50 = sprintf(' %2.3f (Simulation)',simulations.p50wealth);
        sim90 = sprintf(' %2.3f (Simulation)',simulations.p90wealth);
        sim99 = sprintf(' %2.3f (Simulation)',simulations.p99wealth);
    else
        sim10 =[];
        sim25 =[];
        sim50 =[];
        sim90 =[];
        sim99 =[];
    end
    disp(['    10th percentile:' direct10 sim10]);
    disp(['    25th percentile:' direct25 sim25]);
    disp(['    50th percentile:' direct50 sim50]);    
    disp(['    90th percentile:' direct90 sim90]);
    disp(['    99th percentile:' direct99 sim99]);

    %% Income Distribution
    fprintf('\nINCOME DISTRIBUTION: \n')
    % Variance of log gross labor income
    direct = sprintf(' %2.3f (Direct Comp)',results.var_loggrossy);
    if p.Simulate == 1
        sim = sprintf(' %2.3f (Direct Comp)',simulations.var_loggrossy);
    end
    disp(['    Var Log Gross Earnings:' direct sim]);
    
    % Variance of log net labor income
    direct = sprintf(' %2.3f (Direct Comp)',results.var_lognety);
    if p.Simulate == 1
        sim = sprintf(' %2.3f (Simulation)',simulations.var_lognety);
    end
    disp(['    Var Log Net Earnings:' direct sim]);
end