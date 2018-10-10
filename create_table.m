function T = create_table(params,direct_results,norisk_results,sim_results,decomps)

    % Iterate over parameterizations
    names = {};
    tablearray = zeros(45,numel(params));
    for ip = 1:numel(params)
        p = params(ip);

        % Annual and quarterly MPCs
        if p.freq == 1
            mpcs_A = direct_results{ip}.avg_mpc1(:);
            mpcs_Q = NaN(6,1);
        else
            mpcs_A = direct_results{ip}.avg_mpc4(:);
            mpcs_Q = direct_results{ip}.avg_mpc1(:);
        end

        column = [
            direct_results{ip}.beta_annualized      % Annualized beta
            direct_results{ip}.mean_grossy_A        % Mean annual gross labor income
            direct_results{ip}.stdev_loggrossy_A    % Stdev log annual gross income
            direct_results{ip}.stdev_lognety_A      % Stdev log annual net income
            direct_results{ip}.mean_a               % Mean assets
            direct_results{ip}.constrained(1)       % Fraction with a < eps * mean ann gross inc
            direct_results{ip}.constrained(2)
            direct_results{ip}.constrained(3)
            direct_results{ip}.constrained(4)
            direct_results{ip}.constrained(5)
            direct_results{ip}.constrained(6)
            direct_results{ip}.wpercentiles(1)      % Wealth percentiles
            direct_results{ip}.wpercentiles(2)
            direct_results{ip}.wpercentiles(3)
            direct_results{ip}.wpercentiles(4)
            direct_results{ip}.wpercentiles(5)
            direct_results{ip}.wpercentiles(6)
            direct_results{ip}.wpercentiles(7)
            direct_results{ip}.top10share           % Top 10% wealth share
            direct_results{ip}.top1share            % Top 1% wealth share
            direct_results{ip}.wealthgini           % Gini coefficient
            mpcs_A(1)                               % Annual MPCs
            mpcs_A(2)
            mpcs_A(3)
            mpcs_A(4)
            mpcs_A(5)
            mpcs_A(6)
            mpcs_Q(1)                               % Quarterly MPCs (if freq = 4)
            mpcs_Q(2)
            mpcs_Q(3)
            mpcs_Q(4)
            mpcs_Q(5)
            mpcs_Q(6)
            decomps{ip}(1).term1                    % Decomposition around a=0
            decomps{ip}(1).term2  
            decomps{ip}(1).term3               
            decomps{ip}(1).term4
            decomps{ip}(2).term1                    
            decomps{ip}(2).term2
            decomps{ip}(2).term3
            decomps{ip}(2).term4
            decomps{ip}(3).term1                    
            decomps{ip}(3).term2
            decomps{ip}(3).term3
            decomps{ip}(3).term4];                
        
        % Add column name
        names{end+1} = params(ip).name;
        
        % Add this column to table
        tablearray(:,ip) = column;
    end
    
    T = table(tablearray)

end