function decomp2 = decomposition2(params,direct_results)
    % Based on direct_results, this function decomposes the difference in
    % expected MPC from the baseline, Em1 - Em0
    
    Nparams = numel(params);
    decomp2 = cell(1,Nparams);
    
    for ip = 1:Nparams
        decomp2{ip} = struct([]);
        
        % Loop over abar thresholds (for a <= abar, a > bar decomp)
        for ia = 1:numel(params(ip).abars)
            if params(ip).freq == 1
                baseind = 1;
            elseif params(ip).freq == 4
                baseind = 2;
            end

            try
                m1 = direct_results{ip}.mpcs1_a_direct{5};
                m0 = direct_results{baseind}.mpcs1_a_direct{5};
                g1 = direct_results{ip}.agrid_dist;
                g0 = direct_results{baseind}.agrid_dist;

                decomp2{ip}(ia).Em1_less_Em0 = direct_results{ip}.avg_mpc1_agrid(5) ...
                                - direct_results{baseind}.avg_mpc1_agrid(5);
                decomp2{ip}(ia).Em1_less_Em0_check = g1'*m1 - g0'*m0;
                decomp2{ip}(ia).term1 = g0' * (m1 - m0);
                decomp2{ip}(ia).term2 = m0' * (g1 - g0);
                decomp2{ip}(ia).term3 = (m1 - m0)' * (g1 - g0);
                ME = struct();
            catch ME
            end

            if Nparams==1 || numel(fieldnames(ME))>0
                % Either not running more than one param, or error was thrown
                decomp2{ip}(ia).Em1_less_Em0 = NaN;
                decomp2{ip}(ia).Em1_less_Em0_check = NaN;
                decomp2{ip}(ia).term1 = NaN;
                decomp2{ip}(ia).term2 = NaN;
                decomp2{ip}(ia).term3 = NaN;
                clear ME
            end
        end
    end
end