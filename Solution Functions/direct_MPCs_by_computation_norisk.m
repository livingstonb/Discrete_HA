function norisk_mpcs1_a_direct = direct_MPCs_by_computation_norisk(p,norisk,income,prefs,grids)

 %% DIRECTLY COMPUTED 1-PERIOD MPCs (MODEL WITHOUT INCOME RISK)
    norisk_mpcs1_a_direct = cell(1,6);
 
    for im = 0:numel(p.mpcfrac)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.mpcfrac(im)*income.meany1*p.freq;
        end
        
        x_mpc = grids.a.vec + income.meany1 + mpcamount;
        con = zeros(p.nx_KFE,p.nb);
        for ib = 1:p.nb
            con(:,ib) = norisk.coninterp{ib}(x_mpc);
        end
        
        if im == 0
            con_baseline = con;
        else
            % Compute m(a,beta)
            mpcs1_a_beta = (con - con_baseline) / mpcamount;

            % Compute m(x) = E(m(x,beta)|x)
            %       = sum of P(beta|x) * m(x,beta) over all beta
            % beta is exogenous so P(beta|x) = P(beta)
            norisk_mpcs1_a_direct{im} = mpcs1_a_beta * prefs.betadist;
        end
    end

end