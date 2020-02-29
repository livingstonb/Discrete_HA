function norisk_mpcs1_a_direct = direct_MPCs_by_computation_norisk(...
    p, norisk, income, heterogeneity, grids)

 	%% DIRECTLY COMPUTED 1-PERIOD MPCs (MODEL WITHOUT INCOME RISK)
    norisk_mpcs1_a_direct = cell(1,6);
    for im = 1:6
        norisk_mpcs1_a_direct{im} = NaN;
    end
 
    for im = 0:numel(p.shocks)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.shocks(im) * income.meany1 * p.freq;
        end
        
        x_mpc = grids.a.vec + income.meannety1 + mpcamount;
        con = zeros(p.nx_DST, p.nb);
        for ib = 1:p.nb
            con(:,ib) = norisk.coninterp{ib}(x_mpc);
        end
        
        if im == 0
            con_baseline = con;
        else
            % Compute m(a,z)
            mpcs1_a_z = (con - con_baseline) / mpcamount;

            % Compute m(x) = E(m(x,z)|x)
            %       = sum of P(z|x) * m(x,z) over all z
            norisk_mpcs1_a_direct{im} = mpcs1_a_z * heterogeneity.zdist;
        end
    end

end