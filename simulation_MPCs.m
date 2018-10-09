function [avg_mpc1,avg_mpc4,var_mpc1,var_mpc4] = simulation_MPCs(p,xsim,csim,diesim,ynetsim,yPindsim,yFindsim,...
    betaindsim,income,simulationstruct,xgrid)
    % This function is called by simulate.m to compute MPCs. Outputs are
    % cell arrays, each cell associated with one mpcamount.
    
    %% MPCs
    Nmpcamount = numel(p.mpcfrac);
    for im = 1:Nmpcamount
        mpcamount{im} = p.mpcfrac(im) * income.meany * p.freq;
        
        if mpcamount{im} < 0
            set_mpc_one = false(p.Nsim,1);
        end

        xsim_mpc{im}        = zeros(p.Nsim,4);
        ssim_mpc{im}        = zeros(p.Nsim,4);
        asim_mpc{im}        = zeros(p.Nsim,4);
        xsim_mpc{im}(:,1)   = xsim(:,p.Tsim-3) + mpcamount{im};
        csim_mpc{im}        = zeros(p.Nsim,4);
        for it = 1:4
            simT = p.Tsim - 4 + it;
            if it > 1
                xsim_mpc{im}(:,it) = asim_mpc{im}(:,it) + ynetsim(:,simT);
            end

            for iyF = 1:p.nyF
            for ib = 1:p.nb
            for iyP = 1:p.nyP
                below_grid = xsim_mpc{im}(:,it)<xgrid.longgrid_wide(1,iyP,iyF);
                idx = yPindsim(:,simT)==iyP & betaindsim(:,simT)==ib & yFindsim(:,simT)==iyF;
                % if shock is negative, deal with households that wind up
                % below the bottom of their asset grid
                if mpcamount{im} < 0 && it == 1
                    idx_below = idx & below_grid;
                    xsim_mpc{im}(idx_below,it) = xgrid.longgrid_wide(1,iyP,iyF);
                    % will need to set their MPCs to one
                    set_mpc_one = set_mpc_one | idx_below;
                end
                ssim_mpc{im}(idx,it) = simulationstruct.savinterp{iyP,iyF,ib}(xsim_mpc{im}(idx,it));
            end
            end
            end
            ssim_mpc{im}(ssim_mpc{im}(:,it)<p.borrow_lim,it) = p.borrow_lim;
            csim_mpc{im}(:,it) = xsim_mpc{im}(:,it) - ssim_mpc{im}(:,it) - p.savtax*max(ssim_mpc{im}(:,it)-p.savtaxthresh,0);
            
            if it < 4
                asim_mpc{im}(:,it+1) = p.R * ssim_mpc{im}(:,it);
                if p.Bequests == 0
                    % set assets equal to 0 if hh dies at end of this period
                    asim_mpc{im}(diesim(:,simT+1)==1,it+1) = 0;
                end
            end
        end
        
        mpc1{im} = (csim_mpc{im}(:,1) - csim(:,p.Tsim-3))/mpcamount{im};
        mpc2{im} = mpc1{im} + (csim_mpc{im}(:,2) - csim(:,p.Tsim-2))/mpcamount{im};
        mpc3{im} = mpc2{im} + (csim_mpc{im}(:,3) - csim(:,p.Tsim-1))/mpcamount{im};
        mpc4{im} = mpc3{im} + (csim_mpc{im}(:,4) - csim(:,p.Tsim))/mpcamount{im};
        
        
        % set MPCs equal to one for households at the bottom of their
        % cash-on-hand grids
        if mpcamount{im} < 0
            mpc1{im}(set_mpc_one,1) = 1;
        end

        avg_mpc1{im} = mean(mpc1{im});
        avg_mpc4{im} = mean(mpc4{im});
        
        var_mpc1{im} = var(mpc1{im});
        var_mpc4{im} = var(mpc4{im});
    end
end