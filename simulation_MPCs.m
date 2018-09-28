function [avg_mpc1,avg_mpc4] = simulation_MPCs(p,xsim,csim,diesim,ynetsim,yPindsim,yFindsim,...
    betaindsim,income,simulationstruct,xgrid)
    
    %% MPCs
    Nmpcamount = numel(p.mpcfrac);
    for im = 1:Nmpcamount
        mpcamount{im} = p.mpcfrac{im} * income.meany;
        % period Tsim cash-on-hand
        xsim_mpc{im} = xsim(:,p.Tsim) + mpcamount{im};
        
        csim_mpc{im} = zeros(p.Nsim,1);
        
        xsim_mpc{im} = zeros(p.Nsim,4);
        ssim_mpc{im} = zeros(p.Nsim,4);
        xsim_mpc{im}(:,1) = xsim(:,p.Tsim-3) + mpcamount{im};
        csim_mpc{im} = zeros(p.Nsim,4);
        for it = 1:4
            if it > 1
                xsim_mpc{im}(:,it) = p.R * ssim_mpc{im}(:,it-1) + ynetsim(:,p.Tsim-4+it);
            end
            
            % if negative shock, deal with households at bottom of their cash
            % grid by setting MPC to one and cash to bottom of grid
            if mpcamount{im} < 0
                set_mpc_one = false(p.Nsim,4);
                    for ib = 1:p.nb
                    for iyF = 1:p.nyF
                    for iyP = 1:p.nyP
                        idx = yPindsim(:,it)==iyP & yFindsim(:,it)==iyF & betaindsim(:,it)==ib;
                        idx = idx & xsim_mpc{im}(:,it)<=xgrid.orig_wide(1,iyP,iyF,ib);
                        xsim_mpc{im}(idx,it) = xgrid.orig_wide(1,iyP,iyF,ib);
                        
                        % create mask for housholds whose MPC must be set
                        % to one
                        set_mpc_one(:,it) = set_mpc_one(:,it) | idx;
                    end
                    end
                    end
            end
            
            for iyF = 1:p.nyF
            for ib = 1:p.nb
            for iyP = 1:p.nyP
                idx = yPindsim(:,p.Tsim-4+it)==iyP & betaindsim(:,p.Tsim-4+it)==ib & yFindsim(:,p.Tsim-4+it)==iyF;
                ssim_mpc{im}(idx,it) = simulationstruct.savinterp{iyP,iyF,ib}(xsim_mpc{im}(idx,it));
            end
            end
            end
            ssim_mpc{im}(ssim_mpc{im}(:,it)<p.borrow_lim,it) = p.borrow_lim;
            csim_mpc{im}(:,it) = xsim_mpc{im}(:,it) - ssim_mpc{im}(:,it) - p.savtax*max(ssim_mpc{im}(:,it)-p.savtaxthresh,0);
            if p.WealthInherited == 0
                % set saving equal to 0 if hh dies at end of this period. In it+1,
                % household will have x = net income
                ssim_mpc{im}(diesim(:,p.Tsim-4+it)==1,it) = 0;
            end
        end
        
        mpc1{im} = (csim_mpc{im}(:,1) - csim(:,p.Tsim-3))/mpcamount{im};
        mpc2{im} = (csim_mpc{im}(:,2) - csim(:,p.Tsim-2))/mpcamount{im};
        mpc3{im} = (csim_mpc{im}(:,3) - csim(:,p.Tsim-1))/mpcamount{im};
        mpc4{im} = (csim_mpc{im}(:,4) - csim(:,p.Tsim))/mpcamount{im};
        
        
        % set MPCs equal to one for households at the bottom of their
        % cash-on-hand grids
        if mpcamount{im} < 0
            mpc1{im}(set_mpc_one(:,1)) = 1;
            mpc2{im}(set_mpc_one(:,2)) = 1;
            mpc3{im}(set_mpc_one(:,3)) = 1;
            mpc4{im}(set_mpc_one(:,4)) = 1;
        end
        avg_mpc1{im} = mean(mpc1{im});
        avg_mpc4{im} = mean(mpc1{im}+mpc2{im}+mpc3{im}+mpc4{im});
    end
end