function MPCs = simulation_MPCs(p,xsim,csim,diesim,ynetsim,yPindsim,yFindsim,...
                                                                    betaindsim,income,basemodel,xgrid)
    % This function is called by simulate.m to compute MPCs. Outputs are
    % cell arrays, each cell associated with one mpcamount.
    
    Tmax  = p.freq * 4;
    ynetsim_mpc = ynetsim(:,end-Tmax+1:end);
    diesim_mpc = diesim(:,end-Tmax+1:end);
    csim_baseline = csim(:,end-Tmax+1:end);
    yPindsim_mpc = yPindsim(:,end-Tmax+1:end);
    betaindsim_mpc = betaindsim(:,end-Tmax+1:end);
    
    %% MPCs
    Nmpcamount = numel(p.mpcfrac);
    for im = 1:Nmpcamount
        mpcamount{im} = p.mpcfrac(im) * income.meany1 * p.freq;
        
        if mpcamount{im} < 0
            set_mpc_one = false(p.Nsim,1);
        end

        xsim_mpc        = zeros(p.Nsim,4);
        ssim_mpc        = zeros(p.Nsim,4);
        asim_mpc        = zeros(p.Nsim,4);
        xsim_mpc(:,1)   = xsim(:,p.Tsim-Tmax+1) + mpcamount{im};
        csim_mpc        = zeros(p.Nsim,4);
        
        for it = 1:Tmax
            if it > 1
                xsim_mpc(:,it) = asim_mpc(:,it) + ynetsim_mpc(:,it);
            end

            for iyF = 1:p.nyF
            for ib = 1:p.nb
            for iyP = 1:p.nyP
                below_grid = xsim_mpc(:,it)<xgrid.longgrid(1,iyP,iyF);
                idx = yPindsim_mpc(:,it)==iyP & betaindsim_mpc(:,it)==ib & yFindsim(:)==iyF;
                % if shock is negative, deal with households that wind up
                % below the bottom of their asset grid
                if mpcamount{im} < 0 && it == 1
                    idx_below = idx & below_grid;
                    xsim_mpc(idx_below,it) = xgrid.longgrid(1,iyP,iyF);
                    % will need to set their MPCs to one
                    set_mpc_one = set_mpc_one | idx_below;
                end
                ssim_mpc(idx,it) = basemodel.savinterp{iyP,iyF,ib}(xsim_mpc(idx,it));
            end
            end
            end
            ssim_mpc(ssim_mpc(:,it)<p.borrow_lim,it) = p.borrow_lim;
            csim_mpc(:,it) = xsim_mpc(:,it) - ssim_mpc(:,it) - p.savtax*max(ssim_mpc(:,it)-p.savtaxthresh,0);
            
            if it < Tmax
                asim_mpc(:,it+1) = p.R * ssim_mpc(:,it);
                if p.Bequests == 0
                    % set assets equal to 0 if hh dies at end of this period
                    asim_mpc(diesim_mpc(:,it+1)==1,it+1) = 0;
                end
            end
        end
        
        mpcs_1_1 = (csim_mpc(:,1) - csim_baseline(:,1))/mpcamount{im};
        mpcs_1_2 = (csim_mpc(:,2) - csim_baseline(:,2))/mpcamount{im};
        mpcs_1_3 = (csim_mpc(:,3) - csim_baseline(:,3))/mpcamount{im};
        mpcs_1_4 = (csim_mpc(:,4) - csim_baseline(:,4))/mpcamount{im};
        
        
        % set MPCs equal to one for households at the bottom of their
        % cash-on-hand grids
        if mpcamount{im} < 0
            mpcs_1_1(set_mpc_one) = 1;
            mpcs_1_2(set_mpc_one) = 1;
            mpcs_1_3(set_mpc_one) = 1;
            mpcs_1_4(set_mpc_one) = 1;
        end

        MPCs.avg_1_1(im) = mean(mpcs_1_1);
        MPCs.avg_1_2(im) = mean(mpcs_1_2);
        MPCs.avg_1_3(im) = mean(mpcs_1_3);
        MPCs.avg_1_4(im) = mean(mpcs_1_4);
        MPCs.avg_1_1to4(im) = mean(mpcs_1_1+mpcs_1_2+mpcs_1_3+mpcs_1_4);
        
        if p.freq == 4
            mpcs_1_x = cell(1,16);
            for ip = 5:16
                mpcs_1_x{ip} = (csim_mpc(:,ip) - csim_baseline(:,ip))/mpcamount{im};
            end
            
            MPCs.avg_1_5to8 = mean(mpcs_1_x{5}+mpcs_1_x{6}+mpcs_1_x{7}+mpcs_1_x{8});
            MPCs.avg_1_9to12 = mean(mpcs_1_x{9}+mpcs_1_x{10}+mpcs_1_x{11}+mpcs_1_x{12});
            MPCs.avg_1_13to16 = mean(mpcs_1_x{13}+mpcs_1_x{14}+mpcs_1_x{15}+mpcs_1_x{16});
        else
            MPCs.avg_1_5to8 = NaN;
            MPCs.avg_1_9to12 = NaN;
            MPCs.avg_1_13to16 = NaN;
        end
        
        MPCs.var_1_1(im) = var(mpcs_1_1);
        MPCs.var_1_4(im) = var(mpcs_1_4);
    end
end