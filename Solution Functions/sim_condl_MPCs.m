function MPCs = simulation_MPCs(p,income,basemodel,xgrid,a)
    % This function is called to compute E[MPC|a]
    
    Tmax  = p.freq * 4;

    %% Income and patience draws
    yFrand = zeros(p.Nsim,1,'single');
    yPrand = zeros(p.Nsim,Tmax,'single');
    yTrand = zeros(p.Nsim,Tmax,'single');
    dierand = rand(p.Nsim,p.Tsim,'single');
    betarand = rand(p.Nsim,p.Tsim,'single');

    diesim = dierand < p.dieprob;

    yTindsim = zeros(p.Nsim,Tmax,'int8');
    yPindsim = zeros(p.Nsim,Tmax,'int8');

    % get cumulative distributions
    adist = reshape(results.sim.adist,[p.nxlong p.nyF p.nyP p.nb]);
    adist = squeeze(adist(adist(:,1,1,1)==a,:));
    adist = adist / sum(adist(:));
    yFcumdist = cumsum(sum(sum(adist,3),2));
    yPcumdist = cumsum(sum(sum(adist,3),1)');
    betacumdist = cumsum(squeeze(sum(sum(adist,2),1)));


    % initial indices
    [~,yFindsim] = max(yFrand < yFcumdist',[],2);
    yPindsim(:,1) = max(yFrand(:,1) < yPcumdist',[],2);
    yTindsim(:,1) = max(yTrand(:,1) < yTcumdist',[],2);
    betaindsim(:,1) = max(betarand(:,1) < betacumdist',[],2)

    for it = 1:Tmax
    	dead = diesim(:,it) == 1;
    	
    	yPindsim(dead,it) = max(yPrand(dead,it) < yPcumdist',[],2);
    	yTindsim(:,it) = max(yTrand(:,it) < yTcumdist',[],2);
    	betaindsim(:,it) = max(betarand(:,it) < betacumdist',[],2);

    	yPindsim(~dead,it) = max(yPrand(~dead,it) < yPcumtrans(yPindsim(~dead,it-1),:),[],2);
    end
    
    ygrosssim = income.yPgrid(yPindsim) .* income.yTgrid(yTindsim) .* income.yFgrid(yFindsim);

    % net income
    ynetsim = income.lumptransfer + (1-p.labtaxlow)*ygrosssim - p.labtaxhigh*max(ygrosssim-income.labtaxthresh,0);
    
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
        xsim_mpc	 	= zeros(p.Nsim,4);
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
            mpcs_1_2(set_mpc_one) = 0;
            mpcs_1_3(set_mpc_one) = 0;
            mpcs_1_4(set_mpc_one) = 0;
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
            
            MPCs.avg_1_5to8(im) = mean(mpcs_1_x{5}+mpcs_1_x{6}+mpcs_1_x{7}+mpcs_1_x{8});
            MPCs.avg_1_9to12(im) = mean(mpcs_1_x{9}+mpcs_1_x{10}+mpcs_1_x{11}+mpcs_1_x{12});
            MPCs.avg_1_13to16(im) = mean(mpcs_1_x{13}+mpcs_1_x{14}+mpcs_1_x{15}+mpcs_1_x{16});
        else
            MPCs.avg_1_5to8(im) = NaN;
            MPCs.avg_1_9to12(im) = NaN;
            MPCs.avg_1_13to16(im) = NaN;
        end
        
        MPCs.var_1_1(im) = var(mpcs_1_1);
        MPCs.var_1_4(im) = var(mpcs_1_4);
    end
end