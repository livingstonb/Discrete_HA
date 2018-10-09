function [asim1,mpcs1,mpcs4] = direct_MPCs_deterministic(p,prefs,income,norisk,basemodel,xgrid)

    if p.Display == 1
        disp(' Simulating 4 periods to get deterministic MPCs')
    end
    
    % Draw from same stationary distribution as in stochastic case,
    % collapsed over (yP,yF,beta). This is so we can integrate with respect
    % to the probability measure associated with the model with income
    % risk.
    cumdistr = cumsum(basemodel.SSdist_noincrisk(:));
    
    % Number of draws from stationary distribution
    Nsim = 1e5;
    
    % Initial state random draw
    state_rand0 = rand(Nsim,1);
    
    betarand = rand(Nsim,4);
    dierand  = rand(Nsim,4);
    diesim   = dierand < p.dieprob;
    betaindsim  = zeros(Nsim,4);
    diesim      = zeros(Nsim,4);
    
    % Extended xgrid to coincide with above distribution
    xgrid_extended = repmat(xgrid.longgrid,p.nb,1);
    
    % Index of beta for each point in distr
    betainds = kron((1:p.nb)',ones(p.nxlong,1));
    
    % Starting point in asset space
    [~,ind0] = max(bsxfun(@le,state_rand0,cumdistr'),[],2);
    x0 = xgrid_extended(ind0);
    betaindsim0 = betainds(ind0)';
    
    % 1 percent of mean gross annual income
    mpcamount = 0.01 * income.meany * p.freq;


    %% Simulate beta
    for it = 1:4
        if it == 1
            [~,betaindsim(:,it)] = max(bsxfun(@le,betarand(:,it),prefs.betacumtrans(betaindsim0,:)),[],2);
        else
            [~,betaindsim(:,it)] = max(bsxfun(@le,betarand(:,it),prefs.betacumtrans(betaindsim(:,it-1),:)),[],2);
        end
    end
    
    %% Simulate 0th period to get starting assets and cash-in-hand
    sav0 = zeros(Nsim,1);
    for ib = 1:p.nb
        idx = betaindsim0==ib;
        sav0(idx) = basemodel.savinterp{ib}(x0(idx));
    end
    sav0 = max(sav0,p.borrow_lim);
    asim1 = p.R * sav0;
    if p.WealthInherited == 0
        asim1(diesim(:,1)==1) = 0;
    end
    

    %% Simulate decision variables
    
    % First, get the baseline
    for shock = [0,mpcamount]
        xsim = zeros(Nsim,4);
        ssim = zeros(Nsim,4);
        asim = zeros(Nsim,4);
    
        for it = 1:4
            if it == 1
                xsim(:,it) = asim1 + income.meannety + shock;
            else
                xsim(:,it) = asim(:,it) + income.meannety;
            end

            for ib = 1:p.nb
                idx = betaindsim(:,it)==ib;
                ssim(idx,it) = norisk.savinterp{ib}(xsim(idx,it));
            end
            ssim(:,it) = max(ssim(:,it),p.borrow_lim);

            % Assets
            if it < 4
                asim(:,it+1) = p.R * ssim(:,it);
                if p.WealthInherited == 1
                    asim(diesim(:,it+1)==1,it+1) = 0;
                end
            end

        end
        
        if shock == 0
            csim_noshock = xsim - ssim - p.savtax * max(ssim-p.savtaxthresh,0);
        else
            csim = xsim - ssim - p.savtax * max(ssim-p.savtaxthresh,0);
        end
   
    end
    
     %% COMPUTE MPCs
    mpcs1 = (csim(:,1) - csim_noshock(:,1))/mpcamount;
    mpcs2 = mpcs1 + (csim(:,2) - csim_noshock(:,2))/mpcamount;
    mpcs3 = mpcs2 + (csim(:,3) - csim_noshock(:,3))/mpcamount;
    mpcs4 = mpcs3 + (csim(:,4) - csim_noshock(:,4))/mpcamount;
end