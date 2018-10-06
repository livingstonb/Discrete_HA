function [mpc1,mpc4] = direct_MPCs_deterministic(p,prefs,income,norisk)

    if p.Display == 1
        disp(' Simulating 4 periods to get deterministic MPCs')
    end
    
    % Only need to compute MPCs at one point in asset space
    % Number of draws from distribution
    Nsim = 1e5;
    betarand = rand(Nsim,4);
    dierand  = rand(Nsim,4);
    diesim   = dierand < p.dieprob;
    
    % Starting point in asset space
    x0 = income.meannety;
    
    % Consumption at starting point
    c0 = income.meannety;

    % 1 percent of mean gross annual income
    mpcamount = 0.01 * income.meany * p.freq;

    xsim        = zeros(Nsim,4);
    betaindsim  = zeros(Nsim,4);
    diesim      = zeros(Nsim,4);
    ssim        = zeros(Nsim,4);
    asim        = zeros(Nsim,4);

    %% Simulate beta
    for it = 1:4
        if it == 1
            [~,betaindsim(:,it)] = max(bsxfun(@le,betarand(:,it),prefs.betacumdist'),[],2);
        else
            [~,betaindsim(:,it)] = max(bsxfun(@le,betarand(:,it),prefs.betacumtrans(betaindsim(:,it-1),:)),[],2);
        end
    end

    %% Simulate decision variables
    for it = 1:4
        if it == 1
            % Choose lowest point in asset space, x = mean net income
            xsim(:,it) = x0 * ones(Nsim,1) + mpcamount;
        else
            xsim(:,it) = asim(:,it-1) + income.meannety;
        end

        for ib = 1:p.nb
            idx = betaindsim(:,it)==ib;
            ssim(idx,it) = norisk.savinterp{ib}(xsim(idx,it));

            asim(:,it) = p.R * ssim(:,it);
            if p.WealthInherited == 1
                asim(diesim(:,it)==1,it) = 0;
            end
        end
        ssim(ssim(:,it)<p.borrow_lim,it) = p.borrow_lim;
    end
    csim = xsim - ssim - p.savtax * max(ssim-p.savtaxthresh,0);

    %% COMPUTE MPCs
    mpc1 = mean((csim(:,1) - c0)/mpcamount);
    mpc2 = mpc1 + mean((csim(:,2) - c0)/mpcamount);
    mpc3 = mpc2 + mean((csim(:,3) - c0)/mpcamount);
    mpc4 = mpc3 + mean((csim(:,4) - c0)/mpcamount);
end