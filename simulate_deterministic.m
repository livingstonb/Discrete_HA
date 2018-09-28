function norisk_out = simulate_deterministic(norisk,p,income,betacumdist,betacumtrans,xgrid)
    

    %% Simulate beta
    betarand = rand(p.Nsim,p.Tsim);
    betaindsim = zeros(p.Nsim,p.Tsim);
    [~,betaindsim(:,1)] = max(bsxfun(@le,betarand(:,1),betacumdist'),[],2);
    
    for it = 2:p.Tsim
        [~,betaindsim(:,it)] = max(bsxfun(@le,betarand(:,it),betacumtrans(betaindsim(:,it-1),:)),[],2);
        % ilive = diesim(:,it)==0;
        % [~,betaindsim(ilive,it)] = max(bsxfun(@le,betarand(ilive,it),betacumtrans(betaindsim(ilive,it-1),:)),[],2);
        % [~,betaindsim(~ilive,it)] = max(bsxfun(@le,betarand(~ilive,it),betacumdist'),[],2);
    end
    
 %% Simulate savings decisions
    xsim = zeros(p.Nsim,p.Tsim); 
    ssim = zeros(p.Nsim,p.Tsim);
    dierand = rand(p.Nsim,p.Tsim);
    
    diesim = dierand<p.dieprob;
    
    for it = 1:p.Tsim
        if mod(it,50) == 0
            fprintf(' Simulating norisk model, time period %3.0u \n',it);
        end

        % update cash-on-hand
        if it > 1
            xsim(:,it) = p.R * ssim(:,it-1) + income.meannety;
        end
        for ib = 1:p.nb
            idx = betaindsim(:,it) == ib;
            ssim(idx,it) = norisk.savinterp{ib}(xsim(idx,it));
        end
        ssim(ssim(:,it)<p.borrow_lim,it) = p.borrow_lim;
        csim(:,it) = xsim(:,it) - ssim(:,it) - p.savtax * max(ssim(:,it) - p.savtaxthresh,0);
    
        if p.WealthInherited == 0
            % set saving equal to 0 if hh dies at end of this period. In it+1,
            % household will have x = net income
            ssim(diesim(:,it)==1,it) = 0;
        end
    end
    
    %% MPCs
    Risk = 0;
    norisk_out = norisk;
    [norisk_out.avg_mpc1,norisk_out.avg_mpc4] = simulation_MPCs(p,xsim,csim,[],[],[],...
        betaindsim,income,Risk,norisk,xgrid);
    