function [simulations ssim] = simulate(p,income,labtaxthresh,sav,...
    xgrid,lumptransfer,betacumdist,betacumtrans)
    
    
    grossysim = zeros(p.Nsim,p.Tsim);
    ynetsim = zeros(p.Nsim,p.Tsim);
    
    %% Simulate income process
    disp(['Simulating income process...']);
    yTrand = rand(p.Nsim,p.Tsim);
    yPrand = rand(p.Nsim,p.Tsim);
    yFrand = rand(p.Nsim,1);
    dierand = rand(p.Nsim,p.Tsim);
    
    diesim = dierand<p.dieprob;
    
    yTindsim = zeros(p.Nsim,p.Tsim);
    yPindsim = zeros(p.Nsim,p.Tsim);

    [~,yFindsim] = max(bsxfun(@le,yFrand,income.yFcumdist'),[],2);
    yFindsim = repmat(yFindsim,1,p.Tsim);
    
    % simulate yP upon death outside of time loop for speed
    idx = zeros(p.Nsim,p.Tsim);
    for iyP = 1:p.nyP
        if iyP == 1
            idx = yPrand<income.yPcumdist(iyP);
        else
            idx = yPrand<income.yPcumdist(iyP) & yPrand>=income.yPcumdist(iyP-1);
        end
        yPindsim(diesim & idx) = iyP;
    end
        
    % iterate over time periods
    for it = 1:p.Tsim
        [~,yTindsim(:,it)] = max(bsxfun(@le,yTrand(:,it),income.yTcumdist'),[],2);
        
        if it ==1
            [~,yPindsim(diesim(:,it)==0,it)] = max(bsxfun(@le,yPrand(diesim(:,it)==0,it),income.yPcumdist'),[],2);
        else
            [~,yPindsim(diesim(:,it)==0,it)] = max(bsxfun(@le,yPrand(diesim(:,it)==0,it),income.yPcumtrans(yPindsim(diesim(:,it)==0,it-1),:)),[],2);
        end
      
    end
    
    % gross income
    ygrosssim = bsxfun(@times,income.yPgrid(yPindsim).*income.yTgrid(yTindsim),income.yFgrid(yFindsim));
    if p.NormalizeY == 1
        ygrosssim = ygrosssim/income.original_meany;
    end
    
    % net income
    ynetsim = lumptransfer + (1-p.labtaxlow)*ygrosssim - p.labtaxhigh*max(ygrosssim-labtaxthresh,0);
    
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
    
    for iyF = 1:p.nyF
    for ib = 1:p.nb
    for iyP = 1:p.nyP
        savinterp{iyP,ib,iyF} = griddedInterpolant(xgrid.orig_wide(:,iyP,iyF,ib),sav.orig_wide(:,iyP,iyF,ib),'linear');
    end
    end
    end
    
    for it = 1:p.Tsim
        if mod(it,50) == 0
            fprintf(' Simulating, time period %3.0u \n',it);
        end
        % update cash-on-hand
        if it > 1
            xsim(:,it) = p.R * ssim(:,it-1) + ynetsim(:,it);
        end
        
        for iyF = 1:p.nyF
        for ib = 1:p.nb
        for iyP = 1:p.nyP
            idx = yPindsim(:,it)==iyP & betaindsim(:,it)==ib & yFindsim(:,it)==iyF;
            ssim(idx,it) = savinterp{iyP,ib,iyF}(xsim(idx,it));
        end
        end
        end
        
        ssim(ssim(:,it)<p.borrow_lim,it) = p.borrow_lim;
    end
    
    %% Moments
    simulations.mean_s = mean(ssim(:,p.Tsim));
    simulations.mean_x = mean(xsim(:,p.Tsim));
    simulations.frac_constrained = mean(ssim(:,p.Tsim)<=p.borrow_lim);
    simulations.mean_grossy = mean(ygrosssim(:,p.Tsim));
    simulations.mean_nety = mean(ynetsim(:,p.Tsim));
    simulations.mean_loggrossy = mean(log(ygrosssim(:,p.Tsim)));
    simulations.mean_lognety = mean(log(ynetsim(:,p.Tsim)));
    simulations.var_loggrossy = var(log(ygrosssim(:,p.Tsim)));
    simulations.var_lognety = var(log(ynetsim(:,p.Tsim)));
    simulations.frac_less5perc_labincome = mean(ssim(:,p.Tsim)<0.05);
    
    simulations.p10wealth = quantile(ssim(:,p.Tsim),0.1);
    simulations.p25wealth = quantile(ssim(:,p.Tsim),0.25);
    simulations.p50wealth = quantile(ssim(:,p.Tsim),0.5);
    simulations.p90wealth = quantile(ssim(:,p.Tsim),0.9);
    simulations.p99wealth = quantile(ssim(:,p.Tsim),0.99);

end