function [sim_results,assetmeans] = simulate(p,income,model,...
                                        xgrid,prefs)
    % This function runs simulations based on the paratmers in 'p' and the
    % policy functions in 'model'. Output is a structure.
    
    %% Simulate income process
    disp(['Simulating income process...']);
    if p.yTContinuous == 1
        yTrand = randn(p.Nsim,p.Tsim);
    else
        yTrand = rand(p.Nsim,p.Tsim);
    end
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
    
    % simulate yT outside of time loop
    if p.yTContinuous == 1
        logyTsim = - 0.5*p.sd_logyT.^2 + yTrand*p.sd_logyT;
    else
        for iyT = 1:p.nyT
            if iyT == 1
                idx = yTrand<income.yTcumdist(iyT);
            else
                idx = yTrand<income.yTcumdist(iyT) & yTrand>=income.yTcumdist(iyT-1);
            end
            yTindsim(idx) = iyT;
        end
    end
        
    % iterate over time periods
    for it = 1:p.Tsim
        if it ==1
            [~,yPindsim(diesim(:,it)==0,it)] = max(bsxfun(@le,yPrand(diesim(:,it)==0,it),income.yPcumdist'),[],2);
        else
            [~,yPindsim(diesim(:,it)==0,it)] = max(bsxfun(@le,yPrand(diesim(:,it)==0,it),income.yPcumtrans(yPindsim(diesim(:,it)==0,it-1),:)),[],2);
        end
    end
    
    % gross income
    if p.yTContinuous == 1
        ygrosssim = bsxfun(@times,income.yPgrid(yPindsim).*exp(logyTsim),income.yFgrid(yFindsim));
    else
        ygrosssim = bsxfun(@times,income.yPgrid(yPindsim).*income.yTgrid(yTindsim),income.yFgrid(yFindsim));
    end
    
    if p.NormalizeY == 1
        ygrosssim = ygrosssim / (income.original_meany * p.freq);
    end
    
    % net income
    ynetsim = income.lumptransfer + (1-p.labtaxlow)*ygrosssim - p.labtaxhigh*max(ygrosssim-income.labtaxthresh,0);
    
    %% Simulate beta
    betarand = rand(p.Nsim,p.Tsim);
    betaindsim = zeros(p.Nsim,p.Tsim);
    [~,betaindsim(:,1)] = max(bsxfun(@le,betarand(:,1),prefs.betacumdist'),[],2);
    
    for it = 2:p.Tsim
        [~,betaindsim(:,it)] = max(bsxfun(@le,betarand(:,it),prefs.betacumtrans(betaindsim(:,it-1),:)),[],2);
        % ilive = diesim(:,it)==0;
        % [~,betaindsim(ilive,it)] = max(bsxfun(@le,betarand(ilive,it),prefs.betacumtrans(betaindsim(ilive,it-1),:)),[],2);
        % [~,betaindsim(~ilive,it)] = max(bsxfun(@le,betarand(~ilive,it),prefs.betacumdist'),[],2);
    end
    
    %% Simulate savings decisions
    xsim = zeros(p.Nsim,p.Tsim); 
    ssim = zeros(p.Nsim,p.Tsim);
    asim = zeros(p.Nsim,p.Tsim); 
    
    for it = 1:p.Tsim
        if mod(it,50) == 0
            fprintf(' Simulating, time period %3.0u \n',it);
        end
        % update cash-on-hand
        if it > 1
            xsim(:,it) = asim(:,it) + ynetsim(:,it);
        end
        
        for iyF = 1:p.nyF
        for ib = 1:p.nb
        for iyP = 1:p.nyP
            idx = yPindsim(:,it)==iyP & betaindsim(:,it)==ib & yFindsim(:,it)==iyF;
            ssim(idx,it) = model.savinterp{iyP,iyF,ib}(xsim(idx,it));
        end
        end
        end
        
        ssim(ssim(:,it)<p.borrow_lim,it) = p.borrow_lim;
        csim(:,it) = xsim(:,it) - ssim(:,it) - p.savtax * max(ssim(:,it) - p.savtaxthresh,0);
        
        if it < p.Tsim
            asim(:,it+1) = p.R * ssim(:,it);
            if p.WealthInherited == 0
                asim(diesim(:,it+1)==1,it+1) = 0;
            end
        end
    end
    
    %% Moments/important quantities
    sim_results.mean_s          = mean(ssim(:,p.Tsim));
    sim_results.mean_a          = mean(asim(:,p.Tsim));
    sim_results.mean_x          = mean(xsim(:,p.Tsim));
    sim_results.mean_grossy     = mean(ygrosssim(:,p.Tsim));
    sim_results.mean_loggrossy  = mean(log(ygrosssim(:,p.Tsim)));
    sim_results.mean_nety       = mean(ynetsim(:,p.Tsim));
    sim_results.mean_lognety    = mean(log(ynetsim(:,p.Tsim)));
    sim_results.var_loggrossy   = var(log(ygrosssim(:,p.Tsim)));
    sim_results.var_lognety     = var(log(ynetsim(:,p.Tsim)));
    sim_results.wealthgini      = ginicoeff(asim(:,p.Tsim));
    sim_results.grossincgini    = ginicoeff(ygrosssim(:,p.Tsim));
    sim_results.netincgini      = ginicoeff(ynetsim(:,p.Tsim));
    assetmeans = p.R * mean(ssim);
    
    % fraction constrained
    for i = 1:numel(p.epsilon)
        sim_results.constrained(i) = mean(ssim(:,p.Tsim)<=p.borrow_lim+p.epsilon(i)*income.meany*p.freq);
    end
    
    % wealth percentiles
    for i = 1:numel(p.percentiles)
        sim_results.wpercentiles(i) = quantile(asim(:,p.Tsim),p.percentiles(i)/100);
    end
    
    % top shares
    top10w      = quantile(asim(:,p.Tsim),0.9);
    top1w       = quantile(asim(:,p.Tsim),0.99);
    idxtop10    = asim(:,p.Tsim) > top10w;
    idxtop1     = asim(:,p.Tsim) > top1w;
    sim_results.top10share = sum(asim(idxtop10,p.Tsim))/sum(asim(:,p.Tsim));
    sim_results.top1share  = sum(asim(idxtop1,p.Tsim))/sum(asim(:,p.Tsim));
    
    %% MPCs
    
    [sim_results.avg_mpc1,sim_results.avg_mpc4,sim_results.var_mpc1,sim_results.var_mpc4]...
        = simulation_MPCs(p,asim,xsim,csim,diesim,ynetsim,yPindsim,yFindsim,...
                                            betaindsim,income,model,xgrid);


end