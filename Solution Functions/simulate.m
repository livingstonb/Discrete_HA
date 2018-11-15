function [sim_results,assetmeans] = simulate(p,income,model,xgrid,prefs)

    % This function runs simulations based on the paratmers in 'p' and the
    % policy functions in 'model'.
    
    %% Simulate income process
    disp(['Simulating income process...']);
    if p.yTContinuous == 1
        yTrand = randn(p.Nsim,p.Tsim,'single');
    else
        yTrand = rand(p.Nsim,p.Tsim,'single');
    end
    yPrand = rand(p.Nsim,p.Tsim,'single');
    yFrand = rand(p.Nsim,1,'single');
    dierand = rand(p.Nsim,p.Tsim,'single');
    
    diesim = dierand < p.dieprob;
    
    yTindsim = zeros(p.Nsim,p.Tsim,'int8');
    yPindsim = zeros(p.Nsim,p.Tsim,'int8');

    [~,yFindsim] = max(bsxfun(@le,yFrand,income.yFcumdist'),[],2);
    
    % simulate yP upon death outside of time loop for speed
%     for iyP = 1:p.nyP
%         if iyP == 1
%             idx = yPrand<income.yPcumdist(iyP);
%         else
%             idx = (yPrand<income.yPcumdist(iyP)) & (yPrand>=income.yPcumdist(iyP-1));
%         end
%         yPindsim(diesim & idx) = iyP;
%     end
    
    % simulate yT outside of time loop
    if p.yTContinuous == 1 && p.nyT > 1
        lambdarand = rand(p.Nsim,p.Tsim,'single');
        logyTsim = (lambdarand < p.lambdaT) .* (- 0.5*p.sd_logyT.^2 + yTrand*p.sd_logyT);
    elseif p.yTContinuous == 0 && p.nyT > 1
        for iyT = 1:p.nyT
            if iyT == 1
                idx = yTrand<income.yTcumdist(iyT);
            else
                idx = (yTrand<income.yTcumdist(iyT)) & (yTrand>=income.yTcumdist(iyT-1));
            end
            yTindsim(idx) = iyT;
        end
    elseif p.nyT == 0
        logyTsim = 0;
    end
    
        
    % iterate over time periods
    for it = 1:p.Tsim
        [~,yPindsim(diesim(:,it)==1,it)] = max(bsxfun(@le,yPrand(diesim(:,it)==1,it),income.yPcumdist'),[],2);
        if it ==1
            [~,yPindsim(diesim(:,it)==0,it)] = max(bsxfun(@le,yPrand(diesim(:,it)==0,it),income.yPcumdist'),[],2);
        else
            [~,yPindsim(diesim(:,it)==0,it)] = max(bsxfun(@le,yPrand(diesim(:,it)==0,it),income.yPcumtrans(yPindsim(diesim(:,it)==0,it-1),:)),[],2);
        end
    end
    
    % gross income
    if (p.yTContinuous==1) || (p.nyT==0)
        ygrosssim = bsxfun(@times,income.yPgrid(yPindsim).*exp(logyTsim),income.yFgrid(yFindsim));
    else
        ygrosssim = bsxfun(@times,income.yPgrid(yPindsim).*income.yTgrid(yTindsim),income.yFgrid(yFindsim));
    end
    
    % net income
    ynetsim = income.lumptransfer + (1-p.labtaxlow)*ygrosssim - p.labtaxhigh*max(ygrosssim-income.labtaxthresh,0);
    
    %% Simulate beta
    betarand = rand(p.Nsim,p.Tsim);
    betaindsim = zeros(p.Nsim,p.Tsim,'int8');
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
    csim = zeros(p.Nsim,p.Tsim);
    
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
            idx = (yPindsim(:,it)==iyP) & (betaindsim(:,it)==ib) & (yFindsim(:)==iyF);
            ssim(idx,it) = model.savinterp{iyP,iyF,ib}(xsim(idx,it));
        end
        end
        end
        
        ssim(ssim(:,it)<p.borrow_lim,it) = p.borrow_lim;
        csim(:,it) = xsim(:,it) - ssim(:,it) - p.savtax * max(ssim(:,it) - p.savtaxthresh,0);
        
        if it < p.Tsim
            asim(:,it+1) = p.R * ssim(:,it);
            if p.Bequests == 0
                % set saving equal to 0 if hh dies at end of this period. In it+1,
                % household will have x = net income
                asim(diesim(:,it+1)==1,it+1) = 0;
            end
        end
    end
    
    %% Moments/important quantities
    sim_results.mean_s          = mean(ssim(:,p.Tsim));
    sim_results.mean_a          = mean(asim(:,p.Tsim));
    sim_results.mean_x          = mean(xsim(:,p.Tsim));
    
    % Annual statistics
    if p.freq == 1
        sim_results.mean_grossy_A    = mean(ygrosssim(:,p.Tsim));
        sim_results.mean_loggrossy_A = mean(log(ygrosssim(:,p.Tsim)));
        sim_results.mean_nety_A      = mean(ynetsim(:,p.Tsim));
        sim_results.mean_lognety_A   = mean(log(ynetsim(:,p.Tsim)));
        sim_results.var_loggrossy_A  = var(log(ygrosssim(:,p.Tsim)));
        sim_results.var_lognety_A    = var(log(ynetsim(:,p.Tsim)));
        sim_results.wealthgini_A     = ginicoeff(asim(:,p.Tsim));
        sim_results.grossincgini_A   = ginicoeff(ygrosssim(:,p.Tsim));
        sim_results.netincgini_A     = ginicoeff(ynetsim(:,p.Tsim));
    else
        sim_results.mean_grossy_A    = mean(sum(ygrosssim(:,p.Tsim-3:p.Tsim),2));
        sim_results.mean_loggrossy_A = mean(log(sum(ygrosssim(:,p.Tsim),2)));
        sim_results.mean_nety_A      = mean(sum(ynetsim(:,p.Tsim-3:p.Tsim),2));
        sim_results.mean_lognety_A   = mean(log(ynetsim(:,p.Tsim)));
        sim_results.var_loggrossy_A  = var(log(ygrosssim(:,p.Tsim)));
        sim_results.var_lognety_A    = var(log(ynetsim(:,p.Tsim)));
        sim_results.wealthgini_A     = ginicoeff(asim(:,p.Tsim));
        sim_results.grossincgini_A   = ginicoeff(ygrosssim(:,p.Tsim));
        sim_results.netincgini_A     = ginicoeff(ynetsim(:,p.Tsim));
    end

    assetmeans = p.R * mean(ssim);
    
    % fraction constrained
    for i = 1:numel(p.epsilon)
        sim_results.constrained(i) = mean(asim(:,p.Tsim)<=p.borrow_lim+p.epsilon(i)*income.meany1*p.freq);
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
    
    sim_results.mpcs = simulation_MPCs(p,xsim,csim,diesim,ynetsim,...
                            yPindsim,yFindsim,betaindsim,income,model,xgrid);
    sim_results.assetmeans = assetmeans;

end