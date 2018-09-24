function [simulations ssim] = simulate(p,yTcumdist,yFcumdist,...
    yPcumdist,yPcumtrans,yPgrid,yFgrid,yTgrid,labtaxthresh,conm,savm,xgridm,...
    lumptransfer,betacumdist,betacumtrans,original_meany)
    
    
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

    [~,yFindsim] = max(bsxfun(@le,yFrand,yFcumdist'),[],2);
    yFindsim = repmat(yFindsim,1,p.Tsim);
    
    for it = 1:p.Tsim
        [~,yTindsim(:,it)] = max(bsxfun(@le,yTrand(:,it),yTcumdist'),[],2);
        
        ilive = diesim(:,it)==0;
        
        if it ==1
            [~,yPindsim(ilive,it)] = max(bsxfun(@le,yPrand(ilive,it),yPcumdist'),[],2);
        else
            [~,yPindsim(ilive,it)] = max(bsxfun(@le,yPrand(ilive,it),yPcumtrans(yPindsim(ilive,it-1),:)),[],2);
        end
        
        [~,yPindsim(~ilive,it)] = max(bsxfun(@le,yPrand(~ilive,it),yPcumdist'),[],2);
    end
    
    % gross income
    ygrosssim = yPgrid(yPindsim).*yFgrid(yFindsim).*yTgrid(yTindsim);
    if p.NormalizeY == 1
        ygrosssim = ygrosssim/original_meany;
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
        savinterp{iyP,ib,iyF} = griddedInterpolant(xgridm(:,iyP,iyF,ib),savm(:,iyP,iyF,ib),'linear');
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