function [xsim ssim grossysim ynetsim] = simulations(p,yTcumdist,yFcumdist,...
    yPcumdist,yPcumtrans,yPgrid,yFgrid,yTgrid,labtaxthresh,conm,savm,xgridm,...
    lumptransfer,betacumdist,betacumtrans)
    
    
    grossysim = zeros(p.Nsim,p.Tsim);
    ynetsim = zeros(p.Nsim,p.Tsim);
    
    %% Simulate income process
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
    for it = 1:p.Tsim
        ygrosssim(:,it) = ygrosssim(:,it)/mean(ygrosssim(:,it));
    end
    
    % net income
    ynetsim = lumptransfer + (1-p.labtaxlow)*ygrosssim - p.labtaxhigh*max(ygrosssim-labtaxthresh,0);
    
    %% Simulate beta
    betarand = rand(p.Nsim,p.Tsim);
    betaindsim = zeros(p.Nsim,p.Tsim);
    [~,betaindsim(:,1)] = max(bsxfun(@le,betarand(:,1),betacumdist'),[],2);
    
    for it = 2:p.Tsim
        ilive = diesim(:,it)==0;
        [~,betaindsim(ilive,it)] = max(bsxfun(@le,betarand(ilive,it),betacumtrans(betaindsim(ilive,it-1),:)),[],2);
        [~,betaindsim(~ilive,it)] = max(bsxfun(@le,betarand(~ilive,it),betacumdist'),[],2);
    end
    
    %% Simulate savings decisions
    xsim = zeros(p.Nsim,p.Tsim); 
    ssim = zeros(p.Nsim,p.Tsim);
    
    for iyF = 1:p.nyF
    for ib = 1:p.nb
    for iyP = 1:p.nyP
        coninterp{iyP,ib,iyF} = griddedInterpolant(xgridm(:,iyP,iyF),conm(:,iyP,ib,iyF),'linear');
        savinterp{iyP,ib,iyF} = griddedInterpolant(xgridm(:,iyP,iyF),savm(:,iyP,ib,iyF),'linear');
    end
    end
    end
    
    for it = 1:p.Tsim
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
    
    plot(mean(ssim));
    mean(ssim(:,p.Tsim))
    mean(xsim(:,p.Tsim))
    mean(ssim(:,p.Tsim)<=p.borrow_lim)



end