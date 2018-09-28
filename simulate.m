function [simulations,ssim] = simulate(p,income,labtaxthresh,model,...
                                        xgrid,lumptransfer,prefs,results)
    
    
    grossysim = zeros(p.Nsim,p.Tsim);
    ynetsim = zeros(p.Nsim,p.Tsim);
    
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
    ynetsim = lumptransfer + (1-p.labtaxlow)*ygrosssim - p.labtaxhigh*max(ygrosssim-labtaxthresh,0);
    
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
            ssim(idx,it) = model.savinterp{iyP,ib,iyF}(xsim(idx,it));
        end
        end
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
    
    
%     Nmpcamount = numel(p.mpcfrac);
%     for im = 1:Nmpcamount
%         mpcamount{im} = p.mpcfrac{im} * income.meany;
%         % period Tsim cash-on-hand
%         xsim_mpc{im} = xsim(:,p.Tsim) + mpcamount{im};
%         
%         csim_mpc{im} = zeros(p.Nsim,1);
%         
%         xsim_mpc{im} = zeros(p.Nsim,4);
%         ssim_mpc{im} = zeros(p.Nsim,4);
%         xsim_mpc{im}(:,1) = xsim(:,p.Tsim-3) + mpcamount{im};
%         csim_mpc{im} = zeros(p.Nsim,4);
%         for it = 1:4
%             if it > 1
%                 xsim_mpc{im}(:,it) = p.R * ssim_mpc{im}(:,it-1) + ynetsim(:,p.Tsim-4+it);
%             end
%             
%             for iyF = 1:p.nyF
%             for ib = 1:p.nb
%             for iyP = 1:p.nyP
%                 idx = yPindsim(:,p.Tsim-4+it)==iyP & betaindsim(:,p.Tsim-4+it)==ib & yFindsim(:,p.Tsim-4+it)==iyF;
%                 ssim_mpc{im}(idx,it) = savinterp{iyP,ib,iyF}(xsim_mpc{im}(idx,it));
%             end
%             end
%             end
%             ssim_mpc{im}(ssim_mpc{im}(:,it)<p.borrow_lim,it) = p.borrow_lim;
%             csim_mpc{im}(:,it) = xsim_mpc{im}(:,it) - ssim_mpc{im}(:,it) - p.savtax*max(ssim_mpc{im}(:,it)-p.savtaxthresh,0);
%             if p.WealthInherited == 0
%                 % set saving equal to 0 if hh dies at end of this period. In it+1,
%                 % household will have x = net income
%                 ssim_mpc{im}(diesim(:,p.Tsim-4+it)==1,it) = 0;
%             end
%         end
%         
%         mpc1{im} = (csim_mpc{im}(:,1) - csim(:,p.Tsim-3))/mpcamount{im};
%         mpc2{im} = (csim_mpc{im}(:,2) - csim(:,p.Tsim-2))/mpcamount{im};
%         mpc3{im} = (csim_mpc{im}(:,3) - csim(:,p.Tsim-1))/mpcamount{im};
%         mpc4{im} = (csim_mpc{im}(:,4) - csim(:,p.Tsim))/mpcamount{im};
%         simulations.avg_mpc1{im} = mean(mpc1{im});
%         simulations.avg_mpc4{im} = mean(mpc1{im}+mpc2{im}+mpc3{im}+mpc4{im});
%     end
    Risk = 1;
    [simulations.avg_mpc1,simulations.avg_mpc4] = simulation_MPCs(p,xsim,csim,diesim,ynetsim,yPindsim,yFindsim,...
                    betaindsim,income,Risk,model,xgrid);
    
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