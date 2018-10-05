function [mpc1,mpc4,avg_mpc1,avg_mpc4,var_mpc1,var_mpc4] ... 
                                    = direct_MPCs(p,prefs,income,basemodel,xgrid)

    if p.Display == 1
        disp(' Simulating 4 periods to get MPCs')
    end
    
    
    
    % vector of indexes for (yP,yF,beta)
    yPind_trans = repmat(kron((1:p.nyP)',ones(p.nxlong,1)),p.nyF*p.nb,1);
    yFind_trans = repmat(kron((1:p.nyF)',ones(p.nxlong*p.nyP,1)),p.nb,1);
    betaind_trans = kron((1:p.nb)',ones(p.nxlong*p.nyP*p.nyF,1));

    temp = sortrows([basemodel.sav_longgrid yPind_trans...
                            yFind_trans betaind_trans xgrid.longgrid]);
    yPind_trans     = temp(:,2);
    yFind_trans     = temp(:,3);
    betaind_trans   = temp(:,4);
    cashgrid        = temp(:,5);
    
    Nsim = 1e6;
    
    % Construct stationary distribution
    state_rand  = rand(Nsim,4);
    yPrand      = rand(Nsim,4);
    dierand     = rand(Nsim,4);
    yTrand      = rand(Nsim,4);
    ygrosssim   = zeros(Nsim,4);
    ynetsim     = zeros(Nsim,4);
    betaindsim  = ones(Nsim,4);
    yPindsim    = ones(Nsim,4);
    yTindsim    = ones(Nsim,4);
    yFindsim    = ones(Nsim,1);
    x0          = zeros(Nsim,1);

    diesim      = dierand < p.dieprob;
    
    % index of first state
    partitionsize = 1e5;
    Npartition = Nsim/partitionsize;
    for ip = 1:Npartition
        partition = partitionsize*(ip-1)+1:partitionsize*ip;
        [~,ind] = max(bsxfun(@lt,state_rand(partition,1),basemodel.SScumdist'),[],2);
        yPindsim(partition,1)	= yPind_trans(ind);
        yFindsim(partition)     = yFind_trans(ind);
        betaindsim(partition,1)	= betaind_trans(ind);
        
        % initial cash-on-hand
        x0(partition) = cashgrid(ind);
    end
    
    
    
    %% SIMULATE INCOME
    
    % iterate over time periods
    for it = 2:4
        live = (diesim(:,it)==0);
    	[~,yPindsim(live,it)]   = max(bsxfun(@le,yPrand(live,it),income.yPcumtrans(yPindsim(live,it-1),:)),[],2);
        [~,yPindsim(~live,it)]  = max(bsxfun(@le,yPrand(~live,it),income.yPcumdist'),[],2);
        
        [~,betaindsim(:,it)]    = max(bsxfun(@le,betaindsim(:,it),prefs.betacumtrans(betaindsim(:,it-1),:)),[],2);
        [~,yTindsim(:,it)]      = max(bsxfun(@le,yTrand(:,it),income.yTcumdist'),[],2);
    end
    
    % First column is irrelevant
    ygrosssim = income.yPgrid(yPindsim) .*...
            repmat(income.yFgrid(yFindsim),1,4) .* income.yTgrid(yTindsim);
        
    % Normalize so mean annual income == 1
    if p.NormalizeY == 1
        ygrosssim = ygrosssim / (income.original_meany * p.freq);
    end
    
    ynetsim = income.lumptransfer + (1-p.labtaxlow)*ygrosssim...
                        - p.labtaxhigh*max(ygrosssim-income.labtaxthresh,0);
    
    % Loop over mpcfrac sizes, first running simulation as if there was no
    % shock
    mpcamounts = [0 ; p.mpcfrac'*income.meany*p.freq];
    for im = 1:numel(mpcamounts)
        mpcamount = mpcamounts(im);
        im = im - 1;
        
        xsim = zeros(Nsim,4);
        ssim = zeros(Nsim,4);
        asim = zeros(Nsim,4);
        
        %% SIMULATE DECISION VARIABLES UPON MPC SHOCK
        
        set_mpc_one = false(Nsim,1);
        
        for it = 1:4
            % Update cash-on-hand
            if it == 1
                xsim(:,it) = x0 + mpcamount;
            else
                xsim(:,it) = asim(:,it-1) + ynetsim(:,it);
            end
            
            for ib = 1:p.nb
            for iyF = 1:p.nyF
            for iyP = 1:p.nyP
                idx = yPindsim(:,it)==iyP & yFindsim==iyF & betaindsim(:,it)==ib;
                if mpcamount < 0 && it == 1
                    below_grid = xsim(:,it)<xgrid.longgrid_wide(1,iyP,iyF);
                    xsim(set_mpc_one,it) = xgrid.longgrid_wide(1,iyP,iyF);
                    idx_below = idx & below_grid;
                    set_mpc_one = set_mpc_one | idx_below;
                end
                ssim(idx,it) = basemodel.savinterp{iyP,iyF,ib}(xsim(idx,it));
            end
            end
            end
            
            ssim(ssim(:,it)<p.borrow_lim,it) = p.borrow_lim;

            asim(:,it) = p.R * ssim(:,it);
            if p.WealthInherited == 0
                asim(diesim(:,it)==1,it+1) = 0;
            end
        end
        
        csim = xsim - ssim - p.savtax * max(ssim-p.savtaxthresh,0);
        
        %% COMPUTE MPCs
        if im == 0
            % No MPC schock
            csim_noshock = csim;
        else
            mpc1{im} = (csim(:,1) - csim_noshock(:,1)) / mpcamount;
            mpc1{im}(set_mpc_one,1) = 1;
            mpc2{im} = (csim(:,2) - csim_noshock(:,2)) / mpcamount + mpc1{im};
            mpc3{im} = (csim(:,3) - csim_noshock(:,3)) / mpcamount + mpc2{im};
            mpc4{im} = (csim(:,4) - csim_noshock(:,4)) / mpcamount + mpc3{im};
            
            avg_mpc1{im} = mean(mpc1{im});
            avg_mpc4{im} = mean(mpc4{im});
            
            var_mpc1{im} = var(mpc1{im});
            var_mpc4{im} = var(mpc4{im});
        end
    end

end