function [asim1,mpcs1,mpcs4] = direct_MPCs(p,prefs,income,basemodel,xgrid)

    if p.Display == 1
        disp([' Simulating ' num2str(4) ' period(s) to get MPCs'])
    end
    
    % Number of draws from distribution
    Nsim = 1e6;

    % vector of indexes for (yP,yF,beta) consistent with of mean ann inc
    yPind_trans = repmat(kron((1:p.nyP)',ones(p.nxlong,1)),p.nyF*p.nb,1);
    yFind_trans = repmat(kron((1:p.nyF)',ones(p.nxlong*p.nyP,1)),p.nb,1);
    betaind_trans = kron((1:p.nb)',ones(p.nxlong*p.nyP*p.nyF,1));

    cumdist = cumsum(basemodel.SSdist);

    % Construct stationary distribution
    state_rand  = rand(Nsim,4);
    yPrand      = rand(Nsim,4);
    dierand     = rand(Nsim,4);
    betarand    = rand(Nsim,4);
    yTrand      = rand(Nsim,4);
    ygrosssim   = zeros(Nsim,4);
    ynetsim     = zeros(Nsim,4);
    betaindsim  = ones(Nsim,4);
    yPindsim    = ones(Nsim,4);
    yTindsim    = ones(Nsim,4);
    yFindsim    = ones(Nsim,1);
    x0          = zeros(Nsim,1);

    diesim      = dierand < p.dieprob;
    
    % Find (yPgrid,yFgrid,betagrid) indices of draws from stationary distribution
    % Done in partitions to economize on memory
    partitionsize = 1e5;
    Npartition = Nsim/partitionsize;
    xgrid_extended = repmat(xgrid.longgrid,p.nb,1);
    for ip = 1:Npartition
        partition = partitionsize*(ip-1)+1:partitionsize*ip;
        % Location of each draw in SSdist
        [~,ind] = max(bsxfun(@lt,state_rand(partition,1),cumdist'),[],2);
        
        % (yPgrid,yFgrid,betagrid) indices
        yPindsim(partition,1)	= yPind_trans(ind);
        yFindsim(partition)     = yFind_trans(ind);
        betaindsim(partition,1)	= betaind_trans(ind);
        
        % Initial cash-on-hand from stationary distribution
        x0(partition) = xgrid_extended(ind);
    end
    
    
    
    %% SIMULATE INCOME AND BETA
    
    if 4 == 4
        for it = 2:4
            live = (diesim(:,it)==0);
            [~,yPindsim(live,it)]   = max(bsxfun(@le,yPrand(live,it),income.yPcumtrans(yPindsim(live,it-1),:)),[],2);
            [~,yPindsim(~live,it)]  = max(bsxfun(@le,yPrand(~live,it),income.yPcumdist'),[],2);
            [~,betaindsim(:,it)]    = max(bsxfun(@le,betarand(:,it),prefs.betacumtrans(betaindsim(:,it-1),:)),[],2);
            [~,yTindsim(:,it)]      = max(bsxfun(@le,yTrand(:,it),income.yTcumdist'),[],2);
        end
    end
    
    % Column 1 is irrelevant
    ygrosssim = income.yPgrid(yPindsim) .*...
            repmat(income.yFgrid(yFindsim),1,4) .* income.yTgrid(yTindsim);
        
    % Normalize so mean annual income == 1
    if p.NormalizeY == 1
        ygrosssim = ygrosssim / (income.original_meany * p.freq);
    end
    
    % Net income
    ynetsim = income.lumptransfer + (1-p.labtaxlow)*ygrosssim...
                        - p.labtaxhigh*max(ygrosssim-income.labtaxthresh,0);
    
    % Loop over mpcfrac sizes, first running simulation as if there was no
    % shock
    for im = 0:numel(p.mpcfrac)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.mpcfrac(im) * income.meany * p.freq;
        end
        
        xsim = zeros(Nsim,4);
        ssim = zeros(Nsim,4);
        asim = zeros(Nsim,4);
        
        %% SIMULATE DECISION VARIABLES UPON MPC SHOCK
        
        % Location of households that get pushed below bottom of their grid
        % in first period upon shock
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
                    % Bring households pushed below grid back up to grid
                    idx_below = idx & below_grid;
                    xsim(idx_below,it) = xgrid.longgrid_wide(1,iyP,iyF);
                    % Update set_mpc_one
                    set_mpc_one = set_mpc_one | idx_below;
                end
                ssim(idx,it) = basemodel.savinterp{iyP,iyF,ib}(xsim(idx,it));
            end
            end
            end
            
            ssim(ssim(:,it)<p.borrow_lim,it) = p.borrow_lim;

            asim(:,it) = p.R * ssim(:,it);
            if p.WealthInherited == 0
                % Assets discarded
                asim(diesim(:,it)==1,it) = 0;
            end
        end
        
        csim = xsim - ssim - p.savtax * max(ssim-p.savtaxthresh,0);
        
        if im == 5
            asim1 = asim(:,1);
        end
        
        %% COMPUTE MPCs
        if im == 0
            % No MPC schock
            csim_noshock = csim;
        else
            mpcs1{im} = (csim(:,1) - csim_noshock(:,1)) / mpcamount;
            mpcs1{im}(set_mpc_one,1) = 1;  
            mpcs2{im} = (csim(:,2) - csim_noshock(:,2)) / mpcamount + mpcs1{im};
            mpcs3{im} = (csim(:,3) - csim_noshock(:,3)) / mpcamount + mpcs2{im};
            mpcs4{im} = (csim(:,4) - csim_noshock(:,4)) / mpcamount + mpcs3{im};
        end
    end
    
end