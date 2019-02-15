function [MPCs,stdev_loggrossy_A,stdev_lognety_A]... 
                                = direct_MPCs_by_simulation(p,prefs,income,basemodel,xgrid,agrid)
    % This function draws from the stationary distribution of (a,yP,yF,beta) 
    % and simulates 1-4 periods to find MPCs.
    
    if p.Display == 1
        disp([' Simulating ' num2str(p.freq*4) ' period(s) to get MPCs'])
    end
    
    % Number of draws from distribution
    Nsim = p.Nmpcsim;

    % Number of periods to simulate
    Tmax = p.freq * 4;

    % Vector of indexes for (yP,yF,beta) consistent with of mean ann inc
    yPind_trans = repmat(kron((1:p.nyP)',ones(p.nxlong,1)),p.nyF*p.nb,1);
    yFind_trans = repmat(kron((1:p.nyF)',ones(p.nxlong*p.nyP,1)),p.nb,1);
    
    if (numel(p.risk_aver) == 1) && (numel(p.invies) == 1) 
        betaind_trans = kron((1:p.nb)',ones(p.nxlong*p.nyP*p.nyF,1));
        IESind_trans = kron(ones(p.nb,1),ones(p.nxlong*p.nyP*p.nyF,1));
    else
        betaind_trans = kron(ones(p.nb,1),ones(p.nxlong*p.nyP*p.nyF,1));
        IESind_trans = kron((1:p.nb)',ones(p.nxlong*p.nyP*p.nyF,1));
    end

    % Construct stationary distribution
    state_rand  = rand(Nsim,1,'single');
    yPrand      = rand(Nsim,Tmax,'single');
    dierand     = rand(Nsim,Tmax,'single');
    betarand    = rand(Nsim,Tmax,'single');
    IESrand      = rand(Nsim,Tmax,'single');
    yTrand      = rand(Nsim,Tmax,'single');
    betaindsim  = ones(Nsim,Tmax,'int8');
    IESindsim = ones(Nsim,Tmax,'int8');
    yPindsim    = ones(Nsim,Tmax,'int8');
    yTindsim    = ones(Nsim,Tmax,'int8');
    yFindsim    = ones(Nsim,1,'int8');
    
    diesim      = dierand < p.dieprob;
    
    % Find (x,yPgrid,yFgrid,betagrid) indices of draws from stationary distribution
    % Done in partitions to economize on memory
    partitionsize = 1e5;
    Npartition = Nsim/partitionsize;
    cumdist = cumsum(basemodel.adist(:));
   
    a1 = zeros(Nsim,1);
    for ip = 1:Npartition
        partition = partitionsize*(ip-1)+1:partitionsize*ip;
        % Location of each draw in SSdist
        [~,ind] = max(bsxfun(@lt,state_rand(partition),cumdist'),[],2);
        
        % (yPgrid,yFgrid,betagrid) iindices
        yPindsim(partition,1)	= yPind_trans(ind);
        yFindsim(partition)     = yFind_trans(ind);
        betaindsim(partition,1)	= betaind_trans(ind);
        IESindsim(partition,1)  = IESind_trans(ind);
        
        % Initial assets from stationary distribution
        a1(partition) = agrid(ind);
    end
    
    %% SIMULATE INCOME AND BETA
    % Simulate frequency * 4 periods
    for it = 1:Tmax
        live = (diesim(:,it)==0);
        [~,yTindsim(:,it)]      = max(bsxfun(@le,yTrand(:,it),income.yTcumdist'),[],2);
        
        if it > 1
            [~,yPindsim(live,it)]   = max(bsxfun(@le,yPrand(live,it),income.yPcumtrans(yPindsim(live,it-1),:)),[],2);
            [~,yPindsim(~live,it)]  = max(bsxfun(@le,yPrand(~live,it),income.yPcumdist'),[],2);
            [~,betaindsim(:,it)]    = max(bsxfun(@le,betarand(:,it),prefs.betacumtrans(betaindsim(:,it-1),:)),[],2);
            [~,IESindsim(:,it)] = max(bsxfun(@le,IESrand(:,it),prefs.IEScumtrans(IESindsim(:,it-1),:)),[],2);
        end
    end
    
    ygrosssim = income.yPgrid(yPindsim) .*...
            repmat(income.yFgrid(yFindsim),1,Tmax) .* income.yTgrid(yTindsim);

    % Net income
    ynetsim = income.lumptransfer + (1-p.labtaxlow)*ygrosssim...
                        - p.labtaxhigh*max(ygrosssim-income.labtaxthresh,0);
                    
    % Find annual income variances
    stdev_loggrossy_A = std(log(sum(ygrosssim(:,1:p.freq),2)));
    stdev_lognety_A = std(log(sum(ynetsim(:,1:p.freq),2)));
    mean_grossy_A = mean(sum(ygrosssim(:,1:p.freq),2));
    
    
    % Loop over mpcfrac sizes, first running simulation as if there was no
    % shock
    for im = 0:numel(p.mpcfrac)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.mpcfrac(im) * income.meany1 * p.freq;
        end

        ssim = zeros(Nsim,Tmax,'single');
        asim = zeros(Nsim,Tmax,'single');
        xsim = zeros(Nsim,Tmax,'single');

        %% SIMULATE DECISION VARIABLES UPON MPC SHOCK
        
        % Location of households that get pushed below bottom of their grid
        % in first period upon shock
        set_mpc_one = false(Nsim,1);
        
        for it = 1:Tmax
            % Update cash-on-hand          
            if it == 1
                xsim(:,1) = a1 + + ynetsim(:,1) + mpcamount;
            else
                xsim(:,it) = asim(:,it) + ynetsim(:,it);
            end
            
            for ib = 1:p.nb
            for iyF = 1:p.nyF
            for iyP = 1:p.nyP
                if (numel(p.risk_aver) == 1) && (numel(p.invies) == 1) 
                    idx = yPindsim(:,it)==iyP & yFindsim==iyF & betaindsim(:,it)==ib;
                else
                    idx = yPindsim(:,it)==iyP & yFindsim==iyF & IESindsim(:,it)==ib;
                end
                
                if mpcamount < 0 && it == 1
                    below_grid = xsim(:,it)<xgrid.longgrid(1,iyP,iyF);
                    % Bring households pushed below grid back up to grid
                    idx_below = idx & below_grid;
                    xsim(idx_below,it) = xgrid.longgrid(1,iyP,iyF);
                    % Update set_mpc_one
                    set_mpc_one = set_mpc_one | idx_below;
                end
                ssim(idx,it) = basemodel.savinterp{iyP,iyF,ib}(xsim(idx,it));
            end
            end
            end
            
            ssim(:,it) = max(ssim(:,it),p.borrow_lim);

            if it < Tmax
                asim(:,it+1) = p.R * ssim(:,it);
                if p.Bequests == 0
                    % Assets discarded
                    asim(diesim(:,it+1)==1,it+1) = 0;
                end
            end
        end
        
        csim = xsim - ssim - p.savtax * max(ssim-p.savtaxthresh,0);
        
        %% COMPUTE MPCs
        if im == 0
            % No MPC schock
            csim_noshock = csim;
        else
        	% MPC in period 1 out of period 1 shock
            mpcs_1_1 = (csim(:,1) - csim_noshock(:,1)) / mpcamount;
            mpcs_1_1(set_mpc_one,1) = 1;  
            MPCs.avg_1_1(im) = mean(mpcs_1_1);

            % MPC in period 2 out of period 1 shock
            mpcs_1_2 = (csim(:,2) - csim_noshock(:,2)) / mpcamount;
            mpcs_1_2(set_mpc_one) = 0;
            MPCs.avg_1_2(im) = mean(mpcs_1_2);

            % MPC in period 3 out of period 1 shock
            mpcs_1_3 = (csim(:,3) - csim_noshock(:,3)) / mpcamount;
            mpcs_1_3(set_mpc_one) = 0;
            MPCs.avg_1_3(im) = mean(mpcs_1_3);

            % MPC in period 4 out of period 1 shock
            mpcs_1_4 = (csim(:,4) - csim_noshock(:,4)) / mpcamount;
            mpcs_1_4(set_mpc_one) = 0;
            MPCs.avg_1_4(im) = mean(mpcs_1_4);

            % Cumulative MPCs over first 4 periods
            mpcs_1_1to4 = mpcs_1_1+mpcs_1_2+mpcs_1_3+mpcs_1_4;
            MPCs.avg_1_1to4(im) = mean(mpcs_1_1to4);

            if p.freq == 4
            	% MPC in period ip out of period 1 shock
                mpcs_1_x = cell(1,16);
            	for ip = 5:16
            		mpcs_1_x{ip} = (csim(:,ip) - csim_noshock(:,ip)) / mpcamount;
            		mpcs_1_x{ip}(set_mpc_one) = 0;
            	end

            	% Cumulative mean MPC over years 2-4 for quarterly model
            	MPCs.avg_1_5to8(im) = mean(mpcs_1_x{5}+mpcs_1_x{6}+mpcs_1_x{7}+mpcs_1_x{8});
            	MPCs.avg_1_9to12(im) = mean(mpcs_1_x{9}+mpcs_1_x{10}+mpcs_1_x{11}+mpcs_1_x{12});
            	MPCs.avg_1_13to16(im) = mean(mpcs_1_x{13}+mpcs_1_x{14}+mpcs_1_x{15}+mpcs_1_x{16});
                clear mpcs_1_x
            else
            	MPCs.avg_1_5to8(im) = NaN;
            	MPCs.avg_1_9to12(im) = NaN;
            	MPCs.avg_1_13to16(im) = NaN;
            end
        end
    end
    
end