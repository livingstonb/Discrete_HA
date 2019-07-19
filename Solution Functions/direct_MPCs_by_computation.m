function [MPCs,agrid_dist] = direct_MPCs_by_computation(p,basemodel,models,income,prefs,grids,shocksize)
    % This function computes IMPC(s,t) using transition probabilities,
    % where IMPC(s,t) is the MPC in period t out of a period-s shock that
    % is learned about in period 1

    % mpc out of negative shock is slightly innacurate for large shock
    
    NN = p.nx_KFE*p.nyP*p.nyF*p.nb;
    
    if p.Display == 1
        disp('Computing MPCs')
    end

    if numel(p.r) > 1
        r_col = kron(p.r',ones(p.nx_KFE*p.nyP*p.nyF,1));
        r_mat = reshape(r_col,[p.nx_KFE,p.nyP,p.nyF,numel(p.r)]);
    else
        r_mat = p.r;
    end
    
    %% --------------------------------------------------------------------
    % TRANSITION PROBABILITES AND BASELINE CONSUMPTION
    % ---------------------------------------------------------------------
    
    % Each (a,yP,yF) is associated with nyT possible x values, create this
    % grid here
    xgrid_yT = repmat(grids.a.vec,[1 p.nyP p.nyF p.nyT]) + income.netymatKFE;

    % for interpolating back onto agrid
    fspace = fundef({'spli',grids.a.vec,0,1});
    
    % baseline consumption
    con_baseline_yT = get_policy(p,xgrid_yT,basemodel);
    % take expectation wrt yT
    con_baseline = reshape(con_baseline_yT,[],p.nyT) * income.yTdist;
    
    mpcamount = shocksize;

    %% --------------------------------------------------------------------
    % ITERATION OVER (s,t)
    % ---------------------------------------------------------------------
    % period t MPC out of shock in period s, learned about in period 1
    MPCs.avg_s_t = cell(5,5);
    
    for is = 1:5
    for it = 1:5
        MPCs.avg_s_t{is,it} = NaN;
    end
    end
    
    % collect mpcs as a function of state space
    MPCs.mpcs_1_t = cell(1,4);
    
    % which periods to apply shock
    if (shocksize < 0) || (p.mpcshocks_after_period1 == 0)
        IS = 1;
    else
        IS = [1 2 5];
    end
    
    for is = IS
        % iterate over t within s

        if (p.EpsteinZin == 1) && (is > 1)
  			% not interested in this
            continue
        end

        for it = 1:is % in this block, look at t <= s
            if p.Display == 1
                fprintf('    Processing (s,t) = (%d,%d)\n',is,it)
                disp(['    --Time ' datestr(now,'HH:MM:SS')])
            end

            if (shocksize < 0) && (it > 1)
                continue
            end

            % Create transition matrix from period 1 to period
            % t (for last iteration, this is transition from period t to
            % period s)
            if it == 1
                T1t = speye(NN);
            else
                %Ti is transition from t=ii-1 to t=ii
                mpcshock = 0;
                Ti = transition_t_less_s(p,income,xgrid_yT,...
                            models,is,it-1,fspace,mpcshock,r_mat);
                T1t = Ti * T1t;
            end

            % get consumption policy function
            if it == is
                % shock in this period
                x_mpc = xgrid_yT + shocksize;
            else
            	% no shock
                x_mpc = xgrid_yT;
            end

            if shocksize < 0
            	% record which states are pushed below asset grid after negative shock
                % bring back up to asset grid for interpolation
                below_xgrid = false(size(x_mpc));
                for iyT = 1:p.nyT
                    below_xgrid(:,:,:,iyT) = x_mpc(:,:,:,iyT) < grids.x.matrix(1,:,:);

                    x_mpc(:,:,:,iyT) = ~below_xgrid(:,:,:,iyT) .* x_mpc(:,:,:,iyT)...
                                        + below_xgrid(:,:,:,iyT) .* grids.x.matrix(1,:,:);
                end
                below_xgrid = reshape(below_xgrid,[p.nx_KFE p.nyP p.nyF 1 p.nyT]);
                below_xgrid = repmat(below_xgrid,[1 1 1 p.nb 1]);
            end

            % consumption choice given the shock
            con = get_policy(p,x_mpc,models{is,it});

            if (shocksize < 0) && (it == is)
                % make consumption for cases pushed below xgrid equal to consumption
                % at bottom of xgrid - the amount borrowed
                x_before_shock = reshape(grids.x.matrix,[p.nx_KFE p.nyP p.nyF]);
                x_minus_xmin = x_before_shock - grids.x.matrix(1,:,:);
            	con = ~below_xgrid.*con + below_xgrid.*(con_baseline_yT(1,:,:,:,:) + shocksize + x_minus_xmin);
            end

            % expectation over yT
            con = reshape(con,[],p.nyT) * income.yTdist;

            % now compute IMPC(s,t)
            mpcs = ( T1t * con - con_baseline) / mpcamount;

            MPCs.avg_s_t{is,it} = basemodel.adist(:)' * mpcs(:);
            
            if (is == 1) && (it >= 1 && it <= 4)
            	% store is = 1 mpcs
                MPCs.mpcs_1_t{it} = mpcs;
            end

        end

        % transition probabilities from it = is to it = is + 1
        mpcshock = mpcamount;
        T_s1_s = transition_t_less_s(p,income,xgrid_yT,models,is,is,...
                                    fspace,mpcshock,r_mat);
        T1t = T_s1_s * T1t;
        clear T_s1_s

        RHScon = basemodel.statetrans * con_baseline(:);
        LHScon = con_baseline(:);

        for it = is+1:5 % it > is case, policy fcns stay the same in this region
            mpcs = (T1t * LHScon(:) - RHScon(:)) / mpcamount;

            MPCs.avg_s_t{is,it} = basemodel.adist(:)' * mpcs(:);
            RHScon = basemodel.statetrans * RHScon;
            LHScon = basemodel.statetrans * LHScon;
            
            if (is == 1) && (it >= 1 && it <= 4)
                MPCs.mpcs_1_t{it} = mpcs;
            end
        end
    end
    

    %% --------------------------------------------------------------------
    % FIND CUMULATIVE MPCs 1-4, 5-8, 9-12, 13-16
    % ---------------------------------------------------------------------
    % shock in period 1
    MPCs.avg_1_1to4 = cellsum(MPCs.avg_s_t,1,1,4);

    % shock in period 5
    MPCs.avg_5_1to4 = cellsum(MPCs.avg_s_t,5,1,4);

    %% --------------------------------------------------------------------
    % DISTRIBUTION OVER AGRID, P(a)
    % ---------------------------------------------------------------------
    agrid_dist = sum(sum(sum(basemodel.adist,4),3),2);
end

%% ------------------------------------------------------------------------
% SUBFUNCTIONS
% -------------------------------------------------------------------------

function con = get_policy(p,x_mpc,model)
    % Computes consumption policy function after taking expectation to get
    % rid of yT dependence
    
    sav = zeros(p.nx_KFE,p.nyP,p.nyF,p.nb,p.nyT);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        x_iyP_iyF_iyT = x_mpc(:,iyP,iyF,:);
        sav_iyP_iyF_iyT = model.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT,[p.nx_KFE 1 1 1 p.nyT]);
    end
    end
    end
    sav = max(sav,p.borrow_lim);
    x_mpc = reshape(x_mpc,[p.nx_KFE p.nyP p.nyF 1 p.nyT]);
    con = repmat(x_mpc,[1 1 1 p.nb 1]) - sav - p.savtax * max(sav-p.savtaxthresh,0);
end

function T1 = transition_t_less_s(p,income,xgrid_yT,models,is,ii,...
                                            fspace,mpcshock,r_mat)
    % Computes the transition matrix between t=ii and t=ii + 1 given shock in
    % period 'is'
    NN = p.nx_KFE*p.nyP*p.nyF*p.nb;
    x_mpc = xgrid_yT + mpcshock; % cash on hand after receiving shock
   
    sav = zeros(p.nx_KFE,p.nyP,p.nyF,p.nb,p.nyT);

    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        x_iyP_iyF_iyT = reshape(x_mpc(:,iyP,iyF,:),[],1);
        
        if mpcshock < 0 && (ii == is)
        	below_xgrid = x_iyP_iyF_iyT < min(xgrid_yT(1,iyP,iyF,:));
            x_iyP_iyF_iyT(below_xgrid) = min(xgrid_yT(1,iyP,iyF,:));
        end
        
        sav_iyP_iyF_iyT = models{is,ii}.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT);
        if mpcshock < 0 && (ii == is)
            sav_iyP_iyF_iyT(below_xgrid) = 0;
        end

        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT,[p.nx_KFE 1 1 1 p.nyT]);

    end
    end
    end

    sav = max(sav,p.borrow_lim);

    aprime_live = (1+repmat(r_mat,[1,1,1,1,p.nyT])) .* sav; % next period's assets

    % interpolate next period's assets back onto asset grid
    interp = funbas(fspace,aprime_live(:));
    interp = reshape(interp,NN,p.nx_KFE*p.nyT) * kron(speye(p.nx_KFE),income.yTdist);
    if p.Bequests == 1
        interp_death = interp;
    else
        interp_death = sparse(NN,p.nx_KFE);
        interp_death(:,1) = 1;
    end

    % Transition matrix from 'ii' to 'ii'+1
    T1 = sparse(NN,NN);
    for col = 1:p.nyP*p.nyF*p.nb
        newblock_live = bsxfun(@times,kron(income.ytrans_live(:,col),ones(p.nx_KFE,1)),interp);
        newblock_death = bsxfun(@times,kron(income.ytrans_death(:,col),ones(p.nx_KFE,1)),interp_death);
        T1(:,p.nx_KFE*(col-1)+1:p.nx_KFE*col) = (1-p.dieprob)*newblock_live + p.dieprob*newblock_death;
    end
end

function out = cellsum(array,index,start,last)

    out = 0;
    for i = start:last
        out = out + array{index,i};
    end
end

