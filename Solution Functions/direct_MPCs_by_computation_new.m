function [MPCs,agrid_dist] = direct_MPCs_by_computation_new(p,basemodel,models,income,prefs,agrid_short)
    % This function computes IMPC(s,t) using transition probabilities,
    % where IMPC(s,t) is the MPC in period t out of a period-s shock that
    % is learned about in period 1
    
    NN = p.nxlong*p.nyP*p.nyF*p.nb;
    
    if p.Display == 1
        disp('Computing MPCs')
    end
    
    %% --------------------------------------------------------------------
    % TRANSITION PROBABILITES AND BASELINE CONSUMPTION
    % ---------------------------------------------------------------------
    
    % Each (a,yP,yF) is associated with nyT possible x values, create this
    % grid here
    netymat_reshape = reshape(income.netymat,[1 p.nyP p.nyF p.nyT]);
    netymat_reshape = repmat(netymat_reshape,[p.nxlong 1 1 1]);
    xgrid_yT = repmat(agrid_short,[1 p.nyP p.nyF p.nyT]) + netymat_reshape;

    fspace = fundef({'spli',agrid_short,0,1});
    % transition matrix between (yP,yF,beta) states, cond'l on living
    trans_live = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    % transition matrix between (yP,yF,beta) states cond'l on dying
    yPtrans_stationary = repmat(income.yPdist',p.nyP,1);
    trans_death = kron(prefs.betatrans,kron(eye(p.nyF),yPtrans_stationary));
    
    con_baseline = get_policy(p,xgrid_yT,basemodel,income);
    
    mpcamount = 0.01*income.meany1*p.freq;

    %% --------------------------------------------------------------------
    % ITERATION OVER (s,t)
    % ---------------------------------------------------------------------
    MPCs.avg_s_t = cell(16,16);
    maxT = p.freq * 4;
    for is = 1:maxT
        % iterate over t within s
        for it = 1:is % in this block, look at t <= s
            % Create transition matrix from period 1 to period
            % t (for last iteration, this is transition from period t to
            % period s)
            for ii = 1:it
                if ii == 1
                    T1t = speye(NN);
                else
                    %Ti is transition from t=ii-1 to t=ii
                    mpcshock = 0;
                    Ti = transition_t_less_s(p,income,xgrid_yT,...
                        models,is,ii-1,fspace,trans_live,trans_death,mpcshock);
                    T1t = T1t * Ti;
                end
            end

            % get consumption policy function
            if it == is
                % shock is 1% of annual income
                x_mpc = xgrid_yT + 0.01*income.meany1*p.freq;
            else
                x_mpc = xgrid_yT;
            end
            con = get_policy(p,x_mpc,models{is,it},income);

            % now compute IMPC(s,t)
            mpcs = ( T1t * con(:) - con_baseline(:) ) / mpcamount;
            MPCs.avg_s_t{is,it} = basemodel.adist(:)' * mpcs(:);

        end

        % transition probabilities from it = is to it = is + 1
        mpcshock = mpcamount;
        T_s1_s = transition_t_less_s(p,income,xgrid_yT,models,is,is,...
                                    fspace,trans_live,trans_death,mpcshock);

        RHScon = con_baseline(:);
        for it = is+1:maxT % it > is case, policy fcns stay the same in this region
            mpcs = ( T1t * T_s1_s * RHScon(:) - con_baseline(:) ) / mpcamount;
            MPCs.avg_s_t{is,it} = basemodel.adist(:)' * mpcs(:);
            RHScon = basemodel.statetrans * RHScon;
        end
    end
    
    if p.freq == 1
        for is = 1:16
        for it = 1:16
            if (is>4) || (it>4)
                MPCs.avg_s_t{is,it} = NaN;
            end
        end
        end
    end

    % Distribution over agrid, P(a)
    agrid_dist = sum(sum(sum(basemodel.adist,4),3),2);
end

%% ------------------------------------------------------------------------
% SUBFUNCTIONS
% -------------------------------------------------------------------------

function Econ = get_policy(p,x_mpc,model,income)
    % Computes consumption policy function after taking expectation to get
    % rid of yT dependence
    
    sav = zeros(p.nxlong,p.nyP,p.nyF,p.nb,p.nyT);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        x_iyP_iyF_iyT = x_mpc(:,iyP,iyF,:);
        sav_iyP_iyF_iyT = model.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT,[p.nxlong 1 1 1 p.nyT]);
    end
    end
    end
    sav = max(sav,p.borrow_lim);

    % Take expectation over yT
    % con becomes E[con(x,yP,yF,beta)|a,yP,yF,beta]
    Esav = reshape(sav,[],p.nyT) * income.yTdist;
    Esav = reshape(Esav,[p.nxlong p.nyP p.nyF p.nb]);
    Exx = reshape(reshape(x_mpc,[],p.nyT) * income.yTdist,[p.nxlong p.nyP p.nyF]);
    Econ = repmat(Exx,[1 1 1 p.nb]) - Esav - p.savtax * max(Esav - p.savtaxthresh,0);  
end

function T1 = transition_t_less_s(p,income,xgrid_yT,models,is,ii,...
                                            fspace,trans_live,trans_death,mpcshock)
    % Computes the transition matrix between t=ii and t=is
    NN = p.nxlong*p.nyP*p.nyF*p.nb;
    x_mpc = xgrid_yT + mpcshock;
    sav = zeros(p.nxlong,p.nyP,p.nyF,p.nb,p.nyT);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        x_iyP_iyF_iyT = x_mpc(:,iyP,iyF,:);
        sav_iyP_iyF_iyT = models{is,ii}.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT,[p.nxlong 1 1 1 p.nyT]);
    end
    end
    end
    sav = max(sav,p.borrow_lim);

    aprime_live = p.R * sav;
    interp = funbas(fspace,aprime_live(:));
    interp = reshape(interp,NN,p.nxlong*p.nyT) * kron(speye(p.nxlong),income.yTdist);
    if p.Bequests == 1
        interp_death = interp;
    else
        interp_death = sparse(NN,p.nxlong);
        interp_death(:,1) = 1;
    end

    % Transition matrix from 'ii' to 'ii'+1
    T1 = sparse(NN,NN);
    for col = 1:p.nyP*p.nyF*p.nb
        newblock_live = bsxfun(@times,kron(trans_live(:,col),ones(p.nxlong,1)),interp);
        newblock_death = bsxfun(@times,kron(trans_death(:,col),ones(p.nxlong,1)),interp_death);
        T1(:,p.nxlong*(col-1)+1:p.nxlong*col) = (1-p.dieprob)*newblock_live + p.dieprob*newblock_death;
    end
end