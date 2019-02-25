function [MPCs,agrid_dist] = direct_MPCs_by_computation(p,basemodel,models,income,prefs,xgrid,agrid_short,shocksize)
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
    
    % transition matrix between (yP,yF,beta) states cond'l on dying
    yPtrans_stationary = repmat(income.yPdist',p.nyP,1);
    
    % transition matrix between (yP,yF,beta) states, cond'l on living
    if (numel(p.risk_aver) == 1) && (numel(p.invies) == 1)
        trans_live = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
        trans_death = kron(prefs.betatrans,kron(eye(p.nyF),yPtrans_stationary));
    else
        trans_live = kron(prefs.IEStrans,kron(eye(p.nyF),income.yPtrans));
        trans_death = kron(prefs.IEStrans,kron(eye(p.nyF),yPtrans_stationary));
    end
    
    con_baseline_yT = get_policy(p,xgrid_yT,basemodel,income);
    con_baseline = reshape(con_baseline_yT,[],p.nyT) * income.yTdist;
    
    mpcamount = shocksize;

    %% --------------------------------------------------------------------
    % ITERATION OVER (s,t)
    % ---------------------------------------------------------------------
    % period t MPC out of shock in period s, learned about in period 1
    MPCs.avg_s_t = cell(16,16);
    MPCs.mpcs_1_t = cell(1,4);
    maxT = p.freq * 4;
    
    if (shocksize < 0) || (p.mpcshocks_after_period1 == 1)
        IS = 1;
    elseif maxT == 4
        IS = [1 2 3 4];
    else
        IS = [1 2 3 4 5 9 13];
    end
    
    for is = IS
        % iterate over t within s

        if (p.EpsteinZin == 1) && (is > 1)
            continue
        end

        for it = 1:is % in this block, look at t <= s
            if p.Display == 1
                fprintf('    Processing (s,t) = (%d,%d)\n',is,it)
                disp(['    --Time ' datestr(now,'HH:MM:SS')])
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
                    models,is,it-1,fspace,trans_live,trans_death,mpcshock);
                T1t = T1t * Ti;
            end

%             for ii = 1:it
%                 if ii == 1
%                     T1t = speye(NN);
%                 else
%                     if p.Display == 1
%                         fprintf('      Transition from t=%d to t=%d\n',ii-1,ii)
%                         disp(['      --Time ' datestr(now,'HH:MM:SS')])
%                     end
%                     %Ti is transition from t=ii-1 to t=ii
%                     mpcshock = 0;
%                     Ti = transition_t_less_s(p,income,xgrid_yT,...
%                         models,is,ii-1,fspace,trans_live,trans_death,mpcshock);
%                     T1t = T1t * Ti;
%                 end
%             end

            % get consumption policy function
            if it == is
                % shock is 1% of annual income
                x_mpc = xgrid_yT + shocksize;
            else
                x_mpc = xgrid_yT;
            end

            if shocksize < 0 && it == 1
                below_xgrid = false(size(x_mpc));
                for iyT = 1:p.nyT
                    below_xgrid (:,:,:,iyT) = x_mpc(:,:,:,iyT) < xgrid.full(1,:,:);
                end
                below_xgrid = reshape(below_xgrid,[p.nxlong p.nyP p.nyF 1 p.nyT]);
                below_xgrid = repmat(below_xgrid,[1 1 1 p.nb 1]);
            end

            con = get_policy(p,x_mpc,models{is,it},income);

            if shocksize < 0 && it == 1
                % set MPC 1 for points below xgrid
                con(below_xgrid) = con_baseline_yT(below_xgrid) + mpcamount;
            elseif shocksize < 0 && it > 1
                % set MPC 0 for points that were already set to 1
                con(below_xgrid) = con_baseline_yT(below_xgrid);
            end

            con = reshape(con,[],p.nyT) * income.yTdist;

            % now compute IMPC(s,t)
            mpcs = ( T1t * con - con_baseline) / mpcamount;


            MPCs.avg_s_t{is,it} = basemodel.adist(:)' * mpcs(:);
            
            if (is == 1) && (it >= 1 && it <= 4)
                MPCs.mpcs_1_t{it} = mpcs;
            end

        end

        % transition probabilities from it = is to it = is + 1
        mpcshock = mpcamount;
        T_s1_s = transition_t_less_s(p,income,xgrid_yT,models,is,is,...
                                    fspace,trans_live,trans_death,mpcshock);
        T1t = T1t * T_s1_s;
        clear T_s1_s

        RHScon = con_baseline(:);
        for it = is+1:maxT % it > is case, policy fcns stay the same in this region
            mpcs = ( T1t * RHScon(:) - con_baseline(:) ) / mpcamount;
            MPCs.avg_s_t{is,it} = basemodel.adist(:)' * mpcs(:);
            RHScon = basemodel.statetrans * RHScon;
            
            if (is == 1) && (it >= 1 && it <= 4)
                MPCs.mpcs_1_t{it} = mpcs;
            end
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
    else
        for is = [6 7 8 10 11 12 14 15 16]
        for it = 1:16
            MPCs.avg_s_t{is,it} = NaN;
        end
        end
    end

    if p.EpsteinZin == 1
        for is = 2:16
        for it = 1:16
            MPCs.avg_s_t{is,it} = NaN;
        end
        end
    end

    %% --------------------------------------------------------------------
    % FIND CUMULATIVE MPCs 1-4, 5-8, 9-12, 13-16
    % ---------------------------------------------------------------------
    % shock in period 1
    MPCs.avg_1_1to4 = cellsum(MPCs.avg_s_t,1,1,4);
    MPCs.avg_1_5to8 = cellsum(MPCs.avg_s_t,1,5,8);
    MPCs.avg_1_9to12 = cellsum(MPCs.avg_s_t,1,9,12);
    MPCs.avg_1_13to16 = cellsum(MPCs.avg_s_t,1,13,16);

    % shock in period 2
    MPCs.avg_2_1to4 = cellsum(MPCs.avg_s_t,2,1,4);
    MPCs.avg_2_5to8 = cellsum(MPCs.avg_s_t,2,5,8);
    MPCs.avg_2_9to12 = cellsum(MPCs.avg_s_t,2,9,12);
    MPCs.avg_2_13to16 = cellsum(MPCs.avg_s_t,2,13,16);

    % shock in period 3
    MPCs.avg_3_1to4 = cellsum(MPCs.avg_s_t,3,1,4);
    MPCs.avg_3_5to8 = cellsum(MPCs.avg_s_t,3,5,8);
    MPCs.avg_3_9to12 = cellsum(MPCs.avg_s_t,3,9,12);
    MPCs.avg_3_13to16 = cellsum(MPCs.avg_s_t,3,13,16);

    % shock in period 4
    MPCs.avg_4_1to4 = cellsum(MPCs.avg_s_t,4,1,4);
    MPCs.avg_4_5to8 = cellsum(MPCs.avg_s_t,4,5,8);
    MPCs.avg_4_9to12 = cellsum(MPCs.avg_s_t,4,9,12);
    MPCs.avg_4_13to16 = cellsum(MPCs.avg_s_t,4,13,16);

    % shock in period 5
    MPCs.avg_5_1to4 = cellsum(MPCs.avg_s_t,5,1,4);
    MPCs.avg_5_5to8 = cellsum(MPCs.avg_s_t,5,5,8);
    MPCs.avg_5_9to12 = cellsum(MPCs.avg_s_t,5,9,12);
    MPCs.avg_5_13to16 = cellsum(MPCs.avg_s_t,5,13,16);

    % shock in period 9
    MPCs.avg_9_1to4 = cellsum(MPCs.avg_s_t,9,1,4);
    MPCs.avg_9_5to8 = cellsum(MPCs.avg_s_t,9,5,8);
    MPCs.avg_9_9to12 = cellsum(MPCs.avg_s_t,9,9,12);
    MPCs.avg_9_13to16 = cellsum(MPCs.avg_s_t,9,13,16);

    % shock in period 13
    MPCs.avg_13_1to4 = cellsum(MPCs.avg_s_t,13,1,4);
    MPCs.avg_13_5to8 = cellsum(MPCs.avg_s_t,13,5,8);
    MPCs.avg_13_9to12 = cellsum(MPCs.avg_s_t,13,9,12);
    MPCs.avg_13_13to16 = cellsum(MPCs.avg_s_t,13,13,16);

    %% --------------------------------------------------------------------
    % DISTRIBUTION OVER AGRID, P(a)
    % ---------------------------------------------------------------------
    agrid_dist = sum(sum(sum(basemodel.adist,4),3),2);
end

%% ------------------------------------------------------------------------
% SUBFUNCTIONS
% -------------------------------------------------------------------------

function con = get_policy(p,x_mpc,model,income)
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
    x_mpc = reshape(x_mpc,[p.nxlong p.nyP p.nyF 1 p.nyT]);
    con = repmat(x_mpc,[1 1 1 p.nb 1]) - sav - p.savtax * max(sav-p.savtaxthresh,0);

    % Take expectation over yT
    % con becomes E[con(x,yP,yF,beta)|a,yP,yF,beta]
    % Esav = reshape(sav,[],p.nyT) * income.yTdist;
    % Esav = reshape(Esav,[p.nxlong p.nyP p.nyF p.nb]);
    % Exx = reshape(reshape(x_mpc,[],p.nyT) * income.yTdist,[p.nxlong p.nyP p.nyF]);
    % Econ = repmat(Exx,[1 1 1 p.nb]) - Esav - p.savtax * max(Esav - p.savtaxthresh,0);  
end

function T1 = transition_t_less_s(p,income,xgrid_yT,models,is,ii,...
                                            fspace,trans_live,trans_death,mpcshock)
    % Computes the transition matrix between t=ii and t=ii + 1 given shock
    % in is
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

function out = cellsum(array,index,start,last)

    out = 0;
    for i = start:last
        out = out + array{index,i};
    end
end

