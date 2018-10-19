function [avg_mpc1_agrid,mpcs1_a_direct,avg_mpc4_agrid,mpcs4_a_direct,agrid_dist,norisk_mpcs1_a_direct]...
                = direct_MPCs_by_computation(p,basemodel,income,prefs,agrid_short,norisk)

    %% DIRECTLY COMPUTED 1-PERIOD MPCs (MODEL WITH INCOME RISK)

    % Find P(yP,yF,beta|a) = P(a,yP,yF,beta)/P(a)
    Pa = sum(sum(sum(basemodel.adist,4),3),2);
    Pa = repmat(Pa,[1 p.nyP p.nyF p.nb]);
    Pcondl = basemodel.adist ./ Pa;
    Pcondl(Pa == 0) = 0;

    % Each (a,yP,yF) is associated with nyT possible x values, create this
    % grid here
    netymat_reshape = reshape(income.netymat,[1 p.nyP p.nyF p.nyT]);
    netymat_reshape = repmat(netymat_reshape,[p.nxlong 1 1 1]);
    xgrid_yT = repmat(agrid_short,[1 p.nyP p.nyF p.nyT]) + netymat_reshape;

    % For 4-period MPCs
    fspace = fundef({'spli',agrid_short,0,1});
    % transition matrix between (yP,yF,beta) states, cond'l on living
    trans_live = kron(prefs.betatrans,kron(eye(p.nyF),income.yPtrans));
    % transition matrix between (yP,yF,beta) states cond'l on dying
    yPtrans_stationary = repmat(income.yPdist',p.nyP,1);
    trans_death = kron(prefs.betatrans,kron(eye(p.nyF),yPtrans_stationary));
    
    for im = 0:numel(p.mpcfrac)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.mpcfrac(im)*income.meany1*p.freq;
        end

        x_mpc = xgrid_yT + mpcamount;
        sav = zeros(p.nxlong,p.nyP,p.nyF,p.nb,p.nyT);
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            x_iyP_iyF_iyT = x_mpc(:,iyP,iyF,:);
            sav_iyP_iyF_iyT = basemodel.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
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

        if im == 0
            con_baseline = Econ;
        else
            % Compute m(a,yP,yF,beta) = E[m(x,yP,yF,beta)|a,yP,yF,beta]
            mpcs1_a_yP_yF_beta = (Econ - con_baseline) / mpcamount;
            avg_mpc1_agrid(im) = basemodel.adist(:)' * mpcs1_a_yP_yF_beta(:);

            % Compute m(a) = E(m(a,yP,yF,beta)|a)
            %       = sum of P(yP,yF,beta|a) * m(a,yP,yF,beta) over all
            %         (yP,yF,beta)
            mpcs1_a_direct{im} = sum(sum(sum(Pcondl .* mpcs1_a_yP_yF_beta,4),3),2);
        end
        
        % 4-PERIOD MPCS
        % Construct transition matrix from period 1 to period 2
        NN = p.nxlong*p.nyP*p.nyF*p.nb;
        if im>0 && p.freq==4

            aprime_live = p.R * sav;
            interp = funbas(fspace,aprime_live(:));
            interp = reshape(interp,NN,p.nxlong*p.nyT) * kron(speye(p.nxlong),income.yTdist);
            if p.Bequests == 1
                interp_death = interp;
            else
                interp_death = sparse(NN,p.nxlong);
                interp_death(:,1) = 1;
            end
            
            T12 = sparse(NN,NN);
            for col = 1:p.nyP*p.nyF*p.nb
                newblock_live = bsxfun(@times,kron(trans_live(:,col),ones(p.nxlong,1)),interp);
                newblock_death = bsxfun(@times,kron(trans_death(:,col),ones(p.nxlong,1)),interp_death);
                T12(:,p.nxlong*(col-1)+1:p.nxlong*col) = (1-p.dieprob)*newblock_live + p.dieprob*newblock_death;
            end

            mpcs2_a_yP_yF_beta = mpcs1_a_yP_yF_beta(:) + (T12-speye(NN))*con_baseline(:)/mpcamount;
            mpcs3_a_yP_yF_beta = mpcs2_a_yP_yF_beta + (T12*basemodel.statetrans-speye(NN))*con_baseline(:)/mpcamount;
            mpcs4_a_yP_yF_beta = mpcs3_a_yP_yF_beta + (T12*basemodel.statetrans^2-speye(NN))*con_baseline(:)/mpcamount;
            mpcs4_a_yP_yF_beta = reshape(mpcs4_a_yP_yF_beta,[p.nxlong p.nyP p.nyF p.nb]);
            avg_mpc4_agrid(im) = basemodel.adist(:)' * mpcs4_a_yP_yF_beta(:);
            mpcs4_a_direct{im} = sum(sum(sum(Pcondl .* mpcs4_a_yP_yF_beta,4),3),2);
        elseif im>0 && p.freq==1
            avg_mpc4_agrid(im) = NaN;
            mpcs4_a_direct{im} = [];
        end
    end

    % Distribution over agrid, P(a)
    agrid_dist = sum(sum(sum(basemodel.adist,4),3),2);
    
    %% DIRECTLY COMPUTED 1-PERIOD MPCs (MODEL WITHOUT INCOME RISK)
    for im = 0:numel(p.mpcfrac)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.mpcfrac(im)*income.meany1*p.freq;
        end
        
        x_mpc = agrid_short + income.meannety1 + mpcamount;
        con = zeros(p.nxlong,p.nb);
        for ib = 1:p.nb
            con(:,ib) = norisk.coninterp{ib}(x_mpc);
        end
        
        if im == 0
            con_baseline = con;
        else
            % Compute m(a,beta)
            mpcs1_a_beta = (con - con_baseline) / mpcamount;

            % Compute m(x) = E(m(x,beta)|x)
            %       = sum of P(beta|x) * m(x,beta) over all beta
            % beta is exogenous so P(beta|x) = P(beta)
            norisk_mpcs1_a_direct{im} = mpcs1_a_beta * prefs.betadist;
        end
    end
end