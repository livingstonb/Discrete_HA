function [avg_mpc1_agrid,mpcs1_a_direct,agrid_dist,norisk_mpcs1_a_direct]...
            = direct_MPCs_by_computation(p,basemodel,income,prefs,agrid,norisk)

    %% DIRECTLY COMPUTED 1-PERIOD MPCs (MODEL WITH INCOME RISK)
    % First get stationary distribution associated with agrid
    adist = find_stationary_adist(p,basemodel,income,prefs,agrid);

    % Find P(yP,yF,beta|a) = P(a,yP,yF,beta)/P(a)
    Pa = sum(sum(sum(adist,4),3),2);
    Pa = repmat(Pa,[1 p.nyP p.nyF p.nb]);
    Pcondl = adist ./ Pa;
    Pcondl(Pa == 0) = 0;

    % Each (a,yP,yF) is associated with nyT possible x values, create this
    % grid here
    netymat_reshape = reshape(income.netymat,[1 p.nyP p.nyF p.nyT]);
    netymat_reshape = repmat(netymat_reshape,[p.nxlong 1 1 1]);
    xgrid_yT = repmat(agrid,[1 p.nyP p.nyF p.nyT]) + netymat_reshape;

    for im = 0:numel(p.mpcfrac)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.mpcfrac(im)*income.meany1*p.freq;
        end

        x_mpc = xgrid_yT + mpcamount;
        con = zeros(p.nxlong,p.nyP,p.nyF,p.nb,p.nyT);
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            x_iyP_iyF_ib = x_mpc(:,iyP,iyF,:);
            con_iyP_iyF_ib = basemodel.coninterp{iyP,iyF,ib}(x_iyP_iyF_ib(:));
            con(:,iyP,iyF,ib,:) = reshape(con_iyP_iyF_ib,[p.nxlong 1 1 1 p.nyT]);
        end
        end
        end

        % Take expectation over yT
        % con becomes E[con(x,yP,yF,beta)|a,yP,yF,beta]
        con = reshape(con,[],p.nyT) * income.yTdist;
        con = reshape(con,[p.nxlong p.nyP p.nyF p.nb]);

        if im == 0
            con_baseline = con;
        else
            % Compute m(a,yP,yF,beta) = E[m(x,yP,yF,beta)|a,yP,yF,beta]
            mpcs1_a_yP_yF_beta = (con - con_baseline) / mpcamount;
            avg_mpc1_agrid(im) = adist(:)' * mpcs1_a_yP_yF_beta(:);

            % Compute m(a) = E(m(a,yP,yF,beta)|a)
            %       = sum of P(yP,yF,beta|a) * m(a,yP,yF,beta) over all
            %         (yP,yF,beta)
            mpcs1_a_direct{im} = sum(sum(sum(Pcondl .* mpcs1_a_yP_yF_beta,4),3),2);
        end
    end

    % Distribution over agrid, P(a)
    agrid_dist = sum(sum(sum(adist,4),3),2);
    
    %% DIRECTLY COMPUTED 1-PERIOD MPCs (MODEL WITHOUT INCOME RISK)
    for im = 0:numel(p.mpcfrac)
        if im == 0
            mpcamount = 0;
        else
            mpcamount = p.mpcfrac(im)*income.meany1*p.freq;
        end
        
        x_mpc = agrid + income.meannety1 + mpcamount;
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