function [sim_results,direct_results,norisk_results,checks,decomp] ... 
                                                = main(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017
    
    % This is the main function file for this code repository. Given a
    % structure of parameters, p, this script calls functions primarily to 
    % compute policy functions via the method of endogenous grip points, 
    % and to find the implied stationary distribution over the state space.

    direct_results  = struct();
    norisk_results  = struct();
    sim_results     = struct();
    checks          = {};
    
    %% ADJUST PARAMETERS FOR DATA FREQUENCY, OTHER FACTORS

    p.R = 1 + p.r;
    p.R = p.R^(1/p.freq);
    p.r = p.R - 1;
    
    p.savtax        = p.savtax/p.freq;
    p.Tsim          = p.Tsim * p.freq; % Increase simulation time if quarterly
    p.beta0         = p.beta0^(1/p.freq);
    p.dieprob       = 1 - (1-p.dieprob)^(1/p.freq);
    p.betaswitch    = 1 - (1-p.betaswitch)^(1/p.freq);
    p.betaL         = p.betaL^(1/p.freq);
    p.betaH         = 1/((p.R)*(1-p.dieprob));
    
    if p.Annuities == 1
        % Turn off bequests
        p.Bequests = 0;
        p.bequest_weight = 0;
        p.r = p.r + p.dieprob;
    end

    p.N = p.nx*p.nyF*p.nyP*p.nb;
   

    %% LOAD INCOME VARIABLES
    % Create income structure
    income = gen_income_variables(p);
    
    % savtaxthresh should be a multiple of mean gross labor income
    p.savtaxthresh  = p.savtaxthresh * income.meany1;

    %% DISCOUNT FACTOR
    % discount factor distribution
    if  p.nb == 1
        prefs.betadist = 1;
        prefs.betatrans = 1;
    elseif p.nb > 1
        % Equal probability in stationary distribution
        prefs.betadist = ones(p.nb,1) / p.nb; 
        % Probability of switching from beta_i to beta_j, for i=/=j
        betaswitch_ij = p.betaswitch / (p.nb-1);
        % Create matrix with (1-betaswitch) on diag and betaswitch_ij
        % elsewhere
        diagonal = (1-p.betaswitch) * ones(p.nb,1);
        off_diag = betaswitch_ij * ones(p.nb);
        off_diag = off_diag - diag(diag(off_diag));
        prefs.betatrans = off_diag + diag(diagonal);
    end
    prefs.betacumdist = cumsum(prefs.betadist);
    prefs.betacumtrans = cumsum(prefs.betatrans,2);
    
    % Create grid - add beta to grid later since we may iterate
    bw = p.betawidth;
    switch p.nb
        case 1
            prefs.betagrid0 = 0;
        case 2
            prefs.betagrid0 = [-bw/2 bw/2]';
        case 3
            prefs.betagrid0 = [-bw 0 bw]';
        case 4
            prefs.betagrid0 = [-3*bw/2 -bw/2 bw/2 3*bw/2]';
        case 5
            prefs.betagrid0 = [-2*bw -bw 0 bw 2*bw]';
    end
    
    if (max(p.beta0+prefs.betagrid0) >= p.betaH) && (p.IterateBeta == 0)
        error('Max beta on betagrid too high, no stationary distribution')
    end

    %% ASSET GRIDS
    
    % savings grids
    sgrid.orig = linspace(0,1,p.nx)';
    sgrid.orig = sgrid.orig.^(1./p.xgrid_par);
    sgrid.orig = p.borrow_lim + (p.xmax-p.borrow_lim).*sgrid.orig;
    sgrid.short = sgrid.orig;
    sgrid.full = repmat(sgrid.short,[1 p.nyP p.nyF]);
    p.ns = p.nx;

    % xgrids (cash on hand), different min points for each value of (iyP,iyF)
    minyT               = kron(min(income.netymat,[],2),ones(p.nx,1));
    xgrid.orig          = sgrid.full(:) + minyT;
    xgrid.full          = reshape(xgrid.orig,[p.nx p.nyP p.nyF]);
    
    % xgrid for model without income risk
    xgrid.norisk_short  = sgrid.short + income.meannety1;
    xgrid.norisk_longgrid = linspace(0,1,p.nxlong);
    xgrid.norisk_longgrid = xgrid.norisk_longgrid.^(1/p.xgrid_par);
    xgrid.norisk_longgrid = p.borrow_lim + (p.xmax-p.borrow_lim).*xgrid.norisk_longgrid;
    xgrid.norisk_longgrid = xgrid.norisk_longgrid + income.meannety1;
    
    % create longer xgrid
    minyT = kron(min(income.netymat,[],2),ones(p.nxlong,1));
    xgrid.longgrid = linspace(0,1,p.nxlong)';
    xgrid.longgrid = xgrid.longgrid.^(1/p.xgrid_par);
    xgrid.longgrid = p.borrow_lim + (p.xmax - p.borrow_lim)*xgrid.longgrid;
    xgrid.longgrid = repmat(xgrid.longgrid,p.nyP*p.nyF,1);
    xgrid.longgrid = xgrid.longgrid + minyT;
    xgrid.longgrid = reshape(xgrid.longgrid,[p.nxlong p.nyP p.nyF]);
    
    % Create common agrid to compute mpcs along same agrid for all
    % parameterizations. Dimension nxlong x 1
    agrid = linspace(0,1,p.nxlong)';
    agrid = agrid.^(1/p.xgrid_par);
    agrid = p.borrow_lim + (p.xmax - p.borrow_lim) * agrid;
    
    %% UTILITY FUNCTION, BEQUEST FUNCTION
    if p.risk_aver==1
        prefs.u = @(c)log(c);
        prefs.beq = @(a) p.bequest_weight.* log(a+ p.bequest_luxury);
    else    
        prefs.u = @(c)(c.^(1-p.risk_aver)-1)./(1-p.risk_aver);
        prefs.beq = @(a) p.bequest_weight.*((a+p.bequest_luxury).^(1-p.risk_aver)-1)./(1-p.risk_aver);
    end    

    prefs.u1 = @(c) c.^(-p.risk_aver);
    prefs.u1inv = @(u) u.^(-1./p.risk_aver);
    prefs.beq1 = @(a) p.bequest_weight.*(a+p.bequest_luxury).^(-p.risk_aver);

    %% MODEL SOLUTION
    if p.IterateBeta == 1
        
        iterate_EGP = @(x) solve_EGP(x,p,xgrid,sgrid,prefs,income);

        if p.nb == 1
            beta_ub = p.betaH - 1e-4;
        else
            % Don't let highest beta be such that (1-dieprob)*R*beta >= 1
            beta_ub = p.betaH - 1e-4 - max(prefs.betagrid0);
        end
        beta_lb = p.betaL;

        check_evals = @(x,y,z) fzero_checkiter(x,y,z,p.maxiterAY);
        options = optimset('TolX',p.tolAY,'OutputFcn',check_evals);
        [beta_final,~,exitflag] = fzero(iterate_EGP,[beta_lb,beta_ub],options);
        if exitflag ~= 1
            checks{end+1} = 'NoBetaConv';
            return
        end
    else
        % Beta was set in parameters
        beta_final = p.beta0;    
    end
    
    % Get policy functions and stationary distribution for final beta, in
    % 'basemodel' structure
    [~,basemodel] = solve_EGP(beta_final,p,xgrid,sgrid,prefs,income);

    % Report beta and annualized beta
    direct_results.beta_annualized = beta_final^p.freq;
    direct_results.beta = beta_final;
    
    % Set parameter equal to this beta
    p.beta = beta_final;

    if basemodel.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        checks{end+1} = 'NoEGPConv';
        return
    end
    
    %% IMPORTANT MOMENTS
    % Create income grid associated with longgrid
    ymat_onlonggrid = repmat(kron(income.ymat,ones(p.nxlong,1)),p.nb,1);
    netymat_onlonggrid = repmat(kron(income.netymat,ones(p.nxlong,1)),p.nb,1);
    
    direct_results.mean_s = basemodel.mean_s;
    direct_results.mean_a = basemodel.mean_a;
    direct_results.mean_x = repmat(xgrid.longgrid(:)',1,p.nb) * basemodel.SSdist(:);
    if p.Bequests == 1
        direct_results.mean_bequests = p.dieprob * direct_results.mean_s;
    else
        direct_results.mean_bequests = 0;
    end
    
    % One-period income statistics
    direct_results.mean_grossy1 = (ymat_onlonggrid*income.yTdist)' * basemodel.SSdist(:);
    direct_results.mean_loggrossy1 = (log(ymat_onlonggrid)*income.yTdist)' * basemodel.SSdist(:);
    direct_results.mean_nety1 = (netymat_onlonggrid*income.yTdist)' * basemodel.SSdist(:);
    direct_results.mean_lognety1 = (log(netymat_onlonggrid)*income.yTdist)' * basemodel.SSdist(:);
    direct_results.var_loggrossy1 = basemodel.SSdist(:)' * (log(ymat_onlonggrid) - direct_results.mean_loggrossy1).^2 * income.yTdist;
    direct_results.var_lognety1 = basemodel.SSdist(:)' * (log(netymat_onlonggrid)- direct_results.mean_lognety1).^2 * income.yTdist;
    
    direct_results.mean_x_check = direct_results.mean_a + direct_results.mean_nety1;
   
    % Reconstruct yPdist from computed stationary distribution for sanity
    % check
    yPdist_check = reshape(basemodel.SSdist,[p.nxlong p.nyP p.nyF*p.nb]);
    yPdist_check = sum(sum(yPdist_check,3),1)';
    

    %% RECORD PROBLEMS
    if  abs((p.targetAY - direct_results.mean_a/(income.meany1*p.freq))/p.targetAY) > 1e-3
        checks{end+1} = 'BadAY';
    end
    if abs((direct_results.mean_x-direct_results.mean_x_check)/direct_results.mean_x)> 1e-3
        checks{end+1} = 'DistNotStationary';
    end
    if abs((income.meannety1-direct_results.mean_nety1)/income.meannety1) > 1e-3
        checks{end+1} = 'BadNetIncomeMean';
    end
    if norm(yPdist_check-income.yPdist) > 1e-3
        checks{end+1} = 'Bad_yP_Dist';
    end
    if min(basemodel.SSdist(:)) < - 1e-3
        checks{end+1} = 'NegativeStateProbability';
    end

    %% WEALTH DISTRIBUTION
    
    % Create values for fraction constrained at every pt in asset space,
    % defining constrained as s <= epsilon * mean annual gross labor income 
    % + borrowing limit
    [uniquevals, iu] = unique(basemodel.asset_sortvalues,'last');
    wpinterp = griddedInterpolant(uniquevals,basemodel.asset_cumdist(iu),'linear');
    for i = 1:numel(p.epsilon)        
        % create interpolant to find fraction of constrained households
        if p.epsilon(i) == 0
            % Get exact figure
            direct_results.constrained(i) = basemodel.asset_dist(:)' * (basemodel.asset_values(:) == 0);
        else
            direct_results.constrained(i) = wpinterp(p.borrow_lim + p.epsilon(i)*income.meany1*p.freq);
        end
    end
    
    % Wealth percentiles
    wpinterp_inverse = griddedInterpolant(basemodel.asset_cumdist_unique,basemodel.asset_sortvalues(basemodel.asset_uniqueind),'linear');
    direct_results.wpercentiles = wpinterp_inverse(p.percentiles/100);
    
    % Top shares
    % Amount of total assets that reside in each pt on asset space
    totassets = basemodel.asset_dist_sort .* basemodel.asset_sortvalues;
    % Fraction of total assets in each pt on asset space
    cumassets = cumsum(totassets) / direct_results.mean_a;
    cumassets = cumassets(basemodel.asset_uniqueind);
    
    % create interpolant from wealth percentile to cumulative wealth share
    cumwealthshare = griddedInterpolant(basemodel.asset_cumdist_unique,cumassets,'linear');
    direct_results.top10share  = 1 - cumwealthshare(0.9);
    direct_results.top1share   = 1 - cumwealthshare(0.99);
    
    %% EGP FOR MODEL WITHOUT INCOME RISK
    % Deterministic model
    norisk = solve_EGP_deterministic(p,xgrid,sgrid,prefs,income);
    if norisk.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        checks{end+1} = 'NoRiskNoEGPConv';
        return
    end
    
    %% SIMULATIONS
    if p.Simulate == 1
        [sim_results,assetmeans] = simulate(p,income,basemodel,xgrid,prefs);
    else
        assetmeans = [];
    end
    
    %% DIRECTLY COMPUTED 1-PERIOD MPCs (MODEL WITH INCOME RISK)
    % First get stationary distribution associated with agrid
    adist = find_stationary_adist(p,basemodel,income,prefs,agrid);
    
    % Find P(yP,yF,beta|a) = P(a,yP,yF,beta)/P(a)
    Pa = sum(sum(sum(adist,2),3),4);
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
            direct_results.avg_mpc1_agrid(im) = adist(:)' * mpcs1_a_yP_yF_beta(:);
            
            % Compute m(a) = E(m(a,yP,yF,beta)|a)
            %       = sum of P(yP,yF,beta|a) * m(a,yP,yF,beta) over all
            %         (yP,yF,beta)
            direct_results.mpcs1_a_direct{im} = sum(sum(sum(Pcondl .* mpcs1_a_yP_yF_beta,4),3),2);
        end
    end
    
    % Distribution over agrid, P(a)
    direct_results.agrid_dist = sum(sum(sum(adist,4),3),2);
    
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
            norisk_results.mpcs1_a_direct{im} = mpcs1_a_beta * prefs.betadist;
        end
    end
    
    %% MPCs via DRAWING FROM STATIONARY DISTRIBUTION AND SIMULATING
    % Model with income risk
    if p.ComputeDirectMPC == 1          
        [mpcs1,mpcs4,stdev_loggrossy_A,stdev_lognety_A] ...
                            = direct_MPCs(p,prefs,income,basemodel,xgrid);
        basemodel.mpcs1_sim = mpcs1;
        basemodel.mpcs4_sim = mpcs4;
        
        % Find annual mean and standard deviations of income
        if p.freq == 4
            % Direct computations
            direct_results.mean_grossy_A = direct_results.mean_grossy1 * 4;
            % Simulations
            direct_results.stdev_loggrossy_A = stdev_loggrossy_A;
            direct_results.stdev_lognety_A = stdev_lognety_A;     
        else
            % Use direct computations
            direct_results.mean_grossy_A = direct_results.mean_grossy1;
            direct_results.stdev_loggrossy_A = sqrt(direct_results.var_loggrossy1);
            direct_results.stdev_lognety_A = sqrt(direct_results.var_lognety1);
        end
        
        for im = 1:numel(p.mpcfrac)
            direct_results.avg_mpc1_sim(im) = mean(basemodel.mpcs1_sim{im});
            direct_results.var_mpc1_sim(im) = var(basemodel.mpcs1_sim{im});
            if p.freq == 4
                direct_results.avg_mpc4_sim(im) = mean(basemodel.mpcs4_sim{im});
                direct_results.var_mpc4(im) = var(basemodel.mpcs4_sim{im});
            end
        end
    end
    
    
    %% DECOMPOSITION
	decomp = struct([]);
    if p.nb == 1
        m_ra = p.R * (p.beta*p.R)^(-1/p.risk_aver) - 1;
 
        m0 = direct_results.mpcs1_a_direct{5};
        g0 = direct_results.agrid_dist;
        mbc  = norisk_results.mpcs1_a_direct{5};
        for ia = 1:numel(p.abars)
            zidx = agrid < p.abars(ia);
            
            decomp(ia).term1 = m_ra;
            decomp(ia).term2 = (m0(zidx) - m_ra)' * g0(zidx);
            decomp(ia).term3 = (mbc(~zidx) - m_ra)' * g0(~zidx);
            decomp(ia).term4 = (m0(~zidx) - mbc(~zidx))' * g0(~zidx);
        end
    else
        decomp(ia).term1 = NaN;
        decomp(ia).term2 = NaN;
        decomp(ia).term3 = NaN;
        decomp(ia).term4 = NaN;
    end
    
    %% GINI
    % Wealth
    direct_results.wealthgini = direct_gini(basemodel.asset_values,basemodel.asset_dist);
    
    % Gross income
    direct_results.grossincgini = direct_gini(income.ysort,income.ysortdist);
    
    % Net income
    direct_results.netincgini = direct_gini(income.netymat,income.ymatdist);   

    function gini = direct_gini(level,distr)
        % Sort distribution and levels by levels
        sorted = sortrows([level(:),distr(:)]);
        level_sort = sorted(:,1);
        dist_sort  = sorted(:,2);
        S = [0;cumsum(dist_sort .* level_sort)];
        gini = 1 - dist_sort' * (S(1:end-1)+S(2:end)) / S(end);
    end
    
    %% MAKE PLOTS
  
    if p.MakePlots ==1 
        makeplots(p,xgrid,sgrid,basemodel,income,assetmeans);
    end 
    
    %% Print Results
    if p.Display == 1 && 0
        print_statistics(direct_results,sim_results,norisk_results,checks,p,decomp);
    end
    
end