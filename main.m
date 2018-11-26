function [results,checks,decomp] = main(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017
    
    % This is the main function file for this code repository. Given a
    % structure of parameters, p, this script calls functions primarily to 
    % compute policy functions via the method of endogenous grip points, 
    % and to find the implied stationary distribution over the state space.

    results = struct('direct',[],'norisk',[],'sim',[]);
    checks = {};
    
    %% --------------------------------------------------------------------
    % LOAD INCOME VARIABLES, CONSTRUCT BETA DISTRIBUTION
    % ---------------------------------------------------------------------

    % Create income structure
    income = gen_income_variables(p);
    
    % savtaxthresh should be a multiple of mean gross labor income
    p.savtaxthresh  = p.savtaxthresh * income.meany1;

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

    %% --------------------------------------------------------------------
    % ASSET GRIDS
    % ---------------------------------------------------------------------
    
    % savings grids
    sgrid.orig = linspace(0,1,p.nx)';
    sgrid.orig = sgrid.orig.^(1./p.xgrid_par);
    sgrid.orig = p.borrow_lim + (p.xmax-p.borrow_lim).*sgrid.orig;
    sgrid.short = sgrid.orig;
    sgrid.full = repmat(sgrid.short,[1 p.nyP p.nyF]);

    % xgrids (cash on hand), different min points for each value of (iyP,iyF)
    minyT               = kron(min(income.netymat,[],2),ones(p.nx,1));
    xgrid.orig          = sgrid.full(:) + minyT;
    xgrid.full          = reshape(xgrid.orig,[p.nx p.nyP p.nyF]);
    
    % xgrid for model without income risk
    xgrid.norisk_short  = sgrid.short + income.meany1;
    xgrid.norisk_longgrid = linspace(0,1,p.nxlong);
    xgrid.norisk_longgrid = xgrid.norisk_longgrid.^(1/p.xgrid_par);
    xgrid.norisk_longgrid = p.borrow_lim + (p.xmax-p.borrow_lim).*xgrid.norisk_longgrid;
    xgrid.norisk_longgrid = xgrid.norisk_longgrid + income.meany1;
    
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
    % Decomp2 only valid between specs with same nxlong, xgrid_par,
    % borrow_lim, and xmax
    agrid = linspace(0,1,p.nxlong)';
    agrid = agrid.^(1/p.xgrid_par);
    agrid = p.borrow_lim + (p.xmax - p.borrow_lim) * agrid;
    agrid_short = agrid;
    agrid = repmat(agrid,p.nyP*p.nyF*p.nb,1);
    
    %% --------------------------------------------------------------------
    % UTILITY FUNCTION, BEQUEST FUNCTION
    % ---------------------------------------------------------------------
    if p.risk_aver==1
        prefs.u = @(c)log(c);
    else    
        prefs.u = @(c)(c.^(1-p.risk_aver)-1)./(1-p.risk_aver);
        
    end    
    
    if p.bequest_curv == 1
        prefs.beq = @(a) p.bequest_weight.* log(a+ p.bequest_luxury);
    else
        prefs.beq = @(a) p.bequest_weight.*((a+p.bequest_luxury).^(1-p.bequest_curv)-1)./(1-p.bequest_curv);
    end

    prefs.u1 = @(c) c.^(-p.risk_aver);
    prefs.u1inv = @(u) u.^(-1./p.risk_aver);
    prefs.beq1 = @(a) p.bequest_weight.*(a+p.bequest_luxury).^(-p.bequest_curv);

    %% --------------------------------------------------------------------
    % MODEL SOLUTION
    % ---------------------------------------------------------------------

    if p.IterateBeta == 1
        
        mpcshock = 0;
        Iterating = 1;
        if p.EpsteinZin == 1
            iterate_EGP = @(x) solve_EGP_EZ(x,p,xgrid,sgrid,agrid_short,prefs,income,Iterating);
        else
            iterate_EGP = @(x) solve_EGP(x,p,xgrid,sgrid,agrid_short,prefs,income,Iterating,mpcshock,[]);
        end

        if p.nb == 1
            beta_ub = p.betaH;
        else
            % Don't let highest beta be such that (1-dieprob)*R*beta >= 1
            beta_ub = p.betaH  - max(prefs.betagrid0);
        end
        beta_lb = p.betaL;

        % output function that limits number of fzero iterations
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
    Iterating = 0;
    if p.EpsteinZin == 1
        [~,basemodel] = solve_EGP_EZ(beta_final,p,xgrid,sgrid,agrid_short,prefs,income,Iterating);
    else
        mpcshock = 0;
        [~,basemodel] = solve_EGP(beta_final,p,xgrid,sgrid,agrid_short,prefs,income,Iterating,mpcshock,[]);
    end
    results.direct.adist = basemodel.adist;

    % Report beta and annualized beta
    results.direct.beta_annualized = beta_final^p.freq;
    results.direct.beta = beta_final;
    
    if basemodel.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        checks{end+1} = 'NoEGPConv';
        return
    end
    
    %% --------------------------------------------------------------------
    % IMPORTANT MOMENTS
    % ---------------------------------------------------------------------

    results.direct.mean_s = basemodel.xdist(:)' * basemodel.sav_x(:);
    results.direct.mean_a = basemodel.mean_a;
    results.direct.mean_x = basemodel.xdist(:)' * basemodel.xvals(:);
    results.direct.mean_c = basemodel.xdist(:)' * basemodel.con_x(:);
    
    % One-period income statistics
    results.direct.mean_grossy1 = basemodel.xdist(:)' * basemodel.y_x(:);
    results.direct.mean_loggrossy1 = basemodel.xdist(:)' * log(basemodel.y_x(:));
    results.direct.mean_nety1 = basemodel.xdist(:)' * basemodel.nety_x(:);
    results.direct.mean_lognety1 = basemodel.xdist(:)' * log(basemodel.nety_x(:));
    results.direct.var_loggrossy1 = basemodel.xdist(:)' * (log(basemodel.y_x(:)) - results.direct.mean_loggrossy1).^2;
    results.direct.var_lognety1 = basemodel.xdist(:)' * (log(basemodel.nety_x(:)) - results.direct.mean_lognety1).^2;
    
    results.direct.mean_x_check = results.direct.mean_a + results.direct.mean_nety1;
   
    % Reconstruct yPdist and yFdist from computed stationary distribution 
    % for sanity check
    yPdist_check = reshape(basemodel.adist,[p.nxlong p.nyP p.nyF*p.nb]);
    yPdist_check = sum(sum(yPdist_check,3),1)';
    yFdist_check = reshape(basemodel.adist,[p.nxlong*p.nyP p.nyF p.nb]);
    yFdist_check = sum(sum(yFdist_check,3),1)';
    
    %% --------------------------------------------------------------------
    % RECORD PROBLEMS
    % ---------------------------------------------------------------------
    if  abs((p.targetAY - results.direct.mean_a/(income.meany1*p.freq))/p.targetAY) > 1e-3
        checks{end+1} = 'BadAY';
    end
    if abs((results.direct.mean_x-results.direct.mean_x_check)/results.direct.mean_x)> 1e-3
        checks{end+1} = 'DistNotStationary';
    end
    if abs((income.meannety1-results.direct.mean_nety1)/income.meannety1) > 1e-3
        checks{end+1} = 'BadNetIncomeMean';
    end
    if norm(yPdist_check-income.yPdist) > 1e-3
        checks{end+1} = 'Bad_yP_Dist';
    end
    if norm(yFdist_check-income.yFdist) > 1e-3
        checks{end+1} = 'Bad_yF_Dist';
    end
    if basemodel.adiff > 1e-8
        checks{end+1} = sprintf('NoStatConv, adiff = %5.3e',basemodel.adiff);
    end
    if min(basemodel.adist(:)) < - 1e-3
        checks{end+1} = 'LargeNegativeStateProbability';
    elseif min(basemodel.adist(:)) < -1e-8
        checks{end+1} = 'MedNegativeStateProbability';
    elseif min(basemodel.adist(:)) < -1e-13
        checks{end+1} = 'SmallNegativeStateProbability';
    end

    %% --------------------------------------------------------------------
    % WEALTH DISTRIBUTION
    % ---------------------------------------------------------------------

    % Create values for fraction constrained at every pt in asset space,
    % defining constrained as s <= epsilon * mean annual gross labor income 
    % + borrowing limit
    sort_aspace = sortrows([agrid basemodel.adist(:)]);
    sort_agrid = sort_aspace(:,1);
    sort_adist = sort_aspace(:,2);

    sort_acumdist = cumsum(sort_adist);

    [aunique,uind] = unique(sort_agrid,'last');
    wpinterp = griddedInterpolant(aunique,sort_acumdist(uind),'linear');
    for i = 1:numel(p.epsilon)        
        % create interpolant to find fraction of constrained households
        if p.epsilon(i) == 0
            % Get exact figure
            results.direct.constrained(i) = basemodel.adist(:)' * (agrid==0);
        else
            results.direct.constrained(i) = wpinterp(p.borrow_lim + p.epsilon(i)*income.meany1*p.freq);
        end
    end
    
    % Wealth percentiles
    [acumdist_unique,uniqueind] = unique(sort_acumdist,'last');
    wpinterp_inverse = griddedInterpolant(acumdist_unique,sort_agrid(uniqueind),'linear');
    results.direct.wpercentiles = wpinterp_inverse(p.percentiles/100);
    
    % Top shares
    % Amount of total assets that reside in each pt on sorted asset space
    totassets = sort_adist .* sort_agrid;
    % Fraction of total assets in each pt on asset space
    cumassets = cumsum(totassets) / results.direct.mean_a;
    
    % create interpolant from wealth percentile to cumulative wealth share
    cumwealthshare = griddedInterpolant(acumdist_unique,cumassets(uniqueind),'linear');
    results.direct.top10share  = 1 - cumwealthshare(0.9);
    results.direct.top1share   = 1 - cumwealthshare(0.99);
    
    %% --------------------------------------------------------------------
    % EGP FOR MODEL WITHOUT INCOME RISK
    % ---------------------------------------------------------------------
    
    % Deterministic model
    norisk = solve_EGP_deterministic(p,xgrid,sgrid,prefs,income,results.direct);
    if norisk.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        checks{end+1} = 'NoRiskNoEGPConv';
        return
    end
    
    %% --------------------------------------------------------------------
    % SIMULATIONS
    % ---------------------------------------------------------------------
    if p.Simulate == 1
        [results.sim,assetmeans] = simulate(p,income,basemodel,xgrid,prefs);
    else
        assetmeans = [];
    end
    
    %% --------------------------------------------------------------------
    % DIRECTLY COMPUTED MPCs, IMPC(1,j)
    % ---------------------------------------------------------------------
    [MPCs,agrid_dist] = direct_MPCs_by_computation(p,basemodel,income,prefs,agrid_short);
    results.direct.mpcs = MPCs;
    results.direct.agrid_dist = agrid_dist;

    results.norisk.mpcs1_a_direct = direct_MPCs_by_computation_norisk(p,norisk,income,prefs,agrid_short);

    %% --------------------------------------------------------------------
    % DIRECTLY COMPUTED MPCs, IMPC(s,t)
    % ---------------------------------------------------------------------
    
    % mpcmodels(s,t) stores the policy functions associated with the case
    % where the household is currently in period t, but recieved news about
    % the period-s shock in period 1
    mpcmodels = cell(4,4);
    
    % policy functions are the same as baseline when shock is received in
    % the current period
    mpcmodels{1,1} = basemodel;
    mpcmodels{2,2} = basemodel;
    mpcmodels{3,3} = basemodel;
    mpcmodels{4,4} = basemodel;
    
    
    for is = 2:4
    for it = is-1:-1:1
        
        if it == is-1
        	nextmpcshock = 0.01 * income.meany1 * p.freq;
        else
            nextmpcshock = 0;
        end
        Iterating = 0;
        [~,mpcmodels{is,it}] = solve_EGP(results.direct.beta,p,xgrid,sgrid,...
                    agrid_short,prefs,income,Iterating,nextmpcshock,mpcmodels{is,it+1});
    end
    end
    
    [MPCs,agrid_dist] = direct_MPCs_by_computation_new(p,basemodel,mpcmodels,income,prefs,agrid_short);
    
    
    %% --------------------------------------------------------------------
    % MPCs via DRAWING FROM STATIONARY DISTRIBUTION AND SIMULATING
    % ---------------------------------------------------------------------
    % Model with income risk
    if p.freq == 4
        [MPCs,stdev_loggrossy_A,stdev_lognety_A] ...
                            = direct_MPCs_by_simulation(p,prefs,income,basemodel,xgrid,agrid);
        results.direct.mpcs_sim = MPCs;
    end

    % Find annual mean and standard deviations of income
    if p.freq == 4
        % Direct computations
        results.direct.mean_grossy_A = results.direct.mean_grossy1 * 4;
        % Simulations
        results.direct.stdev_loggrossy_A = stdev_loggrossy_A;
        results.direct.stdev_lognety_A = stdev_lognety_A;     
    else
        % Use direct computations
        results.direct.mean_grossy_A = results.direct.mean_grossy1;
        results.direct.stdev_loggrossy_A = sqrt(results.direct.var_loggrossy1);
        results.direct.stdev_lognety_A = sqrt(results.direct.var_lognety1);
    end

    %% --------------------------------------------------------------------
    % DECOMPOSITION 1 (DECOMP OF EM)
    % ---------------------------------------------------------------------
	decomp = struct([]);
    if p.nb == 1 && p.EpsteinZin == 0 && p.bequest_weight == 0 && p.temptation == 0
        m_ra = p.R * (results.direct.beta*p.R)^(-1/p.risk_aver) - 1;
 
        % MPC shock of 0.01 * annual income
        m0 = results.direct.mpcs.mpcs_1_1{5};
        g0 = results.direct.agrid_dist;
        mbc  = results.norisk.mpcs1_a_direct{5};
        for ia = 1:numel(p.abars)
            zidx = agrid_short <= p.abars(ia);
            
            decomp(ia).term1 = m_ra;
            decomp(ia).term2 = (m0(zidx) - m_ra)' * g0(zidx);
            decomp(ia).term3 = (mbc(~zidx) - m_ra)' * g0(~zidx);
            decomp(ia).term4 = (m0(~zidx) - mbc(~zidx))' * g0(~zidx);
        end
    else
        for ia = 1:numel(p.abars)
            decomp(ia).term1 = NaN;
            decomp(ia).term2 = NaN;
            decomp(ia).term3 = NaN;
            decomp(ia).term4 = NaN;
        end
    end
    
    %% --------------------------------------------------------------------
    % GINI
    % ---------------------------------------------------------------------
    % Wealth
    results.direct.wealthgini = direct_gini(agrid,basemodel.adist);
    
    % Gross income
    results.direct.grossincgini = direct_gini(income.ysort,income.ysortdist);
    
    % Net income
    results.direct.netincgini = direct_gini(income.netymat,income.ymatdist);   

    function gini = direct_gini(level,distr)
        % Sort distribution and levels by levels
        sorted = sortrows([level(:),distr(:)]);
        level_sort = sorted(:,1);
        dist_sort  = sorted(:,2);
        S = [0;cumsum(dist_sort .* level_sort)];
        gini = 1 - dist_sort' * (S(1:end-1)+S(2:end)) / S(end);
    end
    
end