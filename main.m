function [sim_results,direct_results,norisk_results,checks,decomp] ... 
                                                = egp_AR1_IID_tax_recode(p)
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
    xgrid.longgrid      = linspace(0,1,p.nxlong)';
    xgrid.longgrid      = xgrid.longgrid.^(1/p.xgrid_par);
    xgrid.longgrid      = p.borrow_lim + (p.xmax - p.borrow_lim)*xgrid.longgrid;
    xgrid.longgrid      = repmat(xgrid.longgrid,p.nyP*p.nyF,1);
    xgrid.longgrid      = xgrid.longgrid + minyT;
    xgrid.longgrid      = reshape(xgrid.longgrid,[p.nxlong p.nyP p.nyF]);
    
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

    %% MPCs FROM DRAWING FROM STATIONARY DISTRIBUTION AND SIMULATING
    % Model with income risk
    if p.ComputeDirectMPC == 1          
        [a1,betaindsim0,mpcs1,mpcs4,stdev_loggrossy_A,stdev_lognety_A,mean_grossy_A] ...
                        = direct_MPCs(p,prefs,income,basemodel,xgrid);
        basemodel.a1 = a1;
        basemodel.betaindsim0 = betaindsim0;
        basemodel.mpcs1 = mpcs1;
        basemodel.mpcs4 = mpcs4;
        
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
            direct_results.avg_mpc1(im) = mean(basemodel.mpcs1{im});
            direct_results.var_mpc1(im) = var(basemodel.mpcs1{im});
            if p.freq == 4
                direct_results.avg_mpc4(im) = mean(basemodel.mpcs4{im});
                direct_results.var_mpc4(im) = var(basemodel.mpcs4{im});
            end
        end
    end

    % norisk model
    % Initial assets distributed as in basemodel.a1, so the mpc
    % distribution here will not match that of norisk_mpcs1, which is based
    % off stationary distribution of cash-on-hand from risky model, not
    % assets
    [norisk.mpcs1,norisk.mpcs4] = ...
                direct_MPCs_deterministic(p,prefs,income,norisk,basemodel);
    norisk_results.avg_mpc1 = mean(norisk.mpcs1);
    norisk_results.avg_mpc4 = mean(norisk.mpcs4);
    
    %% DECOMPOSITION
    decomp = struct([]);
    if p.nb == 1
        m_ra = p.R * (p.beta*p.R)^(-1/p.risk_aver) - 1;
        for ia = 1:numel(p.abars)
            % Use the initial distribution of assets from direct_MPCs(), a1
            cind  = basemodel.a1 <= p.abars(ia); % constrained households
            Ntotal = numel(basemodel.a1);
            decomp(ia).term1 = m_ra;
            decomp(ia).term2 = sum( (basemodel.mpcs1{5} - m_ra) .* cind/Ntotal);
            decomp(ia).term3 = sum( (norisk.mpcs1 - m_ra) .* (~cind)/Ntotal);
            decomp(ia).term4 = sum( (basemodel.mpcs1{5} - norisk.mpcs1) .* (~cind)/Ntotal);
        end
    else
        for ia = 1:numel(p.abars)
            decomp(ia).term1 = NaN;
            decomp(ia).term2 = NaN;
            decomp(ia).term3 = NaN;
            decomp(ia).term4 = NaN;
        end
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