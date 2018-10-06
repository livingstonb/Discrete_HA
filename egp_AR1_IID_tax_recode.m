function [sim_results,direct_results,norisk_results] = egp_AR1_IID_tax_recode(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017
    
    % This is the main function file for this code repository. Given a
    % structure of parameters, p, this script calls functions primarily to 
    % compute policy functions via the method of endogenous grip points, 
    % and to find the implied stationary distribution over the state space.

    direct_results.issues = {};
    norisk_results.issues = {};
    sim_results.issues    = {};
    
    %% ADJUST PARAMETERS FOR DATA FREQUENCY

    p.R = 1 + p.r;
    p.R = p.R^(1/p.freq);
    p.r = p.R - 1;
    
    p.beta0         = p.beta0^(1/p.freq);
    p.dieprob       = 1 - (1-p.dieprob)^(1/p.freq);
    p.betaswitch    = 1 - (1-p.betaswitch)^(1/p.freq);
    p.betaL         = p.betaL^(1/p.freq);
    p.betaH         = 1/((p.R)*(1-p.dieprob));
    p.savtax        = p.savtax/p.freq;
    p.Tsim          = p.Tsim * p.freq;
    
    p.N = p.nx*p.nyF*p.nyP*p.nb;
    if p.freq == 1
        direct_results.frequency = 'Annual';
    elseif p.freq == 4
        direct_results.frequency = 'Quarterly';
    else
        error('Frequency must be 1 or 4')
    end
    
    if p.nb > 3
        assert(p.betawidth2>0,'must have betawidth2 > 0')
        assert(p.betawidth2>p.betawidth1,'must have betawidth2 > betawidth1')
    elseif p.nb > 1
        assert(p.betawidth1>0,'must have betawidth1 > 0')
    end

    %% LOAD INCOME VARIABLES2

    income = gen_income_variables(p);
    
    % savtaxthresh should be a multiple of mean gross labor income
    p.savtaxthresh  = p.savtaxthresh * income.meany;

    %% DISCOUNT FACTOR

    % discount factor distribution
    if  p.nb == 1
        prefs.betadist = 1;
        prefs.betatrans = 1;
    elseif p.nb > 1
        prefs.betadist = ones(p.nb,1) / p.nb;
        % nb = 2: prefs.betatrans = [1-p.betaswitch p.betaswitch; p.betaswitch 1-p.betaswitch]; %transitions on average once every 40 years;
         
        betaswitch_ij = p.betaswitch / (p.nb-1);
        % create matrix with (1-betaswitch) on diag and betaswitch_ij
        % elsewhere
        diagonal = (1-p.betaswitch) * ones(p.nb,1);
        off_diag = betaswitch_ij * ones(p.nb);
        off_diag = off_diag - diag(diag(off_diag));
        prefs.betatrans = off_diag + diag(diagonal);
    end
    prefs.betacumdist = cumsum(prefs.betadist);
    prefs.betacumtrans = cumsum(prefs.betatrans,2);
    
    % create grid - add beta to grid later since we may iterate
    bw = p.betawidth;
    switch p.nb
        case 1
            prefs.betagrid0 = 0;
        case 2
            prefs.betagrid0 = [-bw/2 bw/2];
        case 3
            prefs.betagrid0 = [-bw 0 bw];
        case 4
            prefs.betagrid0 = [-3*bw/2 -bw/2 bw/2 3*bw/2];
        case 5
            prefs.betagrid0 = [-2*bw -bw 0 bw 2*bw];
    end
        

    %% ASSET GRIDS

    % savings grid
    sgrid.orig = linspace(0,1,p.nx)';
    sgrid.orig = sgrid.orig.^(1./p.xgrid_par);
    sgrid.orig = p.borrow_lim + (p.xmax-p.borrow_lim).*sgrid.orig;
    sgrid.short = sgrid.orig;
    sgrid.wide = repmat(sgrid.short,[1 p.nyP p.nyF]);
    p.ns = p.nx;

    % xgrid, indexed by beta,yF,yP,x (N by 1 matrix)
    % cash on hand grid: different min points for each value of (iyP,iyF)
    minyT               = kron(min(income.netymat,[],2),ones(p.nx,1));
    xgrid.orig          = sgrid.wide(:) + minyT;
    xgrid.orig_wide     = reshape(xgrid.orig,[p.nx p.nyP p.nyF]);
    
    % for model without income risk
    xgrid.norisk_short  = sgrid.short + income.meannety;
    xgrid.norisk_longgrid   = linspace(0,1,p.nxlong);
    xgrid.norisk_longgrid   = xgrid.norisk_longgrid.^(1/p.xgrid_par);
    xgrid.norisk_longgrid = p.borrow_lim + (p.xmax-p.borrow_lim).*xgrid.norisk_longgrid;
    xgrid.norisk_longgrid = xgrid.norisk_longgrid + income.meannety;
    
    % create longer xgrid
    minyT = kron(min(income.netymat,[],2),ones(p.nxlong,1));
    xgrid.longgrid      = linspace(0,1,p.nxlong)';
    xgrid.longgrid      = xgrid.longgrid.^(1/p.xgrid_par);
    xgrid.longgrid      = p.borrow_lim + (p.xmax - p.borrow_lim)*xgrid.longgrid;
    xgrid.longgrid      = repmat(xgrid.longgrid,p.nyP*p.nyF,1);
    xgrid.longgrid      = xgrid.longgrid + minyT;
    xgrid.longgrid_wide = reshape(xgrid.longgrid,[p.nxlong p.nyP p.nyF]);
    

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
            beta_ub = p.betaH - 1e-5;
        else
            % don't let highest beta be such that (1-dieprob)*R*beta >= 1
            beta_ub = p.betaH - 1e-5 - max(prefs.betagrid0);
        end
        beta_lb = p.betaL;

        check_evals = @(x,y,z) fzero_checkiter(x,y,z,p.maxiterAY);
        options = optimset('TolX',p.tolAY,'OutputFcn',check_evals);
        [beta_final,~,exitflag] = fzero(iterate_EGP,[beta_lb,beta_ub],options);
        if exitflag ~= 1
            direct_results.issues{end+1} = 'NoBetaConv';
            return
        else
        end
    else
        % Beta set in parameters
        beta_final = p.beta0;    
    end
    
    % Get policy functions and stationary distribution for final beta
    [~,basemodel] = solve_EGP(beta_final,p,xgrid,sgrid,prefs,income);

    % Report beta and annualized beta
    direct_results.beta_annualized = beta_final^p.freq;
    direct_results.beta = beta_final;
    
    % Set parameter equal to this beta
    p.beta = beta_final;

    if basemodel.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        direct_results.issues{end+1} = 'NoEGPConv';
        return
    end
    
    %% IMPORTANT MOMENTS
    
    ymat_onlonggrid = repmat(kron(income.ymat,ones(p.nxlong,1)),p.nb,1);
    netymat_onlonggrid = repmat(kron(income.netymat,ones(p.nxlong,1)),p.nb,1);
    
    direct_results.mean_s = basemodel.mean_s;
    direct_results.mean_a = basemodel.mean_a;
    direct_results.mean_x = repmat(xgrid.longgrid(:)',1,p.nb) * basemodel.SSdist;
    direct_results.mean_grossy = (ymat_onlonggrid*income.yTdist)' * basemodel.SSdist;
    direct_results.mean_loggrossy = (log(ymat_onlonggrid)*income.yTdist)' * basemodel.SSdist;
    direct_results.mean_nety = (netymat_onlonggrid*income.yTdist)' * basemodel.SSdist;
    direct_results.mean_lognety = (log(netymat_onlonggrid)*income.yTdist)' * basemodel.SSdist;
    direct_results.var_loggrossy = basemodel.SSdist' * (log(ymat_onlonggrid) - direct_results.mean_loggrossy).^2 * income.yTdist;
    direct_results.var_lognety = basemodel.SSdist' * (log(netymat_onlonggrid)- direct_results.mean_lognety).^2 * income.yTdist;
    
    direct_results.mean_x_check = direct_results.mean_a + direct_results.mean_nety;
   
    % reconstruct yPdist from computed stationary distribution for error
    % checking
    yPdist_check = reshape(basemodel.SSdist_wide,[p.nxlong p.nyP p.nyF*p.nb]);
    yPdist_check = sum(sum(yPdist_check,3),1)';
    

    %% Store problems
    direct_results.issues = {};
    if  abs((p.targetAY - direct_results.mean_a/(income.meany*p.freq))/p.targetAY) > 1e-3
        direct_results.issues{end+1} = 'BadAY';
    end
    if abs((direct_results.mean_x-direct_results.mean_x_check)/direct_results.mean_x)> 1e-3
        direct_results.issues{end+1} = 'DistNotStationary';
    end
    if abs((income.meannety-direct_results.mean_nety)/income.meannety) > 1e-3
        direct_results.issues{end+1} = 'BadNetIncomeMean';
    end
    if norm(yPdist_check-income.yPdist) > 1e-3
        direct_results.issues{end+1} = 'Bad_yP_Dist';
    end
    if min(basemodel.SSdist) < - 1e-3
        direct_results.issues = [direct_results.issues,'NegativeStateProbability'];
    end

    %% WEALTH DISTRIBUTION
    
    % create values for fraction constrained at every pt in asset space,
    % defining constrained as s <= epsilon * mean annual gross labor income 
    % + borrowing limit
    for i = 1:numel(p.epsilon)        
        % create interpolant to find fraction of constrained households
        [sav_unique,ind] = unique(basemodel.sav_longgrid_sort,'last');
        wpinterp = griddedInterpolant(sav_unique,basemodel.SScumdist(ind),'linear');
        if p.epsilon(i) == 0
            % get exact figure
            direct_results.constrained(i) = basemodel.SSdist' * (basemodel.sav_longgrid == 0);
        else
            direct_results.constrained(i) = wpinterp(p.borrow_lim + p.epsilon(i)*income.meany*p.freq);
        end
    end
    
    % wealth percentiles
    if p.WealthInherited == 1
        basemodel.a_longgrid_sort = p.R * basemodel.sav_longgrid_sort(basemodel.SScumdist_uniqueind);
    else
        basemodel.a_longgrid_sort = (1 - p.dieprob) * p.R * basemodel.sav_longgrid_sort(basemodel.SScumdist_uniqueind);
    end
    direct_results.wpercentiles = interp1(basemodel.SScumdist_unique,basemodel.a_longgrid_sort,p.percentiles/100,'linear');
    
    % top shares
    % fraction of total assets that reside in each pt on asset space
    if p.WealthInherited == 1
        totassets = basemodel.SSdist_sort .* (p.R*basemodel.sav_longgrid_sort);
    else
        totassets = basemodel.SSdist_sort .* ((1-p.dieprob)*p.R*basemodel.sav_longgrid_sort);
    end
    cumassets = cumsum(totassets) / direct_results.mean_a;
    cumassets = cumassets(basemodel.SScumdist_uniqueind);
    
    % create interpolant from wealth percentile to cumulative wealth share
    cumwealthshare = griddedInterpolant(basemodel.SScumdist_unique,cumassets,'linear');
    direct_results.top10share  = 1 - cumwealthshare(0.9);
    direct_results.top1share   = 1 - cumwealthshare(0.99);
    
    %% EGP FOR MODEL WITHOUT INCOME RISK
    % Deterministic model
    [norisk,norisk.EGP_cdiff] = solve_EGP_deterministic(p,xgrid,sgrid,prefs,income);
    
    %% SIMULATIONS
    if p.Simulate == 1
        [sim_results,assetmeans] = simulate(p,income,basemodel,xgrid,prefs);
    else
        sim_results = struct();
        assetmeans = [];
    end
    
    %% MPCS
        
    % 1-period MPCs, model with income risk
    for im = 1:numel(p.mpcfrac)
        mpcamount = p.mpcfrac(im) * income.meany * p.freq;
        xmpc = xgrid.longgrid_wide + mpcamount;
        set_mpc_one = false(p.nxlong,p.nyP,p.nyF,p.nb);
        conmpc = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
        for ib  = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            if mpcamount < 0
                below_grid = xmpc(:,iyP,iyF) < xgrid.longgrid_wide(1,iyP,iyF);
                xmpc(below_grid,iyP,iyF) = xgrid.longgrid_wide(1,iyP,iyF);
                set_mpc_one(below_grid,iyP,iyF,ib) = true;
            end
            conmpc(:,iyP,iyF,ib) = basemodel.coninterp{iyP,iyF,ib}(xmpc(:,iyP,iyF));
        end
        end
        end
        mpcs1 = (conmpc(:) - basemodel.con_longgrid) / mpcamount;
        mpcs1(set_mpc_one(:)) = 1;
        direct_results.avg_mpc1{im} = basemodel.SSdist' * mpcs1;
    end

    % 1-period MPCs, norisk
    norisk_results.mpc1 = zeros(1,numel(p.mpcfrac));
    for im = 1:numel(p.mpcfrac)
        mpcamount = p.mpcfrac(im) * income.meany * p.freq;
        xmpc = xgrid.norisk_longgrid + mpcamount;
        set_mpc_one = false(p.nxlong,p.nb);
        conmpc = zeros(p.nxlong,p.nb);
        congrid = zeros(p.nxlong,p.nb);
        for ib  = 1:p.nb
            if mpcamount < 0
                below_grid = xmpc < xgrid.longgrid(1);
                xmpc(below_grid) = xgrid.longgrid(1);
                set_mpc_one(below_grid,ib) = true;
            end
            conmpc(:,ib) = norisk.coninterp{ib}(xmpc);
            congrid(:,ib) = norisk.coninterp{ib}(xgrid.norisk_longgrid);
        end
        mpcs1 = (conmpc(:) - congrid(:)) / mpcamount;
        mpcs1(set_mpc_one(:)) = 1;
        mpcs1 = reshape(mpcs1,p.nxlong,p.nb);
        % Only need mpc at one point
        norisk_results.mpc1(im) = mpcs1(1,:) * prefs.betadist;
    end
    norisk_results.mpc1 = mean(norisk_results.mpc1);
    
    % 4-period MPCs
    if p.freq == 4
        % model with income risk
        if p.ComputeDirectMPC == 1
            % Simulate the full 4 periods
            Tmax = 4;
            
            [mpcs1,mpcs4,avg_mpc1,avg_mpc4,var_mpc1,var_mpc4] = direct_MPCs(p,prefs,income,basemodel,xgrid,Tmax);
            direct_results.avg_mpc1sim = avg_mpc1;
            direct_results.avg_mpc4 = avg_mpc4;
            direct_results.var_mpc1sim = var_mpc1;
            direct_results.var_mpc4 = var_mpc4;
        end

        % norisk model
        [mpc1,mpc4] = direct_MPCs_deterministic(p,prefs,income,norisk);
        norisk_results.mpc1sim  = mpc1;
        norisk_results.mpc4     = mpc4;
    end
    
    % Check that one-period mpc is the same, computed two diff ways
    if abs(norisk_results.mpc1sim - norisk_results.mpc1) / norisk_results.mpc1sim > 1e-3
        norisk_results.issues = 'MPCsInconsistent';
    end
    
    %% GINI
    % Wealth
    if p.WealthInherited == 0
        % Create asset distribution
        dist_live   = (1-p.dieprob) * basemodel.SSdist;
        dist_death  = p.dieprob * basemodel.SSdist;
        level_live  = p.R * basemodel.sav_longgrid;
        level_death = zeros(p.nxlong*p.nyP*p.nyF*p.nb,1);
        distr  = [dist_live;dist_death];
        level = [level_live;level_death];

        direct_results.wealthgini = direct_gini(level,distr);
    else
        distr = basemodel.SSdist_sort;
        level = p.R * basemodel.sav_longgrid_sort;
        direct_results.wealthgini = direct_gini(level,distr);
    end
    
    % Gross income
    direct_results.grossincgini = direct_gini(income.ysortdist,income.ysort);
    
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
    if p.Display == 1
        print_statistics(direct_results,sim_results,norisk_results,p);
    end
    
end