function [sim_results,direct_results] = egp_AR1_IID_tax_recode(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017
    
    % This is the main function file for this code repository. Given a
    % structure of parameters, p, this script calls functions primarily to 
    % compute policy functions via the method of endogenous grip points, 
    % and to find the implied stationary distribution over the state space.

    direct_results.issues = {};
    sim_results = struct();
    
    %% ADJUST PARAMETERS FOR DATA FREQUENCY

    p.R = 1 + p.r;
    p.R = p.R^(1/p.freq);
    p.r = p.R - 1;
    
    p.beta0         = p.beta0^(1/p.freq);
    p.dieprob       = 1 - (1-p.dieprob)^(1/p.freq);
    p.betaswitch    = 1 - (1-p.betaswitch)^(1/p.freq);
    p.betaL         = p.betaL^(1/p.freq);
    p.betaH         = 1/((p.R)*(1-p.dieprob));
    p.targetAY      = p.targetAY/p.freq;
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


    %% LOAD INCOME

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
    
    % create xgrid for intermediate computations of ergodic distribution
    minyT = kron(min(income.netymat,[],2),ones(p.nxinterm,1));
    xgrid.interm      = linspace(0,1,p.nxinterm)';
    xgrid.interm      = xgrid.interm.^(1/p.xgrid_par);
    xgrid.interm      = p.borrow_lim + (p.xmax - p.borrow_lim)*xgrid.interm;
    xgrid.interm      = repmat(xgrid.interm,p.nyP*p.nyF,1);
    xgrid.interm      = xgrid.interm + minyT;
    xgrid.interm_wide = reshape(xgrid.interm,[p.nxinterm p.nyP p.nyF]);    
    
    % create longer xgrid for last computation of ergodic distribution
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
        % Pass to function that will speed up iteration
        [beta,exitflag] = iterate_beta(p,xgrid,sgrid,prefs,income);
        if exitflag ~= 1
            direct_results.issues{end+1} = 'NoBetaConv';
            return
        end
        % beta associated with chosen frequency
        p.beta = beta;
    else
        p.maxiterAY = 1;
        % beta associated with chosen frequency
        p.beta = p.beta0;
        beta = p.beta0;
    end
    
    % Report beta and annualized beta
    direct_results.beta_annualized = beta^p.freq;
    direct_results.beta = beta;

    % Use final beta to get policy functions and distribution, with a
    % larger grid and higher tolerance for ergodic distribution
    ergodic_tol = 1e-8;
    ergodic_method = 1;
    Intermediate = 0;
    [~,basemodel] = solve_EGP(beta,p,xgrid,sgrid,prefs,...
                           	ergodic_method,ergodic_tol,income,Intermediate);
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
    if  abs((p.targetAY - direct_results.mean_a/income.meany)/p.targetAY) > 1e-3
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
    
    %% MPCS (DIRECT METHOD)
        
    % basemodel
    if p.ComputeDirectMPC == 1
        [mpcs1,mpcs4,avg_mpc1,avg_mpc4,var_mpc4] = direct_mpcs(xgrid,p,income,basemodel,prefs);
        direct_results.avg_mpc1 = avg_mpc1;
        direct_results.avg_mpc4 = avg_mpc4;
        direct_results.var_mpc4 = var_mpc4;
    else
        mpcs1 = []; mpcs4 = [];
    end
    
    % norisk model, get mpcs and create norisk.SSdist
    [mpcs1,mpcs4,avg_mpc1,avg_mpc4,norisk] = direct_mpcs_deterministic(...
                                                xgrid,p,income,norisk,prefs);
    
    %% GINI
    % Wealth
    if p.WealthInherited == 0
        dist_live   = (1-p.dieprob) * basemodel.SSdist_sort;
        dist_death  = p.dieprob * basemodel.SSdist_sort;
        level_live  = p.R * basemodel.sav_longgrid_sort;
        level_death = zeros(p.nxlong*p.nyP*p.nyF*p.nb,1);
        
        distr  = [dist_live;dist_death];
        level = [level_live;level_death];
        
        % Re-sort distribution
        temp = sortrows([level distr]);
        level = temp(:,1);
        distr = temp(:,2);
        direct_results.wealthgini = direct_gini(distr,level);
    else
        distr = basemodel.SSdist_sort;
        level = p.R * basemodel.sav_longgrid_sort;
        direct_results.wealthgini = direct_gini(distr,level);
    end
    
    
    % Gross income
    direct_results.grossincgini = direct_gini(income.ysortdist,income.ysort);
    
    % Net income
    temp = sortrows([income.netymat(:) income.ymatdist(:)]);
    direct_results.netincgini = direct_gini(temp(:,2),temp(:,1));   

    function gini = direct_gini(dist_sort,level_sort)
        % Distribution and levels must be sorted in terms of levels
        S = [0;cumsum(dist_sort .* level_sort)];
        gini = 1 - dist_sort' * (S(1:end-1)+S(2:end)) / S(end);
    end
    
    %% MAKE PLOTS
  
    if p.MakePlots ==1 
        makeplots(p,xgrid,sgrid,basemodel,income,sim_results,assetmeans,...
                                                                mpcs1,mpcs4);
    end 
    
    %% Print Results
    if p.Display == 1
        print_statistics(direct_results,sim_results,p);
    end
    
end