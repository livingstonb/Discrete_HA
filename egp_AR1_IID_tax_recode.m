function [simulations,results] = egp_AR1_IID_tax_recode(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017
    
    % This is the main function file for this code repository. Given a
    % structure of parameters, p, this script calls functions primarily to 
    % compute policy functions via the method of endogenous grip points, 
    % and to find the implied stationary distribution over the state space.
    
    % IMPORTANT STRUCTURES:
    % basemodel - stores objects from the model with income risk
    % norisk - stores objects from the model without income risk
    % simulations - stores simulation results for basemodel
    % sgrid - stores savings grid in various sizes
    % xgrid - stores cash-on-hand grid in various sizes
    % income - stores objects from the income process
    % prefs - stores objects related to preferences
    
    results.issues = {};
    simulations = struct();
    
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

    %% LOAD INCOME

    income = gen_income_variables(p);
    
    % savtaxthresh should be a multiple of mean gross labor income
    p.savtaxthresh  = p.savtaxthresh * income.meany;

    %% DISCOUNT FACTOR

    % discount factor distribution
    if  p.nb == 1
        prefs.betadist = 1;
        prefs.betatrans = 1;
    elseif p.nb ==2 
        prefs.betadist = [0.5;0.5];
        prefs.betatrans = [1-p.betaswitch p.betaswitch; p.betaswitch 1-p.betaswitch]; %transitions on average once every 40 years;
    else
        error('nb must be 1 or 2');
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
    minyT = kron(min(income.netymat,[],2),ones(p.nx,1));
    xgrid.orig = sgrid.wide(:) + minyT;
    xgrid.orig_wide = reshape(xgrid.orig,[p.nx p.nyP p.nyF]);
    xgrid.norisk_short = sgrid.short + income.meannety;
    

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
            results.issues{end+1} = 'NoBetaConv';
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
    
    results.betaannual = beta^p.freq;
    results.betaquarterly = results.betaannual^(1/4);

    % Use final beta to get policy functions and distribution, with a
    % larger grid and higher tolerance for ergodic distribution
    ergodic_tol = 1e-7;
    [~,basemodel,xgrid.longgrid] = solve_EGP(beta,p,xgrid,sgrid,prefs,...
                                            ergodic_tol,income,p.nxlong);
    if basemodel.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        results.issues{end+1} = 'NoEGPConv';
        return
    end
                                        
    xgrid.longgrid_wide = reshape(xgrid.longgrid,[p.nxlong,p.nyP,p.nyF]);
    
    %% IMPORTANT MOMENTS
    
    ymat_onlonggrid = repmat(kron(income.ymat,ones(p.nxlong,1)),p.nb,1);
    netymat_onlonggrid = repmat(kron(income.netymat,ones(p.nxlong,1)),p.nb,1);
    
    results.mean_s = basemodel.sav_longgrid' * basemodel.SSdist;
    results.mean_x = repmat(xgrid.longgrid(:)',1,p.nb) * basemodel.SSdist;
    results.mean_grossy = (ymat_onlonggrid*income.yTdist)' * basemodel.SSdist;
    results.mean_loggrossy = (log(ymat_onlonggrid)*income.yTdist)' * basemodel.SSdist;
    results.mean_nety = (netymat_onlonggrid*income.yTdist)' * basemodel.SSdist;
    results.mean_lognety = (log(netymat_onlonggrid)*income.yTdist)' * basemodel.SSdist;
    results.var_loggrossy = basemodel.SSdist' * (log(ymat_onlonggrid) - results.mean_loggrossy).^2 * income.yTdist;
    results.var_lognety = basemodel.SSdist' * (log(netymat_onlonggrid)- results.mean_lognety).^2 * income.yTdist;
    
    if p.WealthInherited == 1
        results.mean_x_check = p.R*results.mean_s + results.mean_nety;
    else
        results.mean_x_check = p.R*(1-p.dieprob)*results.mean_s + results.mean_nety;
    end
    % reconstruct yPdist from computed stationary distribution for error
    % checking
    yPdist_check = reshape(basemodel.SSdist_wide,[p.nxlong p.nyP p.nyF*p.nb]);
    yPdist_check = sum(sum(yPdist_check,3),1)';
    

    %% Store problems
    results.issues = {};
    if  abs((p.targetAY - results.mean_s/income.meany)/p.targetAY) > 1e-3
        results.issues{end+1} = 'BadAY';
    end
    if abs((results.mean_x-results.mean_x_check)/results.mean_x)> 1e-3
        results.issues{end+1} = 'DistNotStationary';
    end
    if abs((income.meannety-results.mean_nety)/income.meannety) > 1e-3
        results.issues{end+1} = 'BadNetIncomeMean';
    end
    if norm(yPdist_check-income.yPdist) > 1e-3
        results.issues{end+1} = 'Bad_yP_Dist';
    end
    if min(basemodel.SSdist) < - 1e-3
        results.issues = [results.issues,'NegativeStateProbability'];
    end

    %% WEALTH DISTRIBUTION
    
    % create values for fraction constrained at every pt in asset space,
    % defining constrained as a <= epsilon * mean annual gross labor income 
    % + borrowing limit
    for i = 1:numel(p.epsilon)        
        % create interpolant to find fraction of constrained households
        [sav_unique,ind] = unique(basemodel.sav_longgrid_sort,'last');
        wpinterp = griddedInterpolant(sav_unique,basemodel.SScumdist(ind),'linear');
        if p.epsilon(i) == 0
            % get exact figure
            results.constrained(i) = basemodel.SSdist' * (basemodel.sav_longgrid == 0);
        else
            results.constrained(i) = wpinterp(p.borrow_lim + p.epsilon(i)*income.meany*p.freq);
        end
    end
    
    % wealth percentiles
    results.wpercentiles = interp1(basemodel.SScumdist_unique,basemodel.sav_longgrid_sort(basemodel.SScumdist_uniqueind),p.percentiles/100,'linear');
    
    % top shares
    % fraction of total savings that reside in each pt on asset space
    totsav = basemodel.SSdist_sort .* basemodel.sav_longgrid_sort;
    cumsav = cumsum(totsav)/ results.mean_s;
    cumsav = cumsav(basemodel.SScumdist_uniqueind);
    % create interpolant from wealth percentile to cumulative wealth share
    cumwealthshare = griddedInterpolant(basemodel.SScumdist_unique,cumsav,'linear');
    results.top10share = 1 - cumwealthshare(0.9);
    results.top1share = 1 - cumwealthshare(0.99);
    
    %% EGP FOR MODEL WITHOUT INCOME RISK
    % Deterministic model
    if p.SolveDeterministic == 1
        [norisk,norisk.EGP_cdiff] = solve_EGP_deterministic(p,xgrid,sgrid,prefs,income);
    end
    
    %% SIMULATIONS
    % Full model
    if p.Simulate == 1
        simulations = simulate(p,income,basemodel,xgrid,prefs);
    else
        simulations = struct();
    end

    %% MAKE PLOTS
  
    if p.MakePlots ==1 
        makeplots(p,xgrid,sgrid,basemodel,income,simulations);
    end 
    
    %% MPCS
    
    %mpc amounts
    for im = 1:numel(p.mpcfrac)
        mpcamount{im} = p.mpcfrac{im} * income.meany;
        xgrid_mpc{im} = xgrid.longgrid_wide + mpcamount{im};
    end
    
    for im = 1:numel(p.mpcfrac)
        mpc{im} = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
        % iterate over (yP,yF,beta)
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP 
            mpc{im}(:,iyP,iyF,ib) = (basemodel.coninterp{iyP,iyF,ib}(xgrid_mpc{im}(:,iyP,iyF))...
                        - basemodel.con_longgrid_wide(:,iyP,iyF,ib))/mpcamount{im};     
        end
        end
        end
        % average mpc, one-period ahead
        results.avg_mpc1{im} = basemodel.SSdist' * mpc{im}(:);
    end
    
        
    for im = 1:numel(p.mpcfrac)
        results.mpcamount{im} = p.mpcfrac{im};
    end
        
    if p.ComputeDirectMPC == 1
        [results.avg_mpc1_alt,results.avg_mpc4,results.distMPCamount] ...
                = fourperiodmpc(xgrid,p,income,basemodel,prefs);
    end
    
    %% Print Results
    if p.PrintStats == 1
        print_statistics(results,simulations,p);
    end
    
    results.xgrid = xgrid;
end