function [simulations,results] = egp_AR1_IID_tax_recode(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017
    
    % STRUCTURES:
    % basemodel - stores objects from the model with income risk
    % norisk - stores objects from the model without income risk
    % simulations - stores simulation results for basemodel
    % sgrid - stores savings grid in various sizes
    % xgrid - stores cash-on-hand grid in various sizes
    % income - stores objects from the income process
    % prefs - stores objects related to preferences
    
    %% ADJUST PARAMETERS FOR DATA FREQUENCY, OTHER CHOICES

    %p.r = p.r/p.freq;
    p.R = 1 + p.r;
    p.R = p.R^(1/p.freq);
    p.r = p.R - 1;
    
    p.beta0 = p.beta0^(1/p.freq);
    p.dieprob = p.dieprob^(1/p.freq);
    p.betaswitch = p.betaswitch^(1/p.freq);
    p.betaL = p.betaL^(1/p.freq);
    p.betaH = 1/((p.R)*(1-p.dieprob));
    
    if p.Simulate == 0 && p.ComputeSimMPC == 1
        p.ComputeSimMPC = 0;
        disp('SimMPC turned off because Simulate == 0')
    end
    
    p.N = p.nx*p.nyF*p.nyP*p.nb;

    %% LOAD INCOME

    income = gen_income_variables(p);

    %% DISCOUNT FACTOR

    if p.IterateBeta == 0
        p.maxiterAY = 1;
        % final beta
        beta = p.beta0;
        results.beta = beta;
        p.beta = beta;
    end

    %initial discount factor grid
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
    sgrid.norisk_wide = repmat(sgrid.short,[1 p.nb]);
    sgrid.wide = repmat(sgrid.short,[1 p.nyP p.nyF p.nb]);
    p.ns = p.nx;

    % xgrid, indexed by beta,yF,yP,x (N by 1 matrix)
    % cash on hand grid: different min points for each value of (iyP,iyF)
    minyT = repmat(kron(min(income.netymat,[],2),ones(p.nx,1)),p.nb,1);
    xgrid.orig = sgrid.wide(:) + minyT;
    xgrid.orig_wide = reshape(xgrid.orig,[p.nx p.nyP p.nyF p.nb]);
    xgrid.norisk_short = sgrid.short + income.meannety;
    xgrid.norisk_wide = repmat(xgrid.norisk_short,1,p.nb);
    

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
        [beta,exitflag] = iterate_beta(p,xgrid,sgrid,prefs,income);
        if exitflag ~= 1
            results = struct();
            results.issues = {'NoBetaConv'};
            return
        end
        results.beta = beta;
        p.beta = beta;
    end

    % Use final beta to get policy functions and distribution
    ergodic_tol = 1e-7;
    [~,basemodel,xgrid.longgrid] = solve_EGP(beta,p,xgrid,sgrid,prefs,...
                                            ergodic_tol,income,p.nxlong);
                                        
    xgrid.longgrid_wide = reshape(xgrid.longgrid,[p.nxlong,p.nyP,p.nyF,p.nb]);
    
    %% Store important moments
    
    ymat_onlonggrid = repmat(kron(income.ymat,ones(p.nxlong,1)),p.nb,1);
    netymat_onlonggrid = repmat(kron(income.netymat,ones(p.nxlong,1)),p.nb,1);
    
    results.mean_s = basemodel.sav_longgrid' * basemodel.SSdist;
    results.mean_x = xgrid.longgrid(:)' * basemodel.SSdist;
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
    temp = permute(basemodel.SSdist_wide,[2 1 3 4]);
    % reconstruct yPdist from computed stationary distribution for error
    % checking
    yPdist_check = sum(sum(sum(temp,4),3),2);
    

    %% Store problems
    results.issues = {};
    if basemodel.EGP_cdiff > p.tol_iter
        results.issues = [results.issues,'NoEGPConv'];
    end
    if  abs((p.targetAY - results.mean_s)/p.targetAY) > 1e-3
        results.issues = [results.issues,'BadAY'];
    end
    if abs((results.mean_x-results.mean_x_check)/results.mean_x)> 1e-3
        results.issues = [results.issues,'DistNotStationary'];
    end
    if abs((income.meannety-results.mean_nety)/income.meannety) > 1e-3
        results.issues = [results.issues,'BadNetIncomeMean'];
    end
    if norm(yPdist_check-income.yPdist) > 1e-3
        results.issues = [results.issues,'Bad_yP_Dist'];
    end
    if min(basemodel.SSdist) < - 1e-3
        results.issues = [results.issues,'NegativeStateProbability'];
    end

    %% WEALTH DISTRIBUTION
    
    results.frac_constrained = (basemodel.sav_longgrid<=p.borrow_lim)' * basemodel.SSdist;
    results.frac_less5perc_labincome = (basemodel.sav_longgrid<0.05*income.meany)' * basemodel.SSdist;
    % wealth percentiles;
    percentiles = [0.1 0.25 0.5 0.9 0.99];
    [basemodel.SSdist_unique,iu] = unique(basemodel.SScumdist);
    wealthps = interp1(basemodel.SSdist_unique,basemodel.sav_longgrid_sort(iu),percentiles,'linear');
    results.p10wealth = wealthps(1);
    results.p25wealth = wealthps(2);
    results.p50wealth = wealthps(3);
    results.p90wealth = wealthps(4);
    results.p99wealth = wealthps(5);
    
    %% EGP FOR MODEL WITHOUT INCOME RISK
    % Deterministic model
    if p.SolveDeterministic == 1
        [norisk,norisk.EGP_cdiff] = solve_EGP_deterministic(p,xgrid,sgrid,prefs,income);
    end
    
    %% SIMULATIONS
    % Full model
    if p.Simulate == 1
        simulations = simulate(p,income,basemodel,xgrid,prefs,results);
    else
        simulations =[];
    end
    results.simulations = simulations;

    %% MAKE PLOTS
  
    if p.MakePlots ==1 
        makeplots(p,xgrid,sgrid,basemodel,income,simulations);
    end 
    
    %% MPCS
    
    %mpc amounts
    for im = 1:numel(p.mpcfrac)
        mpcamount{im} = p.mpcfrac{im} * income.meany;
        xgrid.mpc{im} = xgrid.longgrid_wide + mpcamount{im};
    end
    
    for im = 1:numel(p.mpcfrac)
        mpc{im} = zeros(p.nxlong,p.nyP,p.nyF,p.nb);
        % iterate over (yP,yF,beta)
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP 
            mpc{im}(:,iyP,iyF,ib) = (basemodel.coninterp{iyP,iyF,ib}(xgrid.mpc{im}(:,iyP,iyF,ib))...
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
        
    if p.ComputeDistMPC == 1
        [results.avg_mpc1_alt,results.avg_mpc4,results.distMPCamount] ...
                = mpc_forward(xgrid,p,income,basemodel,prefs);
    end
    
    %% Print Results
    if p.PrintStats == 1
        print_statistics(results,simulations,p);
    end
    
    results.xgrid = xgrid;
end