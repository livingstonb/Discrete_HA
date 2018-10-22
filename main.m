function [income,sim_results,direct_results,norisk_results,checks,decomp] ... 
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

    % Deal with special cases
    if p.freq == 1
        add_1eneg2 = {'2 RandomBetaHet5 Width0.01 SwitchProb0.1 NoDeath'
                        '2 RandomBetaHet5 Width0.01 SwitchProb0.1 Death'};
        sub_1eneg2 = {'2 BeqWt0.02 BeqLux0.01 BeqCurv0.1'};
        if ismember(p.name,add_1eneg2)
            p.betaH = p.betaH + 1e-2;
        elseif ismember(p.name,sub_1eneg2);
            p.betaH = p.betaH - 1e-2;
        else
            p.betaH = p.betaH - 1e-3;
        end
    elseif p.freq == 4
        add_1eneg3 = {'2 RandomBetaHet5 Width0.005 SwitchProb0.02 NoDeath'};
        add_5eneg3 = {'2 RandomBetaHet5 Width0.005 SwitchProb0.1 NoDeath'
                        '2 RandomBetaHet5 Width0.005 SwitchProb0.1 Death'
                        '2 RandomBetaHet5 Width0.01 SwitchProb0.02 NoDeath'
                        '2 RandomBetaHet5 Width0.01 SwitchProb0.02 Death'};
        add_125eneg2 = {'2 RandomBetaHet5 Width0.01 SwitchProb0.1 NoDeath'
                        '2 RandomBetaHet5 Width0.01 SwitchProb0.1 Death'};
        sub_1eneg5 = {'4 Temptation0.05'};
        if ismember(p.name,add_1eneg3);
            p.betaH = p.betaH + 1e-3;
        elseif ismember(p.name,add_5eneg3);
            p.betaH = p.betaH + 5e-3;
        elseif ismember(p.name,add_125eneg2);
            p.betaH = p.betaH + 1.25e-2;
        elseif ismember(p.name,sub_1eneg5);
            p.betaH = p.betaH - 1e-5;
        else
            p.betaH = p.betaH - 1e-3;
        end
    end
    
    if p.Annuities == 1
        % Turn off bequests
        p.Bequests = 0;
        p.bequest_weight = 0;
        p.r = p.r + p.dieprob;
        p.R = 1 + p.r;
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
    agrid_short = agrid;
    agrid = repmat(agrid,p.nyP*p.nyF*p.nb,1);
    
    %% UTILITY FUNCTION, BEQUEST FUNCTION
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

    %% MODEL SOLUTION
    if p.IterateBeta == 1
        
        Iterating = 1;
        if p.EpsteinZin == 1
            iterate_EGP = @(x) solve_EGP_EZ(x,p,xgrid,sgrid,agrid_short,prefs,income,Iterating);
        else
            iterate_EGP = @(x) solve_EGP(x,p,xgrid,sgrid,agrid_short,prefs,income,Iterating);
        end

        if p.nb == 1
            beta_ub = p.betaH;
        else
            % Don't let highest beta be such that (1-dieprob)*R*beta >= 1
            beta_ub = p.betaH  - max(prefs.betagrid0);
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
    Iterating = 0;
    if p.EpsteinZin == 1
        [~,basemodel] = solve_EGP_EZ(beta_final,p,xgrid,sgrid,agrid_short,prefs,income,Iterating);
    else
        [~,basemodel] = solve_EGP(beta_final,p,xgrid,sgrid,agrid_short,prefs,income,Iterating);
    end
    direct_results.adist = basemodel.adist;

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
    direct_results.mean_s = basemodel.xdist(:)' * basemodel.sav_x(:);
    direct_results.mean_a = basemodel.mean_a;
    direct_results.mean_x = basemodel.xdist(:)' * basemodel.xvals(:);
    direct_results.mean_c = basemodel.xdist(:)' * basemodel.con_x(:);
    
    % One-period income statistics
    direct_results.mean_grossy1 = basemodel.xdist(:)' * basemodel.y_x(:);
    direct_results.mean_loggrossy1 = basemodel.xdist(:)' * log(basemodel.y_x(:));
    direct_results.mean_nety1 = basemodel.xdist(:)' * basemodel.nety_x(:);
    direct_results.mean_lognety1 = basemodel.xdist(:)' * log(basemodel.nety_x(:));
    direct_results.var_loggrossy1 = basemodel.xdist(:)' * (log(basemodel.y_x(:)) - direct_results.mean_loggrossy1).^2;
    direct_results.var_lognety1 = basemodel.xdist(:)' * (log(basemodel.nety_x(:)) - direct_results.mean_lognety1).^2;
    
    direct_results.mean_x_check = direct_results.mean_a + direct_results.mean_nety1;
   
    % Reconstruct yPdist and yFdist from computed stationary distribution 
    % for sanity check
    yPdist_check = reshape(basemodel.adist,[p.nxlong p.nyP p.nyF*p.nb]);
    yPdist_check = sum(sum(yPdist_check,3),1)';
    yFdist_check = reshape(basemodel.adist,[p.nxlong*p.nyP p.nyF p.nb]);
    yFdist_check = sum(sum(yFdist_check,3),1)';
    

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
    if norm(yFdist_check-income.yFdist) > 1e-3
        checks{end+1} = 'Bad_yF_Dist';
    end
    if min(basemodel.adist(:)) < - 1e-3
        checks{end+1} = 'LargeNegativeStateProbability';
    elseif min(basemodel.adist(:)) < -1e-8
        checks{end+1} = 'MedNegativeStateProbability';
    elseif min(basemodel.adist(:)) < -1e-13
        checks{end+1} = 'SmallNegativeStateProbability';
    end

    %% WEALTH DISTRIBUTION
    
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
            direct_results.constrained(i) = basemodel.adist(:)' * (agrid==0);
        else
            direct_results.constrained(i) = wpinterp(p.borrow_lim + p.epsilon(i)*income.meany1*p.freq);
        end
    end
    
    % Wealth percentiles
    [acumdist_unique,uniqueind] = unique(sort_acumdist,'last');
    wpinterp_inverse = griddedInterpolant(acumdist_unique,sort_agrid(uniqueind),'linear');
    direct_results.wpercentiles = wpinterp_inverse(p.percentiles/100);
    
    % Top shares
    % Amount of total assets that reside in each pt on sorted asset space
    totassets = sort_adist .* sort_agrid;
    % Fraction of total assets in each pt on asset space
    cumassets = cumsum(totassets) / direct_results.mean_a;
    
    % create interpolant from wealth percentile to cumulative wealth share
    cumwealthshare = griddedInterpolant(acumdist_unique,cumassets(uniqueind),'linear');
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
        [sim_results,assetmeans] = simulate(p,income,basemodel,xgrid,prefs,agrid);
    else
        assetmeans = [];
    end
    
    %% DIRECTLY COMPUTED 1-PERIOD MPCs
    [avg_mpc1_agrid,mpcs1_a_direct,avg_mpc4_agrid,mpcs4_a_direct,agrid_dist,norisk_mpcs1_a_direct] = ...
                                direct_MPCs_by_computation(p,basemodel,income,prefs,agrid_short,norisk);
    direct_results.avg_mpc1_agrid = avg_mpc1_agrid;
    direct_results.avg_mpc4_agrid = avg_mpc4_agrid;
    direct_results.mpcs1_a_direct = mpcs1_a_direct;
    direct_results.mpcs4_a_direct = mpcs4_a_direct;
    direct_results.agrid_dist = agrid_dist;
    norisk_results.mpcs1_a_direct = norisk_mpcs1_a_direct;
    
    %% MPCs via DRAWING FROM STATIONARY DISTRIBUTION AND SIMULATING
    % Model with income risk
    if p.freq == 4
        [mpcs1,mpcs4,stdev_loggrossy_A,stdev_lognety_A] ...
                            = direct_MPCs_by_simulation(p,prefs,income,basemodel,xgrid,agrid);
        basemodel.mpcs1_sim = mpcs1;
        basemodel.mpcs4_sim = mpcs4;
    else
        basemodel.mpcs1_sim = cell(1,numel(p.mpcfrac));
        basemodel.mpcs4_sim = cell(1,numel(p.mpcfrac));
    end

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
    
    
    
    %% DECOMPOSITION
	decomp = struct([]);
    if p.nb == 1
        m_ra = p.R * (p.beta*p.R)^(-1/p.risk_aver) - 1;
 
        m0 = direct_results.mpcs1_a_direct{5};
        g0 = direct_results.agrid_dist;
        mbc  = norisk_results.mpcs1_a_direct{5};
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
    
    %% GINI
    % Wealth
    direct_results.wealthgini = direct_gini(agrid,basemodel.adist);
    
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