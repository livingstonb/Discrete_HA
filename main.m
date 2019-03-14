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
    p.savtaxthresh  = p.savtaxthresh * income.meany1 * p.freq;

    if (p.nb > 1) && (numel(p.risk_aver) > 1)
    	error('cannot have both beta and risk aversion heterogeneity')
    end

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
    % IES Heterogeneity
    % ---------------------------------------------------------------------
    % same way as introducing beta heterogeneity

    if numel(p.risk_aver) > 1 || ((numel(p.invies) > 1) && (p.EpsteinZin == 1))
    	% use nb as grid size for risk aversion
        if numel(p.risk_aver) > 1
            p.nb = numel(p.risk_aver);
        else
            p.nb = numel(p.invies);
        end

    	prefs.IESdist = ones(p.nb,1) / p.nb;
    	IESswitch_ij = p.IESswitch / (p.nb-1);

    	diagonal = (1-p.IESswitch) * ones(p.nb,1);
    	off_diag = IESswitch_ij * ones(p.nb);
    	off_diag = off_diag - diag(diag(off_diag));
    	prefs.IEStrans = off_diag + diag(diagonal);
        
        prefs.IEScumdist = cumsum(prefs.IESdist);
        prefs.IEScumtrans = cumsum(prefs.IEStrans,2);
    else
        prefs.IESdist = 1;
        prefs.IEStrans = 0;
        prefs.IEScumdist = 1;
        prefs.IEScumtrans = 0;
    end

    %% --------------------------------------------------------------------
    % ASSET GRIDS
    % ---------------------------------------------------------------------
    
    % savings grids
    sgrid.orig = linspace(0,1,p.nx)';
    sgrid.orig = sgrid.orig.^(1./p.xgrid_par);
    sgrid.orig = p.borrow_lim + (p.xmax-p.borrow_lim).*sgrid.orig;
    sgrid.short = sgrid.orig;
    
    % Force grid spacing >= gridspace_min near 0
    for ix = 1:p.nx-1
        if sgrid.short(ix+1) - sgrid.short(ix) < p.gridspace_min
            sgrid.short(ix+1) = sgrid.short(ix) + p.gridspace_min;
        else
            break
        end
    end
    
    sgrid.full = repmat(sgrid.short,[1 p.nyP p.nyF]);

    % xgrids (cash on hand), different min points for each value of (iyP,iyF)
    minyT               = kron(min(income.netymat,[],2),ones(p.nx,1));
    xgrid.orig          = p.R * sgrid.full(:) + minyT;
    xgrid.full          = reshape(xgrid.orig,[p.nx p.nyP p.nyF]);
    
    % xgrid for model without income risk
    xgrid.norisk_short  = sgrid.short + income.meany1;
    for ix = 1:p.nx-1
        if xgrid.norisk_short(ix+1) - xgrid.norisk_short(ix) < p.gridspace_min
            xgrid.norisk_short(ix+1) = xgrid.norisk_short(ix) + p.gridspace_min;
        else
            break
        end
    end
    xgrid.norisk_longgrid = linspace(0,1,p.nxlong);
    xgrid.norisk_longgrid = xgrid.norisk_longgrid.^(1/p.xgrid_par);
    % Force grid spacing >= gridspace_min near 0
    for ix = 1:p.nxlong-1
        if xgrid.norisk_longgrid(ix+1) - xgrid.norisk_longgrid(ix) < p.gridspace_min
            xgrid.norisk_longgrid(ix+1) = xgrid.norisk_longgrid(ix) + p.gridspace_min;
        else
            break
        end
    end
    xgrid.norisk_longgrid = p.borrow_lim + (p.xmax-p.borrow_lim).*xgrid.norisk_longgrid;
    xgrid.norisk_longgrid = xgrid.norisk_longgrid + income.meany1;
    
    % create longer xgrid
    minyT = kron(min(income.netymat,[],2),ones(p.nxlong,1));
    xgrid.longgrid = linspace(0,1,p.nxlong)';
    xgrid.longgrid = xgrid.longgrid.^(1/p.xgrid_par);
    xgrid.longgrid = p.borrow_lim + (p.xmax - p.borrow_lim)*xgrid.longgrid;
    % Force grid spacing >= gridspace_min near 0
    for ix = 1:p.nxlong-1
        if xgrid.longgrid(ix+1) - xgrid.longgrid(ix) < p.gridspace_min
            xgrid.longgrid(ix+1) = xgrid.longgrid(ix) + p.gridspace_min;
        else
            break
        end
    end
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
    % Force grid spacing >= gridspace_min near 0
    for ia = 1:p.nxlong-1
        if agrid(ia+1) - agrid(ia) < p.gridspace_min
            agrid(ia+1) = agrid(ia) + p.gridspace_min;
        else
            break
        end
    end
    agrid_short = agrid;
    agrid = repmat(agrid,p.nyP*p.nyF*p.nb,1);
    
    %% --------------------------------------------------------------------
    % UTILITY FUNCTION, BEQUEST FUNCTION
    % ---------------------------------------------------------------------
    % utility function
    if numel(p.risk_aver) == 1
	    if p.risk_aver==1
	        prefs.u = @(c)log(c);
	    else    
	        prefs.u = @(c)(c.^(1-p.risk_aver)-1)./(1-p.risk_aver);
        end    
	else
		% risk_aver heterogeneity, preferences defined in utility.m
		prefs.u = @(risk_aver,c) utility(risk_aver,c);
	end
    
    % bequest utility
    if p.bequest_curv == 1
        prefs.beq = @(a) p.bequest_weight.* log(a+ p.bequest_luxury);
    else
        prefs.beq = @(a) p.bequest_weight.*((a+p.bequest_luxury).^(1-p.bequest_curv)-1)./(1-p.bequest_curv);
    end

    % 1st derivative of utility wrt c
    if numel(p.risk_aver) == 1
    	prefs.u1 = @(c) c.^(-p.risk_aver);
    else
    	prefs.u1 = @(risk_aver,c) c.^(-risk_aver);
    end

    % inverse of 1st derivative
    if numel(p.risk_aver) == 1
    	prefs.u1inv = @(u) u.^(-1./p.risk_aver);
    else
    	prefs.u1inv = @(risk_aver,u) u.^(-1./risk_aver);
    end

    % first derivative of bequest utility
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

    % Create values for fraction constrained (HtM) at every pt in asset space,
    % defining constrained as a <= epsilon * mean annual gross labor income 
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

            if p.Bequests == 1
                results.direct.s0 = results.direct.constrained(i);
            else
                results.direct.s0 = results.direct.constrained(i) / (1-p.deathrate);
            end
        else
            results.direct.constrained(i) = wpinterp(p.epsilon(i)*income.meany1*p.freq);
        end
    end

    % HtM with different income frequency
    ymat_large = kron(income.netymat,ones(p.nxlong,1));
    ymat_large = repmat(ymat_large,p.nb,1);

    % fraction constrained in terms of own quarterly income
    xrange = 0:0.01:0.5;
    constrained = zeros(numel(xrange),1);
    ic = 0;
    for c = xrange
        ic = ic + 1;

        ind = agrid < (ymat_large * (p.freq/4) * c);
        probabilities = (ind .* basemodel.adist(:)) .* income.yTdist(:)';
        constrained(ic) = sum(probabilities(:));
    end

    constrained_interp = griddedInterpolant(xrange,constrained,'linear');

    % the following are innacurate (at < y_{t+1}/6 instead of at < yt /6)
    % results.direct.HtM_one_sixth_Q = constrained_interp(1/6);
    % results.direct.HtM_one_twelfth_Q = constrained_interp(1/12);


 %    % 1/6 quarterly income (1/24 annual income)
	% one_sixth_quarterly = agrid < (ymat_large * p.freq / 24);
 %    probabilities = (one_sixth_quarterly .* basemodel.adist(:)) .* income.yTdist(:)';
 %    results.direct.HtM_one_sixth_Q = sum(probabilities(:));

	% % 1/12 quarterly income (1/48 annual income)
 %    one_twelfth_quarterly = agrid < (ymat_large * p.freq / 48);
 %    probabilities = (one_twelfth_quarterly .* basemodel.adist(:)) .* income.yTdist(:)';
 %    results.direct.HtM_one_twelfth_Q = sum(probabilities(:));

    
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
    
    % save adist from model
    results.direct.adist = basemodel.adist;
    
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
        [results.sim,assetmeans] = simulate(p,income,basemodel,xgrid,agrid_short,prefs);
    else
        assetmeans = [];
    end
    
    %% --------------------------------------------------------------------
    % MPCS FOR NO-RISK MODEL
    % ---------------------------------------------------------------------
    
    results.norisk.mpcs1_a_direct = direct_MPCs_by_computation_norisk(p,norisk,income,prefs,agrid_short);

    %% --------------------------------------------------------------------
    % DIRECTLY COMPUTED MPCs, IMPC(s,t)
    % ---------------------------------------------------------------------
    if p.mpcshocks_after_period1 == 1
        maxT = p.freq * 4;
    else
        maxT = 1
    end
    mpcmodels = cell(maxT,maxT);
    
    shocks = [-1e-5 -0.01 -0.1 1e-5 0.01 0.1];
    
    % policy functions are the same as baseline when shock is received in
    % the current period
    
    for ishock = 1:6
        
        for is = 1:maxT
            mpcmodels{is,is} = basemodel;
        end

        if p.EpsteinZin == 0
            % mpcmodels(s,t) stores the policy functions associated with the case
            % where the household is currently in period t, but recieved news about
            % the period-s shock in period 1
            model_lagged = cell(maxT-1);

            % get consumption functions conditional on future shock
            % 'lag' is number of periods before shock
            if shocks(ishock) > 0 && (maxT > 1)
                for lag = 1:maxT-1
                    Iterating = 0;

                    if lag == 1
                        % shock is next period
                        nextmpcshock = shocks(ishock) * income.meany1 * p.freq;
                        nextmodel = basemodel;
                    else
                        % no shock next period
                        nextmpcshock = 0;
                        nextmodel = model_lagged{lag-1};
                    end

                    [~,model_lagged{lag}] = solve_EGP(results.direct.beta,p,xgrid,sgrid,...                   
                                    agrid_short,prefs,income,Iterating,nextmpcshock,nextmodel);
                end

                % populate mpcmodels with remaining (s,t) combinations for t < s
                for is = 2:maxT
                for it = is-1:-1:1
                    mpcmodels{is,it} = model_lagged{is-it};
                end
                end
            end

            shocksize = shocks(ishock) * income.meany1 * p.freq;
            [results.direct.mpcs(ishock),results.direct.agrid_dist] ...
                = direct_MPCs_by_computation(p,basemodel,mpcmodels,income,prefs,xgrid,agrid_short,shocksize);
        else
            % epstein-zin preferences, only do (is,it) for is == 1
            
            shocksize = shocks(ishock) * income.meany1 * p.freq;
            [results.direct.mpcs(ishock),~] ...
                = direct_MPCs_by_computation(p,basemodel,mpcmodels,income,prefs,xgrid,agrid_short,shocksize);
        end
    end
    
    %% --------------------------------------------------------------------
    % MPCs via DRAWING FROM STATIONARY DISTRIBUTION AND SIMULATING
    % ---------------------------------------------------------------------
    % Model with income risk
    MPCs = struct();
    for i = 1:3
        [MPC_trials(i),stdev_loggrossy_A(i),stdev_lognety_A(i),inc_constrained(i)] ...
                            = direct_MPCs_by_simulation(p,prefs,income,basemodel,xgrid,agrid);
    end
    
    results.direct.a_sixth_sim = mean([inc_constrained.a_sixth_Q]);
    results.direct.a_twelfth_sim = mean([inc_constrained.a_twelfth_Q]);
    results.direct.x_sixth_sim = mean([inc_constrained.x_sixth_Q]);
    results.direct.x_twelfth_sim = mean([inc_constrained.x_twelfth_Q]);
    
    MPCs.avg_1_1 = (MPC_trials(1).avg_1_1 + MPC_trials(2).avg_1_1 + MPC_trials(3).avg_1_1)/3;
    MPCs.avg_1_2 = (MPC_trials(1).avg_1_2 + MPC_trials(2).avg_1_2 + MPC_trials(3).avg_1_2)/3;
    MPCs.avg_1_3 = (MPC_trials(1).avg_1_3 + MPC_trials(2).avg_1_3 + MPC_trials(3).avg_1_3)/3;
    MPCs.avg_1_4 = (MPC_trials(1).avg_1_4 + MPC_trials(2).avg_1_4 + MPC_trials(3).avg_1_4)/3;
    MPCs.avg_1_1to4 = (MPC_trials(1).avg_1_1to4 + MPC_trials(2).avg_1_1to4 + MPC_trials(3).avg_1_1to4)/3;
    MPCs.avg_1_5to8 = (MPC_trials(1).avg_1_5to8 + MPC_trials(2).avg_1_5to8 + MPC_trials(3).avg_1_5to8)/3;
    MPCs.avg_1_9to12 = (MPC_trials(1).avg_1_9to12 + MPC_trials(2).avg_1_9to12 + MPC_trials(3).avg_1_9to12)/3;
    MPCs.avg_1_13to16 = (MPC_trials(1).avg_1_13to16 + MPC_trials(2).avg_1_13to16 + MPC_trials(3).avg_1_13to16)/3;
    results.direct.mpcs_sim = MPCs;
    
    stdev_loggrossy_A = mean(stdev_loggrossy_A);
    stdev_lognety_A = mean(stdev_lognety_A);

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
    	% RA MPC
        m_ra = p.R * (results.direct.beta*p.R)^(-1/p.risk_aver) - 1;
 
        % MPC shock of 0.01 * annual income
        m0 = results.direct.mpcs(5).mpcs_1_t{1,1}; % mpcs
        g0 = results.direct.adist; % distribution
        g0_norisk = results.direct.agrid_dist;
        mbc  = results.norisk.mpcs1_a_direct{5}; % norisk distribution
        for ia = 1:numel(p.abars)
            zidx = agrid(:) <= p.abars(ia);
            norisk_zidx = agrid_short <= p.abars(ia);
            
            decomp(ia).term1 = m_ra;
            decomp(ia).term2 = (m0(zidx) - m_ra)' * g0(zidx);
            decomp(ia).term3 = (mbc(~norisk_zidx) - m_ra)' * g0_norisk(~norisk_zidx);
            decomp(ia).term4 = m0(~zidx)' * g0(~zidx)- mbc(~norisk_zidx)' * g0_norisk(~norisk_zidx);
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