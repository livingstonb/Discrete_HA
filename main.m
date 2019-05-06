function [results,decomp] = main(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    % Greg Kaplan 2017
    
    % This is the main function file for this code repository. Given a
    % structure of parameters, p, this script calls functions primarily to 
    % compute policy functions via the method of endogenous grip points, 
    % and to find the implied stationary distribution over the state space.

    results = struct('direct',[],'norisk',[],'sim',[]);
    results.checks = {};
    
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
    
    grdHJB = struct('s',struct(),'x',struct());
    grdKFE = struct('x',struct());
    
    % savings grids
    sgrid = linspace(0,1,p.nx)';
    sgrid = sgrid.^(1./p.xgrid_par);
    sgrid = p.borrow_lim + (p.xmax-p.borrow_lim).*sgrid;
    
    % Force grid spacing >= gridspace_min near 0
    for ix = 1:p.nx-1
        if sgrid(ix+1) - sgrid(ix) < p.gridspace_min
            sgrid(ix+1) = sgrid(ix) + p.gridspace_min;
        else
            break
        end
    end
    
    grdHJB.s.vec = sgrid;
    grdHJB.s.matrix = repmat(sgrid,[1 p.nyP p.nyF]);

    % xgrids (cash on hand), different min points for each value of (iyP,iyF)
    minyT = kron(min(income.netymat,[],2),ones(p.nx,1));
    xgrid = p.R * grdHJB.s.matrix(:) + minyT;
    grdHJB.x.matrix = reshape(xgrid,[p.nx p.nyP p.nyF]);
    
    % xgrid for model without income risk
    xgrid_norisk  = grdHJB.s.vec + income.meany1;
    for ix = 1:p.nx-1
        if xgrid_norisk(ix+1) - xgrid_norisk(ix) < p.gridspace_min
            xgrid_norisk(ix+1) = xgrid_norisk(ix) + p.gridspace_min;
        else
            break
        end
    end
    grdHJB.x.vec_norisk = xgrid_norisk;
    
    xgrid_norisk = linspace(0,1,p.nx_KFE);
    xgrid_norisk = xgrid_norisk.^(1/p.xgrid_par);
    % Force grid spacing >= gridspace_min near 0
    for ix = 1:p.nx_KFE-1
        if xgrid_norisk(ix+1) - xgrid_norisk(ix) < p.gridspace_min
            xgrid_norisk(ix+1) = xgrid_norisk(ix) + p.gridspace_min;
        else
            break
        end
    end
    xgrid_norisk = p.borrow_lim + (p.xmax-p.borrow_lim).*xgrid_norisk;
    grdKFE.x.vec_norisk = xgrid_norisk + income.meany1;
    
    % create longer xgrid
    minyT = kron(min(income.netymat,[],2),ones(p.nx_KFE,1));
    xgrid= linspace(0,1,p.nx_KFE)';
    xgrid = xgrid.^(1/p.xgrid_par);
    xgrid = p.borrow_lim + (p.xmax - p.borrow_lim)*xgrid;
    % Force grid spacing >= gridspace_min near 0
    for ix = 1:p.nx_KFE-1
        if xgrid(ix+1) - xgrid(ix) < p.gridspace_min
            xgrid(ix+1) = xgrid(ix) + p.gridspace_min;
        else
            break
        end
    end
    xgrid = repmat(xgrid,p.nyP*p.nyF,1);
    xgrid = xgrid + minyT;
    grdKFE.x.matrix = reshape(xgrid,[p.nx_KFE p.nyP p.nyF]);
    
    % Create common agrid to compute mpcs along same agrid for all
    % parameterizations. Dimension nx_KFE x 1
    % Decomp2 only valid between specs with same nx_KFE, xgrid_par,
    % borrow_lim, and xmax
    agrid = linspace(0,1,p.nx_KFE)';
    agrid = agrid.^(1/p.xgrid_par);
    agrid = p.borrow_lim + (p.xmax - p.borrow_lim) * agrid;
    % Force grid spacing >= gridspace_min near 0
    for ia = 1:p.nx_KFE-1
        if agrid(ia+1) - agrid(ia) < p.gridspace_min
            agrid(ia+1) = agrid(ia) + p.gridspace_min;
        else
            break
        end
    end
    grdKFE.a.vec = agrid;
    grdKFE.a.matrix = repmat(agrid,[1,p.nyP,p.nyF,p.nb]);
    
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
        if p.EpsteinZin == 1
            iterate_EGP = @(x) solve_EGP_EZ(x,p,grdHJB,grdKFE,prefs,income);
        else
            iterate_EGP = @(x) solve_EGP(x,p,grdHJB,grdKFE,prefs,income,mpcshock,[]);
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
            results.checks{end+1} = 'NoBetaConv';
            return
        end
    else
        % Beta was set in parameters
        beta_final = p.beta0;    
    end
    
    % Get policy functions and stationary distribution for final beta, in
    % 'basemodel' structure
    if p.EpsteinZin == 1
        [~,basemodel] = solve_EGP_EZ(beta_final,p,grdHJB,grdKFE,prefs,income);
    else
        mpcshock = 0;
        [~,basemodel] = solve_EGP(beta_final,p,grdHJB,grdKFE,prefs,income,mpcshock,[]);
    end
    results.direct.adist = basemodel.adist;

    % Report beta and annualized beta
    results.direct.beta_annualized = beta_final^p.freq;
    results.direct.beta = beta_final;
    
    if basemodel.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        results.checks{end+1} = 'NoEGPConv';
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
    yPdist_check = reshape(basemodel.adist,[p.nx_KFE p.nyP p.nyF*p.nb]);
    yPdist_check = sum(sum(yPdist_check,3),1)';
    yFdist_check = reshape(basemodel.adist,[p.nx_KFE*p.nyP p.nyF p.nb]);
    yFdist_check = sum(sum(yFdist_check,3),1)';
    
    %% --------------------------------------------------------------------
    % RECORD PROBLEMS
    % ---------------------------------------------------------------------
    if  abs((p.targetAY - results.direct.mean_a/(income.meany1*p.freq))/p.targetAY) > 1e-3
        results.checks{end+1} = 'BadAY';
    end
    if abs((results.direct.mean_x-results.direct.mean_x_check)/results.direct.mean_x)> 1e-3
        results.checks{end+1} = 'DistNotStationary';
    end
    if abs((income.meannety1-results.direct.mean_nety1)/income.meannety1) > 1e-3
        results.checks{end+1} = 'BadNetIncomeMean';
    end
    if norm(yPdist_check-income.yPdist) > 1e-3
        results.checks{end+1} = 'Bad_yP_Dist';
    end
    if norm(yFdist_check-income.yFdist) > 1e-3
        results.checks{end+1} = 'Bad_yF_Dist';
    end
    if basemodel.adiff > 1e-8
        results.checks{end+1} = sprintf('NoStatConv, adiff = %5.3e',basemodel.adiff);
    end
    if min(basemodel.adist(:)) < - 1e-3
        results.checks{end+1} = 'LargeNegativeStateProbability';
    elseif min(basemodel.adist(:)) < -1e-8
        results.checks{end+1} = 'MedNegativeStateProbability';
    elseif min(basemodel.adist(:)) < -1e-13
        results.checks{end+1} = 'SmallNegativeStateProbability';
    end

    %% --------------------------------------------------------------------
    % WEALTH DISTRIBUTION
    % ---------------------------------------------------------------------

    % Create values for fraction constrained (HtM) at every pt in asset space,
    % defining constrained as a <= epsilon * mean annual gross labor income 
    % + borrowing limit
    sort_aspace = sortrows([grdKFE.a.matrix(:) basemodel.adist(:)]);
    sort_agrid = sort_aspace(:,1);
    sort_adist = sort_aspace(:,2);

    sort_acumdist = cumsum(sort_adist);

    [aunique,uind] = unique(sort_agrid,'last');
    wpinterp = griddedInterpolant(aunique,sort_acumdist(uind),'linear');
    for i = 1:numel(p.epsilon)        
        % create interpolant to find fraction of constrained households
        if p.epsilon(i) == 0
            % Get exact figure
            results.direct.constrained(i) = basemodel.adist(:)' * (grdKFE.a.matrix(:)==0);

            if p.Bequests == 1
                results.direct.s0 = results.direct.constrained(i);
            else
            	c = results.direct.constrained(i);
                results.direct.s0 = (c - p.dieprob) / (1 - p.dieprob);
            end
        else
            results.direct.constrained(i) = wpinterp(p.epsilon(i)*income.meany1*p.freq);
        end
    end

    % fraction constrained in terms of own quarterly income
%     
%     xrange = 0:0.01:0.5;
%     constrained = zeros(numel(xrange),1);
%     ic = 0;
%     for c = xrange
%         ic = ic + 1;
% 
%         ind = agrid < (ymat_large * (p.freq/4) * c);
%         probabilities = (ind .* basemodel.adist(:)) .* income.yTdist(:)';
%         constrained(ic) = sum(probabilities(:));
%     end
% 
%     constrained_interp = griddedInterpolant(xrange,constrained,'linear');

    
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
    norisk = solve_EGP_deterministic(p,grdHJB,prefs,income,results.direct);
    if norisk.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        results.checks{end+1} = 'NoRiskNoEGPConv';
        return
    end
    
    %% --------------------------------------------------------------------
    % SIMULATIONS
    % ---------------------------------------------------------------------
    if p.Simulate == 1
        [results.sim,assetmeans] = simulate(p,income,basemodel,grdKFE,prefs);
    else
        assetmeans = [];
    end
    
    %% --------------------------------------------------------------------
    % MPCS FOR NO-RISK MODEL
    % ---------------------------------------------------------------------
    
    results.norisk.mpcs1_a_direct = direct_MPCs_by_computation_norisk(p,norisk,income,prefs,grdKFE);

    %% --------------------------------------------------------------------
    % DIRECTLY COMPUTED MPCs, IMPC(s,t)
    % ---------------------------------------------------------------------
    if p.mpcshocks_after_period1 == 1
        maxT = p.freq * 4 + 1;
    else
        maxT = 1;
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
                                    grdKFE.a.vec,prefs,income,nextmpcshock,nextmodel);
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
                = direct_MPCs_by_computation(p,basemodel,mpcmodels,income,prefs,grdKFE,shocksize);
        else
            % epstein-zin preferences, only do (is,it) for is == 1
            
            shocksize = shocks(ishock) * income.meany1 * p.freq;
            [results.direct.mpcs(ishock),~] ...
                = direct_MPCs_by_computation(p,basemodel,mpcmodels,income,prefs,grdKFE,shocksize);
        end
    end
    
    %% --------------------------------------------------------------------
    % MPCs via DRAWING FROM STATIONARY DISTRIBUTION AND SIMULATING
    % ---------------------------------------------------------------------
    % Model with income risk
    MPCs = struct();
    for i = 1:3
        [MPC_trials(i),stdev_loggrossy_A(i),stdev_lognety_A(i),inc_constrained(i)] ...
                            = direct_MPCs_by_simulation(p,prefs,income,basemodel,grdKFE);
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
        meanm0 = results.direct.mpcs(5).avg_s_t{1,1};
        g0 = results.direct.adist; % distribution
        g0_norisk = results.direct.agrid_dist;
        mbc  = results.norisk.mpcs1_a_direct{5}; % norisk distribution
        meanmbc = mbc(:)' * g0_norisk(:);

        % interpolate to get the integral of m0(a) * g0(a) between a = 0 and 0.05
        m0g0 = m0(:) .* g0(:);
        m0g0 = reshape(m0g0,[p.nx_KFE p.nyP*p.nyF*p.nb]);
        m0g0 = sum(m0g0,2);
        cum_m0g0 = cumsum(m0g0);
        mg0interp = griddedInterpolant(grdKFE.a.vec,cum_m0g0,'linear');

        % interpolate to get the integral of mpc_norisk(a) * g0_norisk(a)
        mbcg0 = mbc(:) .* g0_norisk(:);
        mbcg0 = reshape(mbcg0,p.nx_KFE,[]);
        mbcg0 = sum(mbcg0,2);
        cum_mbcg0 = cumsum(mbcg0);
        mbcg0interp = griddedInterpolant(grdKFE.a.vec,cum_mbcg0,'linear');

        % get interpolant for cumulative dist of g0_a
        g0_a = sum(reshape(g0,p.nx_KFE,[]),2);
        g0interp = griddedInterpolant(grdKFE.a.vec,cumsum(g0_a),'linear');

        % get interpolant for cumdist of g0_norisk_a
        g0_norisk_a = sum(reshape(g0_norisk,p.nx_KFE,[]),2);
        g0ninterp = griddedInterpolant(grdKFE.a.vec,cumsum(g0_norisk_a),'linear');


        for ia = 1:numel(p.abars)
            decomp(ia).term1 = m_ra;

            if p.abars(ia) == 0
                zidx = agrid(:) <= p.abars(ia);
                norisk_zidx = grdKFE.a.vec <= p.abars(ia);

                decomp(ia).term2 = (m0(zidx) - m_ra)' * g0(zidx);
                decomp(ia).term3 = (mbc(~norisk_zidx) - m_ra)' * g0_norisk(~norisk_zidx);
                decomp(ia).term4 = m0(~zidx)' * g0(~zidx)- mbc(~norisk_zidx)' * g0_norisk(~norisk_zidx);
            else
                abar = p.abars(ia);
                decomp(ia).term2 = mg0interp(abar) - m_ra * g0interp(abar);
                decomp(ia).term3 = meanmbc - mbcg0interp(abar) - m_ra * (1-g0ninterp(abar));
                decomp(ia).term4 = (meanm0 - mg0interp(abar)) - (meanmbc - mbcg0interp(abar));
                
            end
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
    results.direct.wealthgini = direct_gini(grdKFE.a.matrix,basemodel.adist);
    
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