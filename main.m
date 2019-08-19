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
    results.Finished = false;

    % throw error if more than one type of heterogeneity are added
    if (p.nbeta > 1) + (numel(p.risk_aver)>1) + (numel(p.r)>1) > 1
        error('only one form of heterogeneity allowed')
    else
        % find a better way to do this...
        p.nb = max([numel(p.nbeta),numel(p.risk_aver),numel(p.r)]);
    end
    
    %% --------------------------------------------------------------------
    % HETEROGENEITY IN PREFERENCES/RETURNS
    % ---------------------------------------------------------------------
    heterogeneity = Prefs_R_Heterogeneity(p);
    
    %% --------------------------------------------------------------------
    % INCOME
    % ---------------------------------------------------------------------
    income = gen_income_variables(p,heterogeneity);

    %% --------------------------------------------------------------------
    % ASSET GRIDS
    % ---------------------------------------------------------------------
    
    % grids for method of EGP
    grdEGP = Grid(p,income,'EGP');

    % grids for finding stationary distribution
    grdDST = Grid(p,income,'DST');

    %% --------------------------------------------------------------------
    % MODEL SOLUTION
    % ---------------------------------------------------------------------

    if p.IterateBeta == 1
        
        mpcshock = 0;
        if p.EpsteinZin == 1
            iterate_EGP = @(x) solve_EGP_EZ(x,p,grdEGP,grdDST,heterogeneity,income);
        else
            iterate_EGP = @(x) solve_EGP(x,p,grdEGP,grdDST,heterogeneity,income,mpcshock,[]);
        end

        if numel(heterogeneity.betadist) == 1
            beta_ub = p.betaH;
        else
            % Don't let highest beta be such that (1-dieprob)*R*beta >= 1
            beta_ub = p.betaH  - max(heterogeneity.betagrid0);
        end
        beta_lb = p.betaL;

        % output function that limits number of fzero iterations
        check_evals = @(x,y,z) fzero_checkiter(x,y,z,p.maxiterAY);
        
        options = optimset('TolX',p.tolAY,'OutputFcn',check_evals);
        [beta_final,~,exitflag] = fzero(iterate_EGP,[beta_lb,beta_ub],options);
        if exitflag ~= 1
            return
        end
    else
        % Beta was set in parameters
        beta_final = p.beta0;    
    end
    
    % Get policy functions and stationary distribution for final beta, in
    % 'basemodel' structure
    if p.EpsteinZin == 1
        [~,basemodel] = solve_EGP_EZ(beta_final,p,grdEGP,grdDST,heterogeneity,income);
    else
        mpcshock = 0;
        [~,basemodel] = solve_EGP(beta_final,p,grdEGP,grdDST,heterogeneity,income,mpcshock,[]);
    end
    results.direct.adist = basemodel.adist;

    % Report beta and annualized beta
    results.direct.beta_annualized = beta_final^p.freq;
    results.direct.beta = beta_final;
    
    if basemodel.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
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

    %% --------------------------------------------------------------------
    % WEALTH DISTRIBUTION
    % ---------------------------------------------------------------------

    % Create values for fraction constrained (HtM) at every pt in asset space,
    % defining constrained as a <= epsilon * mean annual gross labor income 
    % + borrowing limit
    sort_aspace = sortrows([grdDST.a.matrix(:) basemodel.adist(:)]);
    sort_agrid = sort_aspace(:,1);
    sort_adist = sort_aspace(:,2);

    sort_acumdist = cumsum(sort_adist);

    [aunique,uind] = unique(sort_agrid,'last');
    wpinterp = griddedInterpolant(aunique,sort_acumdist(uind),'linear');
    for i = 1:numel(p.epsilon)        
        % create interpolant to find fraction of constrained households
        if p.epsilon(i) == 0
            % Get exact figure
            results.direct.constrained(i) = basemodel.adist(:)' * (grdDST.a.matrix(:)==0);

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
    norisk = solve_EGP_deterministic(p,grdEGP,heterogeneity,income,results.direct);
    if norisk.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        return
    end
    
    %% --------------------------------------------------------------------
    % SIMULATIONS
    % ---------------------------------------------------------------------
    if p.Simulate == 1
        results.sim = simulate(p,income,basemodel,grdDST,heterogeneity);
    end
    
    %% --------------------------------------------------------------------
    % MPCS FOR NO-RISK MODEL
    % ---------------------------------------------------------------------
    
    results.norisk.mpcs1_a_direct = ...
        direct_MPCs_by_computation_norisk(p,norisk,income,heterogeneity,grdDST);

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

                    [~,model_lagged{lag}] = solve_EGP(results.direct.beta,p,grdEGP,...
                        grdDST,heterogeneity,income,nextmpcshock,nextmodel);
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
                = direct_MPCs_by_computation(p,basemodel,mpcmodels,...
                    income,grdDST,shocksize);
        else
            % epstein-zin preferences, only do (is,it) for is == 1
            
            shocksize = shocks(ishock) * income.meany1 * p.freq;
            [results.direct.mpcs(ishock),~] ...
                = direct_MPCs_by_computation(p,basemodel,mpcmodels,...
                    income,heterogeneity,grdDST,shocksize);
        end
    end
    
    %% --------------------------------------------------------------------
    % MPCs via DRAWING FROM STATIONARY DISTRIBUTION AND SIMULATING
    % ---------------------------------------------------------------------
    % Model with income risk
    MPCs = struct();
    for i = 1:3
        [MPC_trials(i),stdev_loggrossy_A(i),stdev_lognety_A(i),inc_constrained(i)] ...
            = direct_MPCs_by_simulation(p,heterogeneity,income,basemodel,grdDST);
    end

    results.direct.a_sixth_sim = mean([inc_constrained.a_sixth_Q]);
    results.direct.a_twelfth_sim = mean([inc_constrained.a_twelfth_Q]);
    results.direct.x_sixth_sim = mean([inc_constrained.x_sixth_Q]);
    results.direct.x_twelfth_sim = mean([inc_constrained.x_twelfth_Q]);
    results.direct.a_lt_015_annual = mean([inc_constrained.a_lt_015_annual]);

    MPCs.avg_1_1 = (MPC_trials(1).avg_1_1 + MPC_trials(2).avg_1_1 + MPC_trials(3).avg_1_1)/3;
    MPCs.avg_1_2 = (MPC_trials(1).avg_1_2 + MPC_trials(2).avg_1_2 + MPC_trials(3).avg_1_2)/3;
    MPCs.avg_1_3 = (MPC_trials(1).avg_1_3 + MPC_trials(2).avg_1_3 + MPC_trials(3).avg_1_3)/3;
    MPCs.avg_1_4 = (MPC_trials(1).avg_1_4 + MPC_trials(2).avg_1_4 + MPC_trials(3).avg_1_4)/3;
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
    if p.nb == 1 && p.EpsteinZin == 0 && p.bequest_weight == 0 && p.temptation == 0 && (numel(p.r)==1)
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
        mg0interp = griddedInterpolant(grdDST.a.vec,cum_m0g0,'linear');

        % interpolate to get the integral of mpc_norisk(a) * g0_norisk(a)
        mbcg0 = mbc(:) .* g0_norisk(:);
        mbcg0 = reshape(mbcg0,p.nx_KFE,[]);
        mbcg0 = sum(mbcg0,2);
        cum_mbcg0 = cumsum(mbcg0);
        mbcg0interp = griddedInterpolant(grdDST.a.vec,cum_mbcg0,'linear');

        % get interpolant for cumulative dist of g0_a
        g0_a = sum(reshape(g0,p.nx_KFE,[]),2);
        g0interp = griddedInterpolant(grdDST.a.vec,cumsum(g0_a),'linear');

        % get interpolant for cumdist of g0_norisk_a
        g0_norisk_a = sum(reshape(g0_norisk,p.nx_KFE,[]),2);
        g0ninterp = griddedInterpolant(grdDST.a.vec,cumsum(g0_norisk_a),'linear');


        for ia = 1:numel(p.abars)
            decomp(ia).term1 = m_ra;

            if p.abars(ia) == 0
                zidx = grdDST.a.matrix(:) <= p.abars(ia);
                norisk_zidx = grdDST.a.vec <= p.abars(ia);

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
    results.direct.wealthgini = direct_gini(grdDST.a.matrix,basemodel.adist);
    
    % Gross income
    results.direct.grossincgini = direct_gini(income.ysort,income.ysortdist);
    
    % Net income
    results.direct.netincgini = direct_gini(income.netymat,income.ymatdist);  

    results.Finished = true; 

    function gini = direct_gini(level,distr)
        % Sort distribution and levels by levels
        sorted = sortrows([level(:),distr(:)]);
        level_sort = sorted(:,1);
        dist_sort  = sorted(:,2);
        S = [0;cumsum(dist_sort .* level_sort)];
        gini = 1 - dist_sort' * (S(1:end-1)+S(2:end)) / S(end);
    end
    
end