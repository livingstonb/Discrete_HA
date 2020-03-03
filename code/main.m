function results = main(p)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    
    % This is the main function file for this code repository. Given a
    % structure of parameters, p, this script calls functions primarily to 
    % compute policy functions via the method of endogenous grip points, 
    % and to find the implied stationary distribution over the state space.

    results = struct('policy',[],'direct',[],'norisk',[],'sim',[],'decomp_meanmpc',[]);
    results.Finished = false;

    % throw error if more than one type of heterogeneity are added
    if (p.nbeta > 1) + (numel(p.risk_aver)>1) + (numel(p.r)>1)...
        + (numel(p.temptation)>1) + (numel(p.invies)>1) > 1

        error('only one form of heterogeneity allowed')
    else
        % find a better way to do this...
        p.nb = max([p.nbeta, numel(p.risk_aver), numel(p.r), numel(p.invies), numel(p.temptation)]);
    end

    %% --------------------------------------------------------------------
    % HETEROGENEITY IN PREFERENCES/RETURNS
    % ---------------------------------------------------------------------
    heterogeneity = setup.Prefs_R_Heterogeneity(p);

    %% --------------------------------------------------------------------
    % INCOME
    % ---------------------------------------------------------------------
    income = setup.Income(p, heterogeneity);

    %% --------------------------------------------------------------------
    % ASSET GRIDS
    % ---------------------------------------------------------------------
    NBL = -min(income.netymat(:)) / max(p.r);
    loose_constraint = p.nbl_adjustment * NBL;
    if p.borrow_lim <= -1e10
        p.set("borrow_lim", loose_constraint, false);
    end
    
    % grids for method of EGP
    grdEGP = setup.Grid(p, income, 'EGP');

    % grids for finding stationary distribution
    grdDST = setup.Grid(p, income, 'DST');

    %% --------------------------------------------------------------------
    % MODEL SOLUTION
    % ---------------------------------------------------------------------
    % Get policy functions and stationary distribution for final beta, in
    % 'basemodel' structure
    if p.EpsteinZin
        egp_ez_solver = solver.EGP_EZ_Solver(p, grdEGP, heterogeneity, income);
        egp_ez_solver.solve(income);
        basemodel = egp_ez_solver.return_model();
    else
        mpcshock = 0;
        basemodel = solver.solve_EGP(...
            p, grdEGP, heterogeneity, income, mpcshock, []);
    end
    basemodel = solver.find_stationary_adist(...
        p, basemodel, income, grdDST, heterogeneity);
    results.direct.adist = basemodel.adist;

    % Report beta and annualized beta
    results.direct.beta_annualized = p.beta0 ^ p.freq;
    results.direct.beta = p.beta0;
    
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
    % pmf(a, yP, yF, ib)
    results.direct.adist = basemodel.adist;
    
    % pmf(a)
    results.direct.agrid_dist = sum(sum(sum(basemodel.adist,4),3),2);
    
    % Support of the distribution
    support = results.direct.agrid_dist > 1e-9;
    results.direct.agrid = grdDST.a.vec;
    results.direct.agrid_support = grdDST.a.vec(support);
    
    % cdf(a) over the support of pmf(a)
    tmp = cumsum(results.direct.agrid_dist);
    tmp = tmp(support);
    tmp = tmp / tmp(end);
    results.direct.agrid_cumdist_support = tmp;
    
    % pmf(a) over the support of pmf(a)
    tmp = results.direct.agrid_dist(support);
    results.direct.agrid_dist_support = tmp / sum(tmp);

    % Create values for fraction constrained (HtM) at every pt in asset space,
    % defining constrained as a <= epsilon * mean annual gross labor income
    wpinterp = griddedInterpolant(...
        results.direct.agrid_support,...
        results.direct.agrid_cumdist_support,...
        'linear');

    results.direct.find_wealth_pctile = @(a) 100 * wpinterp(a);
    for i = 1:numel(p.epsilon)        
        % create interpolant to find fraction of constrained households
        if p.epsilon(i) == 0
            % Get exact figure
            results.direct.constrained(i) = results.direct.agrid_dist(:)' * (grdDST.a.vec==0);

            if p.Bequests
                results.direct.s0 = results.direct.constrained(i);
            else
            	c = results.direct.constrained(i);
                results.direct.s0 = (c - p.dieprob) / (1 - p.dieprob);
            end
        else
            results.direct.constrained(i) = wpinterp(p.epsilon(i));
        end
    end

    results.direct.wealth_lt_1000 = wpinterp(0.0081*2);
    results.direct.wealth_lt_5000 = wpinterp(0.081);
    results.direct.wealth_lt_10000 = wpinterp(0.081*2);
    
    % Wealth percentiles
    wpinterp_inverse = griddedInterpolant(...
        results.direct.agrid_cumdist_support,...
        results.direct.agrid_support,...
        'linear');
    results.direct.wpercentiles = wpinterp_inverse(p.percentiles/100);
    results.direct.median_a = wpinterp_inverse(0.5);
    
    % Top shares
    % Amount of total assets that reside in each pt on sorted asset space
    totassets = results.direct.agrid_dist_support .* results.direct.agrid_support;
    % Fraction of total assets in each pt on asset space
    cumassets = cumsum(totassets) / results.direct.mean_a;
    
    % Create interpolant from wealth percentile to cumulative wealth share
    cumwealthshare_interp = griddedInterpolant(...
        results.direct.agrid_cumdist_support,...
        cumassets,...
        'linear');
    results.direct.top10share = 1 - cumwealthshare_interp(0.9);
    results.direct.top1share = 1 - cumwealthshare_interp(0.99);
    
    %% --------------------------------------------------------------------
    % MPCs FOR MODEL WITHOUT INCOME RISK
    % ---------------------------------------------------------------------
    norisk = struct('complete', false);
    if p.DeterministicMPCs
        % Solve deterministic model
        norisk = solver.solve_EGP_deterministic(...
            p, grdEGP, heterogeneity);

        if norisk.complete
            % Compute MPCs for deterministic model
            results.norisk.mpcs1_a_direct = ...
                statistics.direct_MPCs_by_computation_norisk(...
                    p, norisk, income, heterogeneity, grdDST);
        else
            p.set('DeterministicMPCs', false, true);
        end
    end

    if ~norisk.complete
        % Fill deterministic MPCs with NaNs
        results.norisk_mpcs1_a_direct = cell(1,6);
        for im = 1:6
            results.norisk_mpcs1_a_direct{im} = NaN;
        end
    end


    %% --------------------------------------------------------------------
    % SIMULATIONS
    % ---------------------------------------------------------------------
    if p.Simulate
        results.sim = solver.simulate(...
            p, income, basemodel, grdDST, heterogeneity);
    end

    %% --------------------------------------------------------------------
    % MPCs over cash-on-hand
    % ---------------------------------------------------------------------
    con_base = basemodel.con;
    for ishock = 1:numel(p.shocks)
        shock_size = p.shocks(ishock);

        con_shock = zeros(p.nx, p.nyP, p.nyF, p.nb);
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            cash_shock = grdEGP.x.matrix(:,iyP,iyF,ib) + shock_size;
            con_shock(:,iyP,iyF,ib) = basemodel.coninterp{iyP,iyF,ib}(cash_shock);
        end
        end
        end

        mpcs = (con_shock - con_base) / shock_size;
        results.direct.mpcs_cash{ishock} = mpcs;
    end

    %% --------------------------------------------------------------------
    % DIRECTLY COMPUTED MPCs, IMPC(s,t)
    % ---------------------------------------------------------------------
    if (p.MPCs_news == 1) || (p.MPCs_loan_and_loss == 1)
        disp('Solving for policy functions of anticipated future shocks')
        if p.freq == 4
            maxT = 10;
        else
            maxT = 5;
        end
    else
        maxT = 1;
    end
    mpcmodels = cell(6, maxT, maxT);
    
    shocks = p.shocks;
    
    % policy functions are the same as baseline when shock is received in
    % the current period
    
    for ishock = 1:6
        for is = 1:maxT
            mpcmodels{ishock,is,is} = basemodel;
        end

        if ~p.EpsteinZin
            % mpcmodels{ishock,s,t} stores the policy functions associated with the case
            % where the household is currently in period t, but recieved news about
            % the period-s shock in period 1. Shock was of size shocks(ishock)
            model_lagged = cell(6,maxT-1);

            % get consumption functions conditional on future shock
            % 'lag' is number of periods before shock
            if maxT > 1
                for lag = 1:maxT-1
                    if lag == 1
                        % shock is next period
                        nextmpcshock = shocks(ishock);
                        nextmodel = basemodel;
                    else
                        % no shock next period
                        nextmpcshock = 0;
                        nextmodel = model_lagged{ishock,lag-1};
                    end

                    model_lagged{ishock,lag} = solver.solve_EGP(...
                        p, grdEGP, heterogeneity, income, nextmpcshock, nextmodel);
                end

                % populate mpcmodels with remaining (s,t) combinations for t < s
                for is = 2:maxT
                for it = is-1:-1:1
                    mpcmodels{ishock,is,it} = model_lagged{ishock,is-it};
                end
                end
            end
        end
    end
    
    mpc_finder = statistics.MPCFinder(p, income, grdDST, heterogeneity,...
        basemodel, mpcmodels);
    if p.MPCs
        disp('Computing MPCs')
        mpc_finder.solve(p, grdDST);
    end

    results.direct.mpcs = mpc_finder.mpcs;
    results.direct.loan = mpc_finder.loan;
    results.direct.loss_in_2_years = mpc_finder.loss_in_2_years;
    clear mpc_finder
    
    %% --------------------------------------------------------------------
    % MPCs via DRAWING FROM STATIONARY DISTRIBUTION AND SIMULATING
    % ---------------------------------------------------------------------
    mpc_simulator = statistics.MPCSimulator(...
        p, results.direct.find_wealth_pctile, heterogeneity);
    mpc_simulator.simulate(...
        p, income, grdDST, heterogeneity, basemodel);

    results.direct = mpc_simulator.append_results(results.direct);

    % find annual mean and standard deviations of income
    if p.freq == 4
        % direct computations
        results.direct.mean_grossy_A = results.direct.mean_grossy1 * 4;
        % from simulations
        results.direct.stdev_loggrossy_A = mpc_simulator.stdev_loggrossy_A;
        results.direct.stdev_lognety_A = mpc_simulator.stdev_lognety_A;     
    else
        % direct computations
        results.direct.mean_grossy_A = results.direct.mean_grossy1;
        results.direct.stdev_loggrossy_A = sqrt(results.direct.var_loggrossy1);
        results.direct.stdev_lognety_A = sqrt(results.direct.var_lognety1);
    end
    
    clear mpc_simulator

    %% --------------------------------------------------------------------
    % DECOMPOSITION 1 (DECOMP OF E[mpc])
    % ---------------------------------------------------------------------
    decomp = statistics.Decomp(p, results.direct);

    doDecomposition = (p.nb==1) && (~p.EpsteinZin) && (p.MPCs)...
        && (p.bequest_weight==0) && isequal(p.temptation,0) && (numel(p.r)==1)...
        && (p.DeterministicMPCs);

    if doDecomposition
        mpcs_baseline = results.direct.mpcs(5).mpcs_1_t{1};
        mpcs_baseline = reshape(mpcs_baseline, p.nx_DST, []);
        mpcs_norisk = results.norisk.mpcs1_a_direct{5};
        mpcs_norisk = reshape(mpcs_norisk, p.nx_DST, []);
        decomp.perform_decompositions(mpcs_baseline, mpcs_norisk);
    end

    results.decomp_RA = decomp.results_RA;
    results.decomp_norisk = decomp.results_norisk;
    clear decomp
    
    %% --------------------------------------------------------------------
    % GINI
    % ---------------------------------------------------------------------
    % Wealth
    results.direct.wealthgini = aux.direct_gini(grdDST.a.matrix, basemodel.adist);
    
    % Gross income
    results.direct.grossincgini = aux.direct_gini(income.ysort, income.ysortdist);
    
    % Net income
    results.direct.netincgini = aux.direct_gini(income.netymat, income.ymatdist);  

    results.Finished = true;
    
    %% --------------------------------------------------------------------
    % FIGURES
    % ---------------------------------------------------------------------
    if p.MakePlots
        % plot(grdDST.a.vec, cumsum(results.direct.agrid_dist))
        % xlim([0 0.2])
        
        % Wealth at low yP
        nbins = 100;
        amin = grdDST.a.vec(1);
        amax = {1};
        amax_visible = 0.5;

        iyP = 1:11;
        nyP = numel(iyP);

        pmf_a = results.direct.adist(:,iyP,:,:);
        pmf_a = pmf_a(:) / sum(pmf_a(:));
        pmf_a = reshape(pmf_a, [], nyP);
        pmf_a = sum(pmf_a, 2);

        wealth_plotter = statistics.WealthPlotter(p, grdDST.a.vec, pmf_a);
        [ax, wealth_hist] = wealth_plotter.create_histogram(nbins, amax{:});
        title("Wealth distribution, truncated above")
        ax.XLim = [amin, amax_visible];
        ax.YLim = [0, max(wealth_hist.Values(1:end-1))];

        figpath = fullfile(p.outdir, 'wealth_distribution.jpg');
        saveas(gcf, figpath)
        
        %% MPCs Function
        fontsize = 12;
        mpcs = results.direct.mpcs(5).mpcs_1_t{1};
        mpc_plotter = statistics.MPCPlotter(p, grdDST.a.matrix, mpcs);
        mpc_plotter.fontsize = fontsize;
        mpc_plotter.show_grid = 'on';

        yP_indices = [3, 8];
        zoomed_window = true;
        shock_size = 0.01;
        [ax_main, ax_window] = mpc_plotter.create_mpcs_plot_yPs(...
                    yP_indices, zoomed_window, shock_size);
        ylim_main = ax_main.YLim;

        imedian = find(p.percentiles == 50);
        median_wealth = results.direct.wpercentiles(imedian);
        ax_main = mpc_plotter.add_median_wealth(ax_main, median_wealth);

        ax_main.XLim = [0, 5];
        ax_main.YLim = ylim_main;

        window_max_x = 0.3;
        ax_window.YLim = ax_main.YLim;
        ax_window.XLim = [0, window_max_x];
        xticks(ax_window, [0:0.1:window_max_x])
        yticks(ax_window, [0:0.1:0.3])
        set(ax_window, 'FontSize', fontsize-2)
        ax_window.YTick = ax_main.YTick(1:2:end);

        figpath = fullfile(p.outdir, 'mpc_function_yPs.jpg');
        saveas(gcf, figpath)

        %% MPCs Function For Diff Shock Sizes
        fontsize = 12;
        mpcs = {    results.direct.mpcs(2).mpcs_1_t{1}
                    results.direct.mpcs(3).mpcs_1_t{1}
                    results.direct.mpcs(5).mpcs_1_t{1}
                    results.direct.mpcs(6).mpcs_1_t{1}
               };
           
        for ii = 1:numel(mpcs)
            mpcs{ii} = reshape(mpcs{ii}, [p.nx_DST p.nyP p.nyF p.nb]);
        end

        mpc_plotter = statistics.MPCPlotter(p, grdDST.a.matrix, mpcs);
        mpc_plotter.fontsize = fontsize;
        mpc_plotter.show_grid = 'on';

        iyP = median(1:p.nyP);
        ishocks = [2 3 5 6];
        zoomed_window = true;
        shock_size = 0.01;
        [ax_main, ax_window] = mpc_plotter.create_mpc_plot_shocks(...
                    iyP, zoomed_window, ishocks);
        ylim_main = ax_main.YLim;

        imedian = find(p.percentiles == 50);
        median_wealth = results.direct.wpercentiles(imedian);
        ax_main = mpc_plotter.add_median_wealth(ax_main, median_wealth);

        ax_main.XLim = [0, 5];
        ax_main.YLim = ylim_main;

        window_max_x = 0.3;
        ax_window.YLim = ax_main.YLim;
        ax_window.XLim = [0, window_max_x];
        xticks(ax_window, [0:0.1:window_max_x])
        yticks(ax_window, [0:0.1:0.3])
        set(ax_window, 'FontSize', fontsize-2)
        ax_window.YTick = ax_main.YTick(1:2:end);

        figpath = fullfile(p.outdir, 'mpc_function_shocks.jpg');
        saveas(gcf, figpath)
    end
    
    % convert Params object to structure for saving
    Sparams = aux.to_structure(p);
    save(p.savematpath, 'Sparams', 'results')
end