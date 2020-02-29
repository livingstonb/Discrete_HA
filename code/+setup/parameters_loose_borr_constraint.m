function params = parameters_loose_borr_constraint(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    import aux.set_shared_fields

    % location of baseline income process for quarterly case
    quarterly_b_path = 'input/income_quarterly_b_contyT.mat';
    quarterly_c_path = 'input/income_quarterly_c_contyT.mat';

    quarterly_b_params = struct();
    quarterly_b_params.sd_logyT = sqrt(0.6376);
    quarterly_b_params.lambdaT = 0.25;

    quarterly_c_params = struct();
    quarterly_c_params.sd_logyT = sqrt(1.6243);
    quarterly_c_params.lambdaT = 0.0727;

    quarterly_a_params = struct();
    quarterly_a_params.sd_logyT = sqrt(0.2087);
    quarterly_a_params.sd_logyP = sqrt(0.01080);
    quarterly_a_params.rho_logyP = 0.9881;
    quarterly_a_params.lambdaT = 1;

    annual_params = struct();
    annual_params.sd_logyT = sqrt(0.0494);
    annual_params.sd_logyP = sqrt(0.0422);
    annual_params.rho_logyP = 0.9525;
    annual_params.lambdaT = 1;

    num_neg_pts = 100;
    neg_grid_params = struct();
    neg_grid_params.nx = 300 + num_neg_pts;
    neg_grid_params.nx_neg = num_neg_pts;
    neg_grid_params.nx_DST = 300 + num_neg_pts;
    neg_grid_params.nx_neg_DST = num_neg_pts;
    neg_grid_params.borrow_lim = -1e10;

    shared_params = struct();
    shared_params.nyT = 3;

    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    beta_annual = 0.984099818277828;
    beta_quarterly = 0.984368055782111;
    baseline_nbl_adjustment = 0.95;
    
    % Annual
    params(1) = setup.Params(1, 'baseline_A', '');
    params(end) = set_shared_fields(params(end), annual_params);
    params(end) = set_shared_fields(params(end), shared_params);
    params(end).beta0 = beta_annual;
    params(end).other = true;

    % % Annual with borrowing
    % params(end+1) = setup.Params(1, 'baseline_A_with_borrowing', '');
    % params(end).beta0 = beta_annual;
    % params(end).nbl_adjustment = 0.95;
    % params(end).borrow_lim = -1e10;
    % params(end).nx = 500 + num_neg_pts;
    % params(end).nx_neg = num_neg_pts;
    % params(end).nx_DST = 400 + num_neg_pts;
    % params(end).nx_neg_DST = num_neg_pts;

    % Quarterly
    params(end+1) = setup.Params(4, 'baseline_Q', quarterly_b_path);
    params(end) = set_shared_fields(params(end), quarterly_b_params);
    params(end) = set_shared_fields(params(end), shared_params);
    params(end).beta0 = beta_quarterly;
    params(end).other = false;

    % % Quarterly with borrowing
    % num_neg_pts = 20;
    % params(end+1) = setup.Params(4, 'baseline_Q_with_borrowing', quarterly_b_path);
    % params(end) = set_shared_fields(params(end), quarterly_b_params);
    % params(end).beta0 = beta_quarterly;
    % params(end).nbl_adjustment = 0.95;
    % params(end).borrow_lim = -1e10;
    % params(end).nx = 500 + num_neg_pts;
    % params(end).nx_neg = num_neg_pts;
    % params(end).nx_DST = 400 + num_neg_pts;
    % params(end).nx_neg_DST = num_neg_pts;

    %----------------------------------------------------------------------
    % VARYING THE PROXIMITY TO THE NBL
    %----------------------------------------------------------------------
    nbl_factors = [0.8 0.85 0.9 0.95 0.99];

    % Use same discount rate as baseline
    for ii = 1:numel(nbl_factors)
        name = sprintf('baseline_Q_with_borrowing_nbl_test%d', ii);
        params(end+1) = setup.Params(4, name, quarterly_b_path);
        params(end) = set_shared_fields(params(end), quarterly_b_params);
        params(end) = set_shared_fields(params(end), neg_grid_params);
        params(end) = set_shared_fields(params(end), shared_params);
        params(end).beta0 = beta_quarterly;
        params(end).nbl_adjustment = nbl_factors(ii);
        params(end).other = false;
    end

    % Calibrating to match mean wealth
    for ii = 1:numel(nbl_factors)
        name = sprintf('baseline_Q_with_borrowing_nbl_test%d_calibrated', ii);
        params(end+1) = setup.Params(4, name, quarterly_b_path);
        params(end) = set_shared_fields(params(end), quarterly_b_params);
        params(end) = set_shared_fields(params(end), neg_grid_params);
        params(end) = set_shared_fields(params(end), shared_params);
        params(end).beta0 = beta_quarterly;
        params(end).nbl_adjustment = nbl_factors(ii);
        params(end).other = true;
    end

    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS, DO NOT CHANGE
    %----------------------------------------------------------------------
    params.set_index();

    % get list of all names
    all_names = cell2table({params.name}');
    
    % select by number if there is one, otherwise select by names,
    % otherwise use all
    if numel(runopts.number) == 1
        params = params(runopts.number);
    elseif numel(runopts.number) > 1
        error('runopts.number must have 1 or zero elements')
    else
        params = setup.Params.select_by_names(params, runopts.names_to_run);
    end

    params.make_adjustments();
    params.set_run_parameters(runopts);

    %----------------------------------------------------------------------
    % ATTACH CALIBRATOR
    %----------------------------------------------------------------------
    if ~params.other
        params.set("calibrate", false, true);
    end

    if params.calibrate
        heterogeneity = setup.Prefs_R_Heterogeneity(params);
        new_betaH = params.betaH - max(heterogeneity.betagrid0);
        params.set("betaH", new_betaH, true);

        calibrator = aux.mean_wealth_calibrator(params);
        calibrator.set_handle(params);
        params.set("calibrator", calibrator, true);
    end
end