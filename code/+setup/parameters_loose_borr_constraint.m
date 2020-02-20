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

    num_neg_pts = 20;
    neg_grid_params = struct();
    neg_grid_params.nx = 500 + num_neg_pts;
    neg_grid_params.nx_neg = num_neg_pts;
    neg_grid_params.nx_DST = 400 + num_neg_pts;
    neg_grid_params.nx_neg_DST = num_neg_pts;
    neg_grid_params.borrow_lim = -1e10;

    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    beta_annual = 0.984099818277828;
    beta_quarterly = 0.984323067410858;
    baseline_nbl_adjustment = 0.95;
    
    % Annual
    params(1) = setup.Params(1, 'baseline_A', '');
    params(end) = set_shared_fields(params(end), annual_params);
    params(1).beta0 = beta_annual;

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
    params(end).beta0 = beta_quarterly;

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
    nbl_factors = 0.95;
    for ii = 1:numel(nbl_factors)
        name = sprintf('baseline_Q_with_borrowing_nbl_test%d', ii);
        params(end+1) = setup.Params(4, name, quarterly_b_path);
        params(end) = set_shared_fields(params(end), quarterly_b_params);
        params(end) = set_shared_fields(params(end), neg_grid_params);
        params(end).beta0 = beta_quarterly;
        params(end).nbl_adjustment = nbl_factors(ii);
    end

    %----------------------------------------------------------------------
    % ADJUST TO QUARTERLY VALUES, DO NOT CHANGE
    %----------------------------------------------------------------------
    params = setup.Params.adjust_if_quarterly(params);

    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS, DO NOT CHANGE
    %----------------------------------------------------------------------

    params.set_run_parameters(runopts);

    % creates ordered 'index' field
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

    %----------------------------------------------------------------------
    % SET SHARED PARAMETERS
    %----------------------------------------------------------------------
    params.calibrate = false;
    params.nyT = 3;
end