function params = parameters_loose_borr_constraint(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    % location of baseline income process for quarterly case
    QIncome = 'input/income_quarterly_b_truncated.mat';
    
    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    beta_annual = 0.984108034755346;
    beta_quarterly = 0.984363510593659;
    baseline_nbl_adjustment = 0.95;
    num_neg_pts = 20;
    
    % Annual
    params(1) = setup.Params(1, 'baseline_A', '');
    params(1).beta0 = beta_annual;

    % Annual with borrowing
    params(end+1) = setup.Params(1, 'baseline_A_with_borrowing', '');
    params(end).beta0 = beta_annual;
    params(end).nbl_adjustment = 0.95;
    params(end).borrow_lim = -1e10;
    params(end).nx = 500 + num_neg_pts;
    params(end).nx_neg = num_neg_pts;
    params(end).nx_DST = 400 + num_neg_pts;
    params(end).nx_neg_DST = num_neg_pts;

    % Quarterly
    params(end+1) = setup.Params(4, 'baseline_Q', QIncome);
    params(end).beta0 = beta_quarterly;
    params(end).Nsim = 1e5;

    % Quarterly with borrowing
    num_neg_pts = 20;
    params(end+1) = setup.Params(4, 'baseline_Q_with_borrowing', QIncome);
    params(end).beta0 = beta_quarterly;
    params(end).nbl_adjustment = 0.95;
    params(end).borrow_lim = -1e10;
    params(end).nx = 500 + num_neg_pts;
    params(end).nx_neg = num_neg_pts;
    params(end).nx_DST = 400 + num_neg_pts;
    params(end).nx_neg_DST = num_neg_pts;

    %----------------------------------------------------------------------
    % VARYING THE PROXIMITY TO THE NBL
    %----------------------------------------------------------------------
    nbl_factors = 0.8:0.01:0.98
    for ii = 1:numel(nbl_factors)
        name = sprintf('baseline_Q_with_borrowing_nbl_test%d', ii);
        params(end+1) = setup.Params(4, name, QIncome);
        params(end).beta0 = beta_quarterly;
        params(end).nbl_adjustment = nbl_factors(ii);
        params(end).borrow_lim = -1e10;
        params(end).nx = 500 + num_neg_pts;
        params(end).nx_neg = num_neg_pts;
        params(end).nx_DST = 400 + num_neg_pts;
        params(end).nx_neg_DST = num_neg_pts;
    end
end

function new_params = new_specification(freq, name, qincome, num_neg_pts, nbl_adj)
    new_params = setup.Params(freq, name, qincome);
    
    new_params.nx = 500 + num_neg_pts;
    new_params.nx_neg = num_neg_pts;
    new_params.nx_DST = 400 + num_neg_pts;
    new_params.nx_neg_DST = num_neg_pts;

    if num_neg_pts > 0
        new_params.nbl_adjustment = nbl_adj;
        new_params.borrow_lim = -1e10;
    end
end