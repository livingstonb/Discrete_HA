function [params, all_names] = parameters_con_adj_costs2(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    import aux.set_shared_fields
    all_names = {};

    shared_params = struct();
    shared_params.annual_inc_dollars = 61728;
    
    dollars = [-500, -2500, -5000, 500, 2500, 5000];
    shared_params.shocks = dollars ./ shared_params.annual_inc_dollars;
    shared_params.lumptransfer = 0.0081 * 2.0 * 4;
    shared_params.ResetIncomeUponDeath = false;
    shared_params.target_value = 0.210464;
    
    shared_params.xgrid_par = 0.1;
    shared_params.xgrid_term1wt = 0.02;
    shared_params.xgrid_term1curv = 0.9;
    
    shared_params.shocks_labels = {};
    for ishock = 1:6
        val = dollars(ishock);
        if val < 0
            shared_params.shocks_labels{ishock} = sprintf('-$%g', abs(val));
        else
            shared_params.shocks_labels{ishock} = sprintf('$%g', abs(val));
        end
    end

    % location of baseline income process for quarterly case
    quarterly_b_path = 'input/income_con_adj_costs.mat';

    params = setup.Params(4, 'Quarterly', quarterly_b_path);
    params = set_shared_fields(params, shared_params);
    params.beta0 = 0.984363510593659;
    params.target_value = 3.2;
    
    params.set_run_parameters(runopts);
    params.make_adjustments();

    %----------------------------------------------------------------------
    % ATTACH CALIBRATOR
    %----------------------------------------------------------------------
    if params.calibrate
        heterogeneity = setup.Prefs_R_Heterogeneity(params);

        if (params.nbeta > 1) && isequal(heterogeneity.ztrans, eye(params.nbeta))
            new_betaH = params.betaH - max(heterogeneity.betagrid0);
            params.set("betaH", new_betaH, true);
        end

        calibrator = aux.mean_wealth_calibrator(params);
        calibrator.set_handle(params);
        params.set("calibrator", calibrator, true);
    end
end