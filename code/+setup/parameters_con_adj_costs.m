function [params, all_names] = parameters_con_adj_costs(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    import aux.set_shared_fields

    shared_params = struct();
    shared_params.annual_inc_dollars = 61728;
    shared_params.shocks = [-500, -2500, -5000, 500, 2500, 5000] ...
        ./ shared_params.annual_inc_dollars;

    % location of baseline income process for quarterly case
    quarterly_b_path = 'input/income_quarterly_b_contyT.mat';

    quarterly_b_params = struct();
    quarterly_b_params.sd_logyT = sqrt(0.6376);
    quarterly_b_params.lambdaT = 0.25;
    quarterly_b_params.nyT = 11;
    
    %----------------------------------------------------------------------
    % EXPERIMENTS
    %----------------------------------------------------------------------

    % Quarterly
    params = setup.Params(4, 'baseline', quarterly_b_path);
    params = set_shared_fields(params, quarterly_b_params);
    params = set_shared_fields(params, shared_params);
    params.lumptransfer = 0.0081 * 2.0 * 4.0;
    params.beta0 = 0.867871450985079;
    params.nx = 150;
    params.nx_DST = 150;
    params.xgrid_par = 0.3;
    params.xmax = 20;
    params.target_value = 0.23;
    params.ResetIncomeUponDeath = false;

    % params = setup.Params(4,'wealth3.2',QIncome);
    % params.targetAY = 3.2;
    % params.lumptransfer = 0.0081 * 2.0 * 4.0;
    % params.shocks = shocks;

    % params(2) = setup.Params(4,'wealth0.3',QIncome);
    % params(2).targetAY = 0.3;
    % params(2).lumptransfer = 0.0081 * 2.0 * 4.0;
    % params(2).shocks = shocks;

    % ii = 3;
    % for bwidth = [0.0289]
    % 	name = sprintf('beta_heterog_width%1.4f', bwidth);
    % 	params(ii) = setup.Params(4,name,QIncome);
    % 	params(ii).betawidth = bwidth;
    %     params(ii).nbeta = 3;
    %     params(ii).targetAY = 3.2;
    %     params(ii).lumptransfer = 0.0081 * 2.0 * 4.0;
    %     params(ii).shocks = shocks;
    %     params(ii).beta0 = 0.82;
    % 	ii = ii + 1;
    % end
    

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

        calibrator = lt_1000_calibrator(params);
        calibrator.set_handle(params);
        params.set("calibrator", calibrator, true);
    end
end

function calibrator = lt_1000_calibrator(p)
    import solver.Calibrator

    param_name = {'beta0'};
    stat_name = {'wealth_lt_1000'};
    stat_target = p.target_value;

    calibrator = Calibrator(p, param_name,...
        stat_name, stat_target);

    beta_bounds = [p.betaL, p.betaH];
    calibrator.set_param_bounds(beta_bounds);
end