function [params, all_names] = parameters_other(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    import aux.set_shared_fields

    % location of baseline income process for quarterly case
    quarterly_b_path = 'input/income_quarterly_b_contyT.mat';

    quarterly_b_params = struct();
    quarterly_b_params.sd_logyT = sqrt(0.6376);
    quarterly_b_params.lambdaT = 0.25;
    
    %----------------------------------------------------------------------
    % EXPERIMENTS
    %----------------------------------------------------------------------

    shocks = [-0.0081, -0.0405, -0.081, 0.0081, 0.0405, 0.081];

    % Quarterly
    params = setup.Params(4, 'target_assets_lt_1000_no_adj_costs', quarterly_b_path);
    params = set_shared_fields(params, quarterly_b_params);
    params.lumptransfer = 0.0081 * 2.0 * 4.0;
    params.beta0 = 0.867871450985079;
    params.shocks = shocks;
    params.nx = 100;
    params.nx_DST = 100;
    params.xgrid_par = 0.3;
    params.xmax = 50;
    params.gridspace_min = 0.0001;

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
    % if params.calibrate
    %     heterogeneity = setup.Prefs_R_Heterogeneity(params);
    %     new_betaH = params.betaH - max(heterogeneity.betagrid0);
    %     params.set("betaH", new_betaH, true);

    %     calibrator = aux.mean_wealth_calibrator(params);
    %     calibrator.set_handle(params);
    %     params.set("calibrator", calibrator, true);
    % end
end
