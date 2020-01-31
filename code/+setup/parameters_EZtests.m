function params = parameters_EZtests(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    % location of baseline income process for quarterly case
    QIncome = 'input/quarterly_b.mat';
    
    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    
%     % Annual
%     params(1) = setup.Params(1,'baseline_A','');
%     
    % Quarterly
    params(1) = setup.Params(4,'baseline_Q',QIncome);
    params(end).beta0 = 0.984363421;
    
    %----------------------------------------------------------------------
    % EPSTEIN ZIN
    %----------------------------------------------------------------------
    
    % EZ baseline
    params(end+1) = setup.Params(4, 'EZ_baseline',QIncome);
    params(end).EpsteinZin = 1;
    params(end).invies = 1;
    params(end).risk_aver = 1;
    params(end).beta0 = 0.984363421;
    
    % Higher risk aversion, holding IES constant
    % turn off beta iteration!
    params(end+1) = setup.Params(4, 'EZ_riskaver2',QIncome);
    params(end).EpsteinZin = 1;
    params(end).invies = 1;
    params(end).risk_aver = 2;
    params(end).beta0 = 0.984363421;
    
    % Even higher risk aversion, holding IES constant
    % turn off beta iteration!
    params(end+1) = setup.Params(4, 'EZ_riskaver4',QIncome);
    params(end).EpsteinZin = 1;
    params(end).invies = 1;
    params(end).risk_aver = 4;
    params(end).beta0 = 0.984363421;
    
    % Higher invies, holding risk aversion constant
    % turn off beta iteration!
    params(end+1) = setup.Params(4, 'EZ_invies2',QIncome);
    params(end).EpsteinZin = 1;
    params(end).invies = 2;
    params(end).risk_aver = 1;
    params(end).beta0 = 0.984363421;
    
    % Even higher invies, holding risk aversion constant
    % turn off beta iteration!
    params(end+1) = setup.Params(4, 'EZ_invies4',QIncome);
    params(end).EpsteinZin = 1;
    params(end).invies = 4;
    params(end).risk_aver = 1;
    params(end).beta0 = 0.984363421;
    
    % Higher invies and risk aversion
    % turn off beta iteration!
    params(end+1) = setup.Params(4, 'EZ_invies2riskaver2',QIncome);
    params(end).EpsteinZin = 1;
    params(end).invies = 2;
    params(end).risk_aver = 2;
    params(end).beta0 = 0.984363421;
    
    % Even higher invies and risk aversion
    % turn off beta iteration!
    params(end+1) = setup.Params(4, 'EZ_invies2riskaver2',QIncome);
    params(end).EpsteinZin = 1;
    params(end).invies = 4;
    params(end).risk_aver = 4;
    params(end).beta0 = 0.984363421;
    
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
    
    % select by number if there is one, otherwise select by names,
    % otherwise use all
    if numel(runopts.number) == 1
        params = setup.Params.select_by_number(params,runopts.number);
    elseif numel(runopts.number) > 1
        error('runopts.number must have 1 or zero elements')
    else
        params = setup.Params.select_by_names(params,runopts.names_to_run);
        params.set_index(); % index within .mat file
    end