function params = parameters_experiment(runopts,IncomeProcess)

    params = MPCParams(4,'test',IncomeProcess);
    params.set_grid(2000,2000,0.3);
    
    %----------------------------------------------------------------------
    % ADJUST PARAMETERS FOR FREQUENCY
    %----------------------------------------------------------------------
    
    params = MPCParams.adjust_if_quarterly(params);

    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    if runopts.Server==1 || runopts.Display==0
        params.set_display(0);
        params.set_simulate(0);
    else
        params.set_display(1);
        params.set_simulate(runopts.Simulate);
    end
    
    
    if runopts.fast == 1
        params.set_fast();
    end
    
    params.set_index();
end