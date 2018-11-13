function params = parameters_experiment(runopts,IncomeProcess)

    params = MPCParams(4,'test',IncomeProcess);
    params.betaswitch = 1/50;
    params.nb = 1;
    
    %----------------------------------------------------------------------
    % ADJUST PARAMETERS FOR FREQUENCY
    %----------------------------------------------------------------------
    
    params = MPCParams.adjust_if_quarterly(params);

    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    if runopts.Server==1 || runopts.Display==0
        params.set_display(0);
    else
        params.set_display(1);
    end
    
    if runopts.fast == 1
        params.set_fast();
    end
    
    params.set_index();
end