function params = parameters_experiment(runopts)

    params = MPCParams(1,'test');
    params.nxlong = 2000;
    params.betaswitch = 1/50;
    params.nb = 1;
    
    params.Simulate = 1;
    
    %----------------------------------------------------------------------
    % ADJUST PARAMETERS FOR FREQUENCY
    %----------------------------------------------------------------------
    
    params = MPCParams.adjust_if_quarterly(params);

    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    if runopts.Server==1 || runopts.Display==0
        [params.Display] = deal(0);
    else
        [params.Display] = deal(1);
    end
    
    if runopts.fast == 1
        params.set_fast();
    end
    
    params.set_index();
end