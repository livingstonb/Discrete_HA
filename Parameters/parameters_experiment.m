function params = parameters_experiment(runopts,IncomeProcess)
    % Used to run experiments outside of batch mode

    params = MPCParams(4,'test',IncomeProcess);
    params.set_grid(2000,2000,0.3);
    
    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    params = MPCParams.adjust_if_quarterly(params);
    params.set_run_parameters(runopts);
    params.set_index();
end