function params = parameters_experiment(runopts,IncomeProcess)
    % Used to run experiments outside of batch mode

    params = MPCParams(4,'test',IncomeProcess);
    
    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    params = MPCParams.adjust_if_quarterly(params);
    params.set_run_parameters(runopts);
    params.set_index();
end