function params = parameters_grid_tests3(runopts,IncomeProcess)
    % This function solves the model with many different grid parameters
    % Uses baseline quarterly specification
    
    %----------------------------------------------------------------------
    % DIFFERENT GRIDS
    %----------------------------------------------------------------------
    
    counter = 1;
    for xmax = [1000 500 7500 1250 1500]
        name = ['xmax',num2str(xmax)];
        params(counter) = MPCParams(4,name,IncomeProcess);
        params(counter).xmax = xmax;
        counter = counter + 1;
    end
    
   
    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    params = MPCParams.adjust_if_quarterly(params);
    params.set_run_parameters(runopts);

    % creates ordered 'index' field
    params.set_index();
    
    % alternative income processes
    for ip = 1:numel(params)
        if isempty(params(ip).IncomeProcess)
            params(ip).IncomeProcess = IncomeProcess;
        end
    end
end