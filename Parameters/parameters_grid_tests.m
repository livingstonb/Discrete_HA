function params = parameters_grid_tests(runopts,selection,IncomeProcess)
    % This function solves the model with many different grid parameters
    
    %----------------------------------------------------------------------
    % DIFFERENT GRIDS
    %----------------------------------------------------------------------
    
    counter = 0;
    for nx = [100 150 200]
    for nxlong = [100 200 400 500]
    for curv = [0.2 0.3 0.4]
        counter = counter + 1;
        name = ['nx',num2str(nx),'_nxlong',num2str(nxlong),'_curv',num2str(curv)];
        params(counter) = MPCParams(4,name,IncomeProcess);
        params(counter).set_grid(nx,nxlong,curv);
    end
    end
    end
    
    %----------------------------------------------------------------------
    % ADJUST PARAMETERS FOR FREQUENCY
    %----------------------------------------------------------------------
    
    params = MPCParams.adjust_if_quarterly(params);
    
    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    params.set_display(runopts.Display);
    params.set_simulate(runopts.Simulate);
    
    params.set_index();
    
    if numel(selection.number) == 1
        params = MPCParams.select_by_number(params,selection.number);
        % index by individual .mat filerunopts
    elseif numel(selection.number) > 1
        error('selection.number must have 1 or zero elements')
    else
        params = MPCParams.select_by_names(params,selection.names_to_run);
        params.set_index(); % index within .mat file
    end
    
    for ip = 1:numel(params)
        if isempty(params(ip).IncomeProcess)
            params(ip).IncomeProcess = IncomeProcess;
        end
    end
end