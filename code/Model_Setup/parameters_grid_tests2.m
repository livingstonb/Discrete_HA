function params = parameters_grid_tests2(runopts,IncomeProcess)
    % This function solves the model with many different grid parameters
    % Uses baseline quarterly specification
    
    %----------------------------------------------------------------------
    % DIFFERENT GRIDS
    %----------------------------------------------------------------------
    runopts.Simulate = 1;
    
    nxlong = 10;
    counter = 0;
    for nx = [200 500 1000]
    for curv = [0.2]
        counter = counter + 1;
        name = ['nx',num2str(nx),'_curv',num2str(curv)];
        params(counter) = Params(4,name,IncomeProcess);
        params(counter).set_grid(nx,nxlong,curv);
    end
    end
    
    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    params = Params.adjust_if_quarterly(params);
    params.set_run_parameters(runopts);

    % creates ordered 'index' field
    params.set_index();

    % select by number if there is one, otherwise select by names,
    % otherwise use all
    if numel(selection.number) == 1
        params = Params.select_by_number(params,runopts.number);
    elseif numel(runopts.number) > 1
        error('selection.number must have 1 or zero elements')
    else
        params = Params.select_by_names(params,runopts.names_to_run);
        params.set_index(); % index within .mat file
    end
    
    % alternative income processes
    for ip = 1:numel(params)
        if isempty(params(ip).IncomeProcess)
            params(ip).IncomeProcess = IncomeProcess;
        end
    end
end