function params = parameters_grid_tests(runopts,IncomeProcess)
    % This function solves the model with many different grid parameters
    % Uses baseline quarterly specification
    
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