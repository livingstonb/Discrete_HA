function params = parameters_grid_tests3(runopts,IncomeProcess)
    % This function solves the model with many different grid parameters
    % Uses baseline quarterly specification
    
    %----------------------------------------------------------------------
    % DIFFERENT GRIDS
    %----------------------------------------------------------------------
    
%     counter = 1;
%     for xmax = [1000 500 7500 1250 1500]
%         name = ['xmax',num2str(xmax)];
%         params(counter) = Params(4,name,IncomeProcess);
%         params(counter).xmax = xmax;
%         counter = counter + 1;
%     end
%     

    counter = 1;
    for xmax = [200 1000]
        for curv = [0.15 0.2 0.25]
            name = ['curv',num2str(curv),' xmax',num2str(xmax)];
            params(counter) = setup.Params(4,name,IncomeProcess);
            params(counter).xmax = xmax;
            params(counter).xgrid_par = curv;
            counter = counter + 1;
        end
    end
    
   
    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS, DO NOT CHANGE
    %----------------------------------------------------------------------
    
    params = setup.Params.adjust_if_quarterly(params);
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