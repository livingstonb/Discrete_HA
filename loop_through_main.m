function [T_annual,T_quarter,results,checks,exceptions] = loop_through_main(runopts,IncomeProcess,selection)
    %% Add paths
    if runopts.Server == 0
        path = '/Users/Brian/Documents/GitHub/MPCrecode';
    else
        path = '/home/livingstonb/GitHub/MPCrecode';
        savetablepath_annual = ['/home/livingstonb/output/table_annual' selection.suffix '.xls'];
        savetablepath_quarterly = ['/home/livingstonb/output/table_quarterly' selection.suffix '.xls'];
        savematpath = ['/home/livingstonb/output/variables' selection.suffix '.mat'];
    end
    addpath([path '/Auxiliary Functions']);
    addpath([path '/MPC Functions']);
    addpath([path '/Output Functions']);
    addpath([path '/EGP']);
    cd(path);

    %% PARAMETERIZATIONS
    if runopts.Batch == 0
        params = parameters_experiment(runopts); 
    else
        params = parameters(runopts,selection);
    end

    [params.IncomeProcess] = deal(IncomeProcess);

    %% CALL MAIN FUNCTION
    Nparams = size(params,2);
    exceptions = cell(1,Nparams); % ME objects on any exceptions thrown
    checks     = cell(1,Nparams); % Information on failed sanity checks
    decomps    = cell(1,Nparams); 

    for ip = 1:Nparams
        tic
        disp(['Trying parameterization ' params(ip).name])
        if runopts.Server == 0
            [results(ip),checks{ip},decomps{ip}] = main(params(ip));
            exceptions{ip} = [];
        else
            try
                % Main function
                [results(ip),checks{ip},decomps{ip}] = main(params(ip));
                exceptions{ip} = []; % main function completed
            catch ME
                checks{ip} = 'EXCEPTION_THROWN';
                exceptions{ip} = ME;
            end
            save(savematpath,'results','checks','exceptions','params','decomps');         
        end
        disp(['Finished parameterization ' params(ip).name])
        toc
    end

    %% DECOMPOSITIONS - COMPARISONS WITH BASELINE
    decomp2 = decomposition2(params,results,exceptions);

    %% CREATE TABLE

    [T_annual,T_quarter] = create_table(params,results,...
                                    decomps,checks,exceptions,decomp2);

    if runopts.Server == 1
        if ~isempty(T_annual)
            writetable(T_annual,savetablepath_annual,'WriteRowNames',true);
        end
        if ~isempty(T_quarter)
            writetable(T_quarter,savetablepath_quarterly,'WriteRowNames',true);
        end
    end
end
