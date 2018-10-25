function [T_annual,T_quarter,results,checks,exceptions] = loop_through_main(runopts,IncomeProcess,selection)
    

    %% PARAMETERIZATIONS
    if runopts.Batch == 0
        params = parameters_experiment(runopts,IncomeProcess); 
    else
        params = parameters(runopts,selection,IncomeProcess);
    end

    %% CALL MAIN FUNCTION
    Nparams = size(params,2);
    exceptions = cell(1,Nparams); % ME objects on any exceptions thrown
    checks     = cell(1,Nparams); % Information on failed sanity checks
    decomps    = cell(1,Nparams); 

    for ip = 1:Nparams
        tic
        if params(ip).freq == 1
            msgfreq = 'annual';
        else
            msgfreq = 'quarterly';
        end
        fprintf('\n Trying %s parameterization "%s"\n',msgfreq,params(ip).name)
        if runopts.TryCatch == 1 || runopts.Server == 1
            try
                % Main function
                [results(ip),checks{ip},decomps{ip}] = main(params(ip));
                exceptions{ip} = []; % main function completed
            catch ME
                checks{ip} = 'EXCEPTION_THROWN';
                exceptions{ip} = ME;
                results(ip) = struct('direct',[],'norisk',[],'sim',[]);
            end
            
            if runopts.Server == 1
                save(runopts.savematpath,'results','checks','exceptions','params','decomps');    
            end
        else
            [results(ip),checks{ip},decomps{ip}] = main(params(ip));
            exceptions{ip} = [];
        end
        disp(['Finished parameterization ' params(ip).name])
        toc
    end

    %% DECOMPOSITIONS - COMPARISONS WITH BASELINE
    for ip = 1:numel(params)
        decomp2(ip) = DecompTwo(params(ip));
    end
    decomp2 = DecompTwo.decompose(decomp2,params,results);

    %% CREATE TABLE

    [T_annual,T_quarter] = create_table(params,results,...
                                    decomps,checks,exceptions,decomp2);

    if runopts.Server == 1
        if ~isempty(T_annual)
            writetable(T_annual,runopts.savetablepath_annual,'WriteRowNames',true);
        end
        if ~isempty(T_quarter)
            writetable(T_quarter,runopts.savetablepath_quarterly,'WriteRowNames',true);
        end
    end
end
