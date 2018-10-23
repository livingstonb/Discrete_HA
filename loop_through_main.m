function [T_annual,T_quarter] = loop_through_main(runopts,IncomeProcess,selection)
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

    direct_results = cell(1,Nparams); % Results from direct computations
    norisk_results = cell(1,Nparams); % Results from norisk model
    sim_results    = cell(1,Nparams); % Results from simulations
    exceptions     = cell(1,Nparams); % ME objects on any exceptions thrown
    checks         = cell(1,Nparams); % Information on failed sanity checks
    decomps        = cell(1,Nparams); 
    income         = cell(1,Nparams);

    if runopts.Batch == 0
        tic
        [income{1},SR,DR,NR,checks{1},decomps{1}] = main(params(1));
        toc
        direct_results{1}  = DR;
        norisk_results{1}  = NR;
        sim_results{1}     = SR;      
    else
        for ip = 1:Nparams
            tic
            disp(['Trying parameterization ' params(ip).name])
            try
                % Main function
                [income{ip},SR,DR,NR,checks{ip},decomps{ip}] = main(params(ip));
                direct_results{ip}  = DR;
                norisk_results{ip}  = NR;
                sim_results{ip}     = SR;
                exceptions{ip} = []; % main function completed
            catch ME
                checks{ip} = 'EXCEPTION_THROWN';
                exceptions{ip} = ME;
            end
            disp(['Finished parameterization ' params(ip).name])
                toc

            if runopts.Server == 1
            save(savematpath,'sim_results','direct_results','norisk_results',...
                                'checks','exceptions','params','decomps');
            end
        end
    end

    %% DECOMPOSITIONS - COMPARISONS WITH BASELINE
    decomp2 = decomposition2(params,direct_results,exceptions);

    %% CREATE TABLE

    [T_annual,T_quarter] = create_table(params,direct_results,decomps,checks,exceptions,decomp2);

    if runopts.Server == 1
        if ~isempty(T_annual)
            writetable(T_annual,savetablepath_annual,'WriteRowNames',true);
        end
        if ~isempty(T_quarter)
            writetable(T_quarter,savetablepath_quarterly,'WriteRowNames',true);
        end
    end
end
