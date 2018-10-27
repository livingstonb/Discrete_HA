clear;
close all;

%% RUN OPTIONS
runopts.Batch = 1; % use parameters.m, not parameters_experiment.m
runopts.Display = 1;
runopts.Server = 1; % use server paths and limit display
runopts.TryCatch = 0; % use try-catch block in main loop (auto-on if Server=1)
runopts.fast = 0; % specify very small asset and income grids for speed
runopts.localdir = '/Users/Brian/Documents/GitHub/MPCrecode';

% empty string if not loading from file
IncomeProcess = 'IncomeGrids/quarterly_b.mat';

% select only a subset of experiments
selection.names_to_run = {}; % cell array of strings, {} to run all
selection.suffix = ''; % string, added to filenames
selection.frequencies = []; % 1,4 ([] to run all)
selection.nb = []; % 1,5 ([] to run all)

%% Add paths
if runopts.Server == 0
    runopts.path = runopts.localdir;
else
    runopts.path = '/home/livingstonb/GitHub/MPCrecode';
    runopts.savetablepath_annual = ['/home/livingstonb/output/table_annual' selection.suffix '.xls'];
    runopts.savetablepath_quarterly = ['/home/livingstonb/output/table_quarterly' selection.suffix '.xls'];
    runopts.savematpath = ['/home/livingstonb/output/variables' selection.suffix '.mat'];
end
addpath([runopts.path '/Classes']);
addpath([runopts.path '/Auxiliary Functions']);
addpath([runopts.path '/Solution Functions']);
addpath([runopts.path '/Output Functions']);
addpath([runopts.path '/Parameters']);
cd(runopts.path);

%% PARAMETERIZATIONS
if runopts.Batch == 0
    params = parameters_experiment(runopts,IncomeProcess); 
else
    params = parameters(runopts,selection,IncomeProcess);
end

Sparams = MPCParams.to_struct(params);

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
            save(runopts.savematpath,'results','checks','exceptions','Sparams','decomps');    
        end
    else
        [results(ip),checks{ip},decomps{ip}] = main(params(ip));
        exceptions{ip} = [];
    end
    disp(['Finished parameterization ' params(ip).name])
    toc
end

%% DECOMPOSITIONS - COMPARISONS WITH BASELINE
decomp2 = decomposition2(params,results);

%% CREATE TABLE/SAVE

if runopts.Server == 0
    [T_annual,T_quarter] = create_table(params,results,...
                                    decomps,checks,exceptions,decomp2);
else
    save(runopts.savematpath,'Sparams','results','decomps','checks','exceptions','decomp2')
    exit
end
