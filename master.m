clear;
close all;

%% RUN OPTIONS
runopts.Batch = 1; % use parameters.m, not parameters_experiment.m
runopts.Display = 1;
runopts.Server = 0; % use server paths and limit display
runopts.TryCatch = 0; % use try-catch block in main loop (auto-on if Server=1)
runopts.fast = 1; % specify very small asset and income grids for speed
runopts.localdir = '/Users/Brian/Documents/GitHub/MPCrecode';

% empty string if not loading from file
IncomeProcess = 'IncomeVariables/quarterly_b.mat';

% select only a subset of experiments
selection.names_to_run = {}; % cell array of strings, {} to run all
selection.suffix = ''; % string, added to filenames
selection.frequencies = [4]; % 1,4 ([] to run all)
selection.nb = [5]; % 1,5 ([] to run all)

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
addpath([runopts.path '/MPC Functions']);
addpath([runopts.path '/Output Functions']);
addpath([runopts.path '/EGP']);
cd(runopts.path);

%% PARAMETERIZATIONS
if runopts.Batch == 0
    params = parameters_experiment(runopts); 
else
    params = parameters(runopts,selection);
end

% Set income process for specifications where it has not been set
unset = find(cellfun(@(x) isempty(x),{params.IncomeProcess}));
[params(unset).IncomeProcess] = deal(IncomeProcess);

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

%% CREATE TABLE/SAVE

if runopts.Server == 0
    [T_annual,T_quarter] = create_table(params,results,...
                                    decomps,checks,exceptions,decomp2);
else
    save(runopts.savematpath,'params','results','decomps','checks','exceptions','decomp2')
    exit
end