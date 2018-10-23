
clear;
close all;

%% RUN OPTIONS
Batch = 0;
Server = 0;
Fast = 1;

% empty string if not loading from file
% IncomeProcess = 'IncomeVariables/quarterly_a.mat';
IncomeProcess = '';

% select only a subset of experiments
names_to_run = {}; % empty cell array to run all names
suffix = '';
Frequencies = [1 4]; % [1 4], 1, or 4

%% Add paths
if Server == 0
    path = '/Users/Brian/Documents/GitHub/MPCrecode';
    savetablepath_annual = '/Users/Brian/Documents/table_annual.xls';
    savetablepath_quarterly = '/Users/Brian/Documents/table_quarterly.xls';
    savematpath = '/Users/Brian/Documents/variables.mat';
else
    path = '/home/livingstonb/GitHub/MPCrecode';
    savetablepath_annual = ['/home/livingstonb/output/table_annual' suffix '.xls'];
    savetablepath_quarterly = ['/home/livingstonb/output/table_quarterly' suffix '.xls'];
    savematpath = ['/home/livingstonb/output/variables' suffix '.mat'];
end
addpath([path '/Auxiliary Functions']);
addpath([path '/MPC Functions']);
addpath([path '/Output Functions']);
addpath([path '/EGP']);
cd(path);

%% PARAMETERIZATIONS
if Batch == 0
	params = parameters_experiment(Fast); 
else
	params = parameters(Fast);
end

params = MPCParams.select_by_names(params,names_to_run);
if numel(Frequencies) == 1
    params = MPCParams.select_by_freq(params,Frequencies);
end

[params.IncomeProcess] = deal(IncomeProcess);
params.set_index();

%% CALL MAIN FUNCTION
Nparams = size(params,2);

direct_results = cell(1,Nparams); % Results from direct computations
norisk_results = cell(1,Nparams); % Results from norisk model
sim_results    = cell(1,Nparams); % Results from simulations
exceptions     = cell(1,Nparams); % ME objects on any exceptions thrown
checks         = cell(1,Nparams); % Information on failed sanity checks
decomps        = cell(1,Nparams); 
income         = cell(1,Nparams);

if Batch == 0
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
            
        if Server == 1
        save(savematpath,'sim_results','direct_results','norisk_results',...
                            'checks','exceptions','params','decomps');
        end
    end
end

%% DECOMPOSITIONS - COMPARISONS WITH BASELINE
decomp2 = decomposition2(params,direct_results,exceptions);

%% SAVE VARIABLES AND CREATE TABLE
                                             
[T_annual,T_quarter] = create_table(params,direct_results,decomps,checks,exceptions,decomp2);
                            
if Server == 1
    if ~isempty(T_annual)
        writetable(T_annual,savetablepath_annual,'WriteRowNames',true);
    end
    if ~isempty(T_quarter)
        writetable(T_quarter,savetablepath_quarterly,'WriteRowNames',true);
    end
    exit
end
