
% basedir = '/home/livingstonb/GitHub/Discrete_HA';
basedir = '/Users/Brian-laptop/Documents/GitHub/Discrete_HA';
addpath([basedir '/code/Auxiliary_Functions']);
addpath([basedir '/code/Solve_Model']);
addpath([basedir '/code/Model_Setup']);

FROM_MATFILE = true;
if ~FROM_MATFILE
    clearvars -except params results decomps FROM_MATFILE
else
    clearvars -except FROM_MATFILE
end

if FROM_MATFILE
    % User must set basedir and date, where variablesX.mat files
    % are stored in <basedir>/<date>

    %% Read .mat files into a cell array
    fulldir = '/Users/Brian-laptop/Documents/midway2_output/discrete_time/8_9_19/';
%     fulldir = [basedir '/Output/'];
    
    decomp_baseline = cell(1,1);
    decomp3 = cell(1,1);
    
    ind = 0;
    for run = 1:999
        runstr = num2str(run);
        fpath = [fulldir,'variables',runstr,'.mat'];
        if exist(fpath,'file')
            ind = ind+1;
            
            S = load(fpath);
            params(ind) = S.Sparams;
            results(ind) = S.results;
            decomps(ind) = S.decomps;
            
            [decomp_baseline{ind},decomp_repagent{ind}] = baseline_repagent_decomps(params,results);
        else
            continue
        end
    end
    
    [T_annual,T_quarter] = create_table(...
    	params,results,decomps,decomp_baseline,decomp_repagent);
    
else
    
    [T_annual,T_quarter] = create_table(params,results,...
                                            decomps,[],[]);
end

writetable(T_quarter,[fulldir 'T_quarter.xlsx'],'WriteRowNames',true);
writetable(T_annual,[fulldir 'T_annual.xlsx'],'WriteRowNames',true);
