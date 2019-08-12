

FROM_MATFILE = true;
if ~FROM_MATFILE
    clearvars -except params results decomps FROM_MATFILE
else
    clearvars -except FROM_MATFILE
end

if FROM_MATFILE
    % User must set basedir and date, where variablesX.mat files
    % are stored in <basedir>/<date>

    %% Select directories
    basedir = '/Users/Brian-laptop/Documents/midway2_output/discrete_time';
    date = '8_9_19';

    %% Read .mat files into a cell array
    fulldir = [basedir '/' date '/'];
    
%     params = struct();
%     results = struct();
%     decomps = struct();
%     decomp2 = struct();
%     decomp3 = struct();
decomp2 = cell(1,1);
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
            
            [decomp2,decomp3] = decomposition2(params,results);
        else
            continue
        end
    end
    
    [T_annual,T_quarter] = create_table(params,results,decomps,decomp2,decomp3);
    
else
    
    [T_annual,T_quarter] = create_table(params,results,...
                                            decomps,[],[]);
end