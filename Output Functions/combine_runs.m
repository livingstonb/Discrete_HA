clearvars -except params results decomps checks

FROM_MATFILE = false;

if FROM_MATFILE
    % User must set basedir and date, where variablesX.mat files
    % are stored in <basedir>/<date>

    %% Select directories
    basedir = '/media/hdd/Other/midway2_output/continuous_time';

    %% Read .mat files into a cell array
    fulldir = [basedir '/' date '/'];
    
    params = struct();
    results = struct();
    decomps = struct();
    decomp2 = struct();
    decomp3 = struct();
    
    ind = 0;
    for run = 1:999
        runstr = num2str(run);
        fpath = [fulldir,'output_',runstr,'.mat'];
        if exist(fpath,'file')
            ind = ind+1;
            
            S = load(fpath);
            params(ind) = S.Sparams;
            results(ind) = S.results;
            decomps(ind) = S.decomps;
            
            [decomp2(ind),decomp3(ind)] = decomposition2(params,results);
        else
            continue
        end
    end
    
else
    
    NoDecomps = true;
    [T_annual,T_quarter] = create_table(params,results,...
                                            decomps,NoDecomps,[],[]);
end