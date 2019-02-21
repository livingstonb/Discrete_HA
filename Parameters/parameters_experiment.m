function params = parameters_experiment(runopts,IncomeProcess)
    % Used to run experiments outside of batch mode

    ibw = 0.01;
    bs = 1/50;
    deathind = ' NoDeath'
    name = ['2 RandomBetaHet5 Width' num2str(ibw) ' SwitchProb' num2str(bs) deathind];
    params(1) = MPCParams(4,name,IncomeProcess);
    params(1).nb = 5;
    params(1).betawidth = ibw;
    params(1).betaswitch = bs;
    params(1).dieprob = 0;
   
    params(1).nx = 500;
    params(1).xgrid_par = 0.2;
    params(1).nxlong = 400;

    %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS
    %----------------------------------------------------------------------
    
    params = MPCParams.adjust_if_quarterly(params);
    params.set_run_parameters(runopts);
    params.set_index();
end