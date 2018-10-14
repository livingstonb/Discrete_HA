function objective = iterate_bequests(beqweight,params,target)
    newparams = params;
    newparams.bequest_weight = beqweight;
   
    % Call main
    [~,direct_results] = main(newparams);
    fprintf('\n\n\n MEAN BEQUESTS =  %6.4f\n\n\n',direct_results.mean_bequests);
    fprintf('\n\n\n BEQUEST WEIGHT =  %6.4f\n\n\n',beqweight);
    objective = direct_results.mean_bequests - target;

end