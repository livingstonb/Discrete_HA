function calibrator = mean_wealth_calibrator(p)
	import solver.Calibrator

	param_name = {'beta0'};
	stat_name = {'mean_a'};
	stat_target = p.target_value;

	calibrator = Calibrator(p, param_name,...
        stat_name, stat_target);

	beta_bounds = [p.betaL, p.betaH];
    calibrator.set_param_bounds(beta_bounds);
end