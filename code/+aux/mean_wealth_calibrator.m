function [fn_handle, x0] = mean_wealth_calibrator(p)
	import solver.Calibrator

	param_name = 'beta0';
	stat_name = 'mean_a';
	stat_target = p.target_value;
	beta_init = p.beta0;

	calibrator = Calibrator(p, param_name,...
        stat_name, stat_target);

	beta_bounds = [p.betaL, p.betaH];
    calibrator.set_param_bounds(beta_bounds);

    x0 = calibrator.convert_to_solver_input(beta_init);
    fn_handle = @(x) calibrator.fn_handle(x, p);
end