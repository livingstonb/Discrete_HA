function AYdiff = iterate_EGP(...
	x,p,grdEGP,grdDST,heterogeneity,income,mpcshock)
	% this function provides a function handle that can
	% be used to find beta iteratively using a routine
	% such as fzero

	% AYdiff is the absolute difference between the ratio
	% of mean wealth to annual income and the target for
	% for that ratio

	if p.EpsteinZin == 1
		egp_ez_solver = solver.EGP_EZ_Solver(x,p,grdEGP,heterogeneity,income);
		egp_ez_solver.solve(income);

		% get policy fns, etc...
		model = egp_ez_solver.return_model();
	else
		% get policy fns, etc...
		model = solver.solve_EGP(x,p,grdEGP,heterogeneity,...
							income,mpcshock,[]);
	end

	AYdiff = solver.find_stationary_adist(p,model,income,grdDST,heterogeneity);

end