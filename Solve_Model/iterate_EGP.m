function AYdiff = iterate_EGP(x,p,grdEGP,grdDST,heterogeneity,income,mpcshock)

	if p.EpsteinZin == 1
		egp_ez_solver = EGP_EZ_Solver(x,p,grdEGP,heterogeneity,income);
		egp_ez_solver.solve(income);
		model = egp_ez_solver.return_model();
	else
		model = solve_EGP(x,p,grdEGP,heterogeneity,income,mpcshock,[]);
	end

	AYdiff = find_stationary_adist(p,model,income,grdDST);

end