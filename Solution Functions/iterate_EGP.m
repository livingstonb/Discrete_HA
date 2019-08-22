function AYdiff = iterate_EGP(x,p,grdEGP,grdDST,heterogeneity,income,mpcshock)

	if p.EpsteinZin == 1
		model = solve_EGP_EZ(x,p,grdEGP,heterogeneity,income);
	else
		model = solve_EGP(x,p,grdEGP,heterogeneity,income,mpcshock,[]);
	end

	AYdiff = find_stationary_adist(p,model,income,heterogeneity,grdDST);

end