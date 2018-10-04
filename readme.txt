
————————————————————————————————————————————————————————————————————————————
M-FILES
————————————————————————————————————————————————————————————————————————————

master.m - Set parameters and call main function file.
	
	Calls:		egp_AR1_IID_tax_recode.m

egp_AR1_IID_tax_recode.m - Main function file. Given a structure of parameters 
	from master.m, this script calls various other functions to generate
	the `sim_results’ and `direct_results’ structures and makes plots.
	
	Called by: 	master.m
	
	Calls: 		gen_income_variables.m
			iterate_beta.m
			solve_EGP.m
			solve_EGP_deterministic.m
			simulate.m
			direct_mpcs.m
			makeplots.m
			print_statistics.m

gen_income_variables.m - Function that creates the `income’ structure, which 
	contains various fields related to the distributions of gross and net income.

	Called by:	none
	
	Calls:		egp_AR1_IID_tax_recode.m

iterate_beta.m - Function used to help find the beta associated with the desired
	wealth/income ratio. This function runs the optimization in two steps:
	(1) pass the function solve_EGP() to fzero() with a small grid and 
	large convergence tolerance, and (2) once beta has been found, restrict
	the search interval to within 0.01 of beta, and pass the function
	solve_EGP() to fzero() with a larger grid and smaller convergence tolerance.

	Called by: 	egp_AR1_IID_tax_recode.m
	
	Calls:		passes solve_EGP() to fzero()

solve_EGP.m - This function performs the method of endogenous grid functions to find
	saving and consumption policy functions. Also calls find_stationary() to
	get the stationary distribution. This function returns the deviation of the
	wealth/income ratio from the target as `AYdiff’, and returns the structure
	`model’ which contains various important objects used elsewhere in the code.

	Called by:	egp_AR1_IID_tax_recode.m

	Calls:		find_stationary.m

find_stationary.m - This function finds the stationary distribution and transition matrix
	for the passed arguments. Also returns the policy functions interpolated onto
	the grid `xgridinput’. All outputs are associated with `xgridinput’. The argument
	`ergodic_method’ allows the user to specify an iterative procedure to find
	the ergodic distribution (ergodic_method=1) or a direct method (ergodic_method=2).
	The argument `ergodic_tol’ specifies the tolerance of the iterative procedure
	if ergodic_method is equal to 1.

	Called by: 	egp_AR1_IID_tax_recode.m
	
	Calls:		ergodicdist.m

	
solve_EGP_deterministic.m - This function performs the same role as solve_EGP, except
	for the model without income risk.

	Called by: 	egp_AR1_IID_tax_recode.m
	
	Calls: 		none

simulate.m - This function performs simulations based on the policy functions found in
	solve_EGP.m. Also computes MPCs by calling simulation_MPCs(). Returns the
	structure `sim_results’ and the column vector `assetmeans’, where the latter
	shows mean assets as a function of time (to check convergence).

	Called by:	egp_AR1_IID_tax_recode.m
	
	Calls:		simulation_MPCs.m

direct_mpcs.m - This function finds MPCs by computing expected consumption functions in
	response to an income shock. Returns one- and four-period MPCs distributed 
	according to basemodel.SSdist. Also returns average MPCs associated with these
	MPC distributions.

	Called by: 	egp_AR1_IID_tax_recode.m

	Calls:		none

direct_mpcs_determinisitic.m - This function finds MPCs by computing expected consumption 
	functions in response to an income shock. Returns one- and four-period MPCs distributed 
	according to norisk.SSdist. Also returns average MPCs associated with these
	MPC distributions, and returns the norisk model with new fields.

	Called by: 	egp_AR1_IID_tax_recode.m

	Calls:		none


makeplots.m - Makes plots accorded to the passed structures.

	Called by: 	egp_AR1_IID_tax_recode.m

	Calls:		none

print_statistics.m - Prints important statistics to the screen.

	Called by: 	egp_AR1_IID_tax_recode.m

	Calls:		none

————————————————————————————————————————————————————————————————————————————
STRUCTURES LOCAL TO egp_AR1_IID_tax_recode.m
————————————————————————————————————————————————————————————————————————————
p -			stores parameters
basemodel - 		stores objects from the model with income risk
norisk 	- 		stores objects from the model without income risk
simulations - 		stores simulation results associated with basemodel
sgrid - 		stores savings grid in various sizes
xgrid - 		stores cash-on-hand grid in various sizes
income - 		stores objects from the income process
prefs - 		stores objects related to preferences
direct_results - 	stores important results associated with direct methods
sim_results - 		stores important results associated with simulations