
————————————————————————————————————————————————————————————————————————————
M-FILES
————————————————————————————————————————————————————————————————————————————

master.m - Set parameters and call main function file.
	
	Calls:		main.m
			parameters.m

parameters.m - Creates structure array containing all of our desired
	parameterizations. Not yet complete.

main.m - Main function file. Given a structure of parameters 
	from master.m, this script calls various other functions to generate
	the `sim_results’ and `direct_results’ structures and makes plots.
	
	Called by: 	master.m
	
	Calls: 		gen_income_variables.m
			solve_EGP.m
			solve_EGP_deterministic.m
			simulate.m
			direct_MPCs.m
			direct_MPCs_deterministic.m
			makeplots.m
			print_statistics.m

gen_income_variables.m - Function that creates the `income’ structure, which 
	contains various fields related to the distributions of gross and net income.

	Called by:	none
	
	Calls:		main.m


solve_EGP.m - This function performs the method of endogenous grid functions to find
	saving and consumption policy functions. Also calls find_stationary() to
	get the stationary distribution. This function returns the deviation of the
	wealth/income ratio from the target as `AYdiff’, and returns the structure
	`model’ which contains various important objects used elsewhere in the code.

	Called by:	main.m

	Calls:		find_stationary.m

find_stationary.m - This function finds the stationary distribution
	for the passed arguments. Also returns the policy functions interpolated onto
	the grid `xgridinput’. All outputs are associated with `xgridinput’.

	Called by: 	main.m
	
	Calls:		eigs()

	
solve_EGP_deterministic.m - This function performs the same role as solve_EGP, except
	for the model without income risk.

	Called by: 	main.m
	
	Calls: 		none

simulate.m - This function performs simulations based on the policy functions found in
	solve_EGP.m. Also computes MPCs by calling simulation_MPCs(). Returns the
	structure `sim_results’ and the column vector `assetmeans’, where the latter
	shows mean assets as a function of time (to check convergence).

	Called by:	main.m
	
	Calls:		simulation_MPCs.m

simulation_MPCs.m - Computes MPCs from the simulated data generated in simulate().

	Called by:	simulate.m
	
	Calls:		none

direct_MPCs.m - Computes MPCs by drawing from the stationary distribution of assets (not
	cash-on-hand) and simulating 1-4 periods.

	Called by: 	main.m

	Calls:  	none

direct_MPCs_deterministic.m - Computes MPCs by using the sample for initial assets found in
	direct_MPCs(), via simulating 1-4 periods.

	Called by: 	main.m

	Calls:  	none

makeplots.m - Makes plots accorded to the passed structures.

	Called by: 	main.m

	Calls:		none

print_statistics.m - Prints important statistics to the screen.

	Called by: 	main.m

	Calls:		none

————————————————————————————————————————————————————————————————————————————
STRUCTURES LOCAL TO main.m
————————————————————————————————————————————————————————————————————————————
p -			stores parameters
basemodel - 		stores objects from the model with income risk
norisk 	- 		stores objects from the model without income risk
simulations - 		stores simulation objects associated with basemodel
sgrid - 		stores savings grid in various sizes
xgrid - 		stores cash-on-hand grid in various sizes
income - 		stores objects from the income process
prefs - 		stores objects related to preferences
direct_results - 	stores results associated with direct methods for basemodel
sim_results - 		stores results associated with simulations for basemodel
norisk_results 		stores results associated with direct methods for norisk model