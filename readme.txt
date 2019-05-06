---------------------------------------------------------------------------------
master.m
---------------------------------------------------------------------------------
This is the main script. Set options in runopts, and set Simulate to 1 to if you 
want to run simulations in addition to direct methods. QIncome locates the income 
income process files, although the code defaults to generating the income processes 
automatically.

To run all parameterizations in get_params, leave selection.names_to_run = {}. 
Otherwise, enter one or more strings in a cell array containing the names of the 
parameterizations you would like to run. The results for the i-th parameterization 
that you run will be stored in results(i).

---------------------------------------------------------------------------------
Parameters
---------------------------------------------------------------------------------
The baseline parameters are stored in the class file MPCParams.m. The numbers are 
meant to be changed, the rest of the code should be left the same.

To create parameterizations, see the parameters.m file. Each parameterization is 
created by instantiating the MPCParams class with the syntax 
MPCParams(frequency,name,incomeprocess), where frequency is 1 or 4, name is a string,
and income process is a string locating the desired income process files (using
the empty string generates the process via code).

If you create an alternate parameters.m script, you must include the two bottom sections
titled "ADJUST TO QUARTERLY VALUES" and "CALL METHODS/CHAGNE SELECTED PARAMETERS"
and they must be at the bottom of the script.

---------------------------------------------------------------------------------
Income process
---------------------------------------------------------------------------------
The income process is either loaded from or created within gen_income_variables.m.
The income variables are stored in the 'income' structure.

---------------------------------------------------------------------------------
Variable notes
---------------------------------------------------------------------------------
'grdHJB' contains asset and cash-on-hand grids for solving the HJB. The grid length
is nx. 'grdKFE' contains grids for solving the KFE. The grid length is nx_KFE.

---------------------------------------------------------------------------------
Warning about beta iteration
---------------------------------------------------------------------------------
Sometimes the lower and/or upper bound for beta do not work. Depending on the model,
you may have to change these, and you can do so by changing betaH0 and betaL
in parameters.m.

---------------------------------------------------------------------------------
Viewing results
---------------------------------------------------------------------------------
Results are stored in the 'results' structure, which contains the fields 'direct', 
'norisk', and 'sim', which store results for the direct methods, simulation, and 
the model with no income risk.