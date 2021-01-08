# Discrete Time Heterogenous Agent Model

The algorithms used to solve the model were in large part
originally written by Greg Kaplan and modified by Brian Livingston
(livingstonb@uchicago.edu).
The *aux_lib* directory contains code provided by Mario Miranda and Paul Fackler
via the CompEcon toolbox. In several other cases, this repository
contains code written by others, with citations where possible.

## Using the master script

The code is ultimately executed from the script *master.m*.
From this script, first set options by assigning values to the *runopts* structure in the OPTIONS section.
The *mode* field must be set to the filename of the parameters script in the *code/+setup* directory excluding the file extension.

### The parameters
Within the default parameters script, a parameterization is assigned to each structure within a structure array. The easiest way to select a specific parameterization to run is to set the *number* field of runopts equal to the index of the desired parameterization within the structure array.
To run your desired parameterization, either modify *code/+setup/parameters.m* or create a new parameters file (see *code/+setup/parameters_template.m* for a stripped-down example of a parameters file).
All parameter defaults are set in the class file *code/+setup/Params.m* and any values set in the selected parameters file will override the default value of the given variable.
If the user needs to wrap the model in a non-linear solver to match user-specified moments, this should be done in the parameters file. See the next section for details.

### Calibration
This repository uses the class defined in *code/+solver/DHACalibrator.m* to assist with matching moments. The class is implemented by creating a DHACalibrator instance in the parameters file, which is then assigned to an attribute of the Params object returned by the parameters file.






DHACalibrator(params, calibration.variables,...
            calibration.target_names, calibration.target_values)

calibrations = struct();
    calibrations.variables = {'beta0'};
    calibrations.target_names = {'median_a'};
    calibrations.target_values = [scf.median_totw];





### Calibrating parameters to match the data

As discussed in the paper, we've selected parameters in part to match
statistics observed in the data (e.g. mean wealth equal to 3.2 times mean annual income).
The *master_replication.m* script partially omits our calibration routines, which
involved solving the model iteratively with a root-finding algorithm.
To recreate those routines, see *master.m* for a taste of how
this was done.

### Output

The output table can be obtained as the *results_table* variable
upon completion of the *master_replication.m* script.




% (4) If convergence fails, betaH0 and/or betaL may need to be adjusted.
% betaL is the lower bound picked for beta during iteration and
% betaH0 is the adjustment factor to the upper bound. The code will
% guess a theoretical upper bound, and then will add betaH0 to
% to that value.

% RUNNING ON THE SERVER: To run in batch on the server, use 
% code/batch/server.sbatch as a template. That script sends an array to SLURM 
% that runs all of the requested parameters in parameters.m. Make sure
% that the range of numbers in the slurm array match the number of 
% parameters in the parameters script. Output files
% are stored in the Output directory

% OUTPUT: Results are stored in the 'results' structure. Its 'direct' property
% contains results found from computing the stationary distribution
% using non-simulation numerical methods. The 'sim' property contains results
% found from simulation, if the option is turned on.