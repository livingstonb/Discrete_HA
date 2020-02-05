# Discrete Time Heterogenous Agent Model

The algorithms used to solve the model were in large part
originally written by Greg Kaplan and modified by Brian Livingston
(livingstonb@uchicago.edu).
The *aux_lib* directory contains code provided by Mario Miranda and Paul Fackler
via the CompEcon toolbox. In several other cases, this repository
contains code written by others, with citations where possible.

## Replicating the experiments in Fuster, Kaplan, and Zafar (2020)

The script *master_replication.m* is provided to solve for the calibrations
presented in the paper and compute statistics.

### The calibration

The desired calibration must be selected prior to running the code.
This is done in the *CHOOSE CALIBRATION* section in *master_replication.m*.
The parameters specific to each calibration can be viewed in the subsequent
sections of that script.
Note than any parameters not explicitly set in *master_replication.m*
will take their default values, which are declared in
*code/+setup/Params.m*.

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