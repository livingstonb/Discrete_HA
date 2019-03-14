import os
import re

# This script takes .mat files labeled variablesXX.mat with XX any number
# greather than or equal to 1, and aggregates them into a matlab table

# user should only need to enter options into 'Startup' section

# ---------------------------------------------------------------------
# Startup
# ---------------------------------------------------------------------

# location of Discrete_HA code directory
repos = '/home/brian/Documents/GitHub/Discrete_HA'

# location of variablesXX.mat files
MWout = '/media/hdd/Other/midway2_output/discrete_time/3_12_19'

# option for using table based on simulations if available
# (may not currently work)
sim = False 

# ---------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------

                    
def gen_mfile_aggregator(MWout,sim,repos):
    # create m-file to aggregate .mat files and create a table

    matdir = MWout
    matfiles = os.listdir(matdir)
    pattern = re.compile('variables([0-9]+).mat')
    matfiles = list(filter(lambda x: re.fullmatch(pattern,x),matfiles))
    
    matfiles = sorted(matfiles)
    # pull out index of each matfile
    mindex = list(map(lambda x: re.fullmatch(pattern,x).group(1),matfiles))
    mindex = list(map(lambda x: int(x),mindex))
    matfiles = list(map(lambda x: os.path.join(matdir,x),matfiles))
    matfiles_new = []

    # sort by mindex
    matfiles = sorted(matfiles,key=lambda x:mindex[matfiles.index(x)])
    
    output_fns_dir = repos + '/Output Functions'
    soln_fns_dir = repos + '/Solution Functions'

    with open(os.path.join(MWout,'aggregate.m'),'w') as newmfile:
        newmfile.write("clear\n% Aggregates specification##.mat's into one .mat file\n\n")
        for i,matfile in enumerate(matfiles):
            newmfile.write("matfiles{{{0:d}}} = '{1:s}';\n".format(i+1,matfile))
        
        if sim:
            tablefn = 'create_table_sim(params,results,decomps,checks);'
            decomp2 = ''
        else:
            tablefn = 'create_table(params,results,decomps,checks,decomp2,decomp3);'
            decomp2 = '[decomp2,decomp3] = decomposition2(params,results);'

        lines = ['',
                 'for im = 1:length(matfiles)',
                 '    if isempty(matfiles{im})',
                 '        continue',
                 '    end',
                 '',
                 '    S = load(matfiles{im});',
                 '    params(im) = S.Sparams;',
                 '    results(im) = S.results;',
                 '    decomps(im) = S.decomps;',
                 '    checks(im) = S.checks;',
                 'end',
                 '',
                 f"addpath('{output_fns_dir}');",
                 f"addpath('{soln_fns_dir}');",
                 '',
                 decomp2,
                 '',
                 f'[T_annual,T_quarter] = {tablefn}']
        newmfile.write('\n'.join(lines))
            
 
# ---------------------------------------------------------------------
# Function calls
# ---------------------------------------------------------------------

# generate m-file that aggreagates .mat files and creates a table
gen_mfile_aggregator(MWout,sim,repos)

