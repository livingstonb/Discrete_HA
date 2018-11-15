import os
import re

# This script 

# ---------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------
                    
def gen_mfile_aggregator(MWout):
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
    for i in range(1,max(mindex)+1):
        # first use mindex = 1, then mindex = 2...
        findmindex = [j for j in range(len(mindex)) if mindex[j] == i]
        for ind in findmindex:
            if ~findmindex:
                # run j must have failed, .mat file not found
                matfiles_new.append('')
            matfiles_new.append(matfiles[ind])
    matfiles = matfiles_new
    
    with open(os.path.join(MWout,'aggregate.m'),'w') as newmfile:
        newmfile.write("clear\n% Aggregates specification##.mat's into one .mat file\n\n")
        for i,matfile in enumerate(matfiles):
            newmfile.write("matfiles{{{0:d}}} = '{1:s}';\n".format(i+1,matfile))
        
        lines = ['',
                 'for im = 1:len(matfiles)',
                 '    if isempty(matfiles{im})',
                 '        continue',
                 '    end',
                 '',
                 '    S = load(matfiles{im});',
                 '    params(im) = S.Sparams(im);',
                 '    results(im) = S.results(im);',
                 '    decomps(im) = S.decomps(im);',
                 '    checks(im) = S.checks(im);',
                 '    exceptions(im) = S.exceptions(im);',
                 'end',
                 '',
                 "addpath('/Users/Brian/Documents/GitHub/MPCrecode/Output Functions');"
                 "addpath('/Users/Brian/Documents/GitHub/MPCrecode/Solution Functions');"
                 '',
                 'decomp2 = decomposition2(params,results);',
                 '',
                 '[T_annual,T_quarter] = '
                 'create_table(params,results,decomps,checks,exceptions,decomp2);']
        newmfile.write('\n'.join(lines))
            
        

# ---------------------------------------------------------------------
# Startup
# ---------------------------------------------------------------------

# location of .mat output files
MWout = '/Users/Brian/Documents/midway2_output/discrete_time/11_14_18_fixedincomebug/matlab_output'

# ---------------------------------------------------------------------
# Function calls
# ---------------------------------------------------------------------
# create batch m-files
#(nb1end,nb5end) = gen_mfiles(mfile,args)

# create batch sbatch scripts
#gen_sbatch(mfile,args,nb1end,nb5end)

# generate m-file that aggreagates .mat files and creates a table
gen_mfile_aggregator(MWout)

