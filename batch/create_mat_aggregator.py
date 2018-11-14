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
    matfiles_new = []*len(mindex)

    for i in range(1,max(mindex)+1):
        # first use mindex = 1, then mindex = 2...
        findmindex = [j for j in range(len(mindex)) if mindex[j] == i]
        for ind in findmindex:
            matfiles_new.append(matfiles[ind])
    matfiles = matfiles_new
    
    with open(os.path.join(MWout,'aggregate.m'),'w') as newmfile:
        newmfile.write("clear\n% Aggregates spec##.mat's into one .mat file\n\n")
        for i,matfile in enumerate(matfiles):
            newmfile.write("matfiles{{{0:d}}} = '{1:s}';\n".format(i+1,matfile))
        
        lines = ['\nspec = 0;',
                 'for im = 1:numel(matfiles)',
                 '    S = load(matfiles{im});',
                 '    for is = 1:numel(S)',
                 '        spec = spec + 1;',
                 '        params(spec) = S.Sparams(is);',
                 '        results(spec) = S.results(is);',
                 '        decomps(spec) = S.decomps(is);',
                 '        checks(spec) = S.checks(is);',
                 '        exceptions(spec) = S.exceptions(is);',
                 '    end',
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
MWout = '/Users/Brian/Documents/midway2_output/discrete_time/11_14_18'

# ---------------------------------------------------------------------
# Function calls
# ---------------------------------------------------------------------
# create batch m-files
#(nb1end,nb5end) = gen_mfiles(mfile,args)

# create batch sbatch scripts
#gen_sbatch(mfile,args,nb1end,nb5end)

# generate m-file that aggreagates .mat files and creates a table
gen_mfile_aggregator(MWout)

