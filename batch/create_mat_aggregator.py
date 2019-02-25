import os
import re

# This script 

# ---------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------

                    
def gen_mfile_aggregator(MWout,sim):
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
                 "addpath('/home/brian/Documents/GitHub/Discrete_HA/Output Functions');",
                 "addpath('/home/brian/Documents/GitHub/Discrete_HA/Solution Functions');",
                 '',
                 decomp2,
                 '',
                 f'[T_annual,T_quarter] = {tablefn}']
        newmfile.write('\n'.join(lines))
            
        

# ---------------------------------------------------------------------
# Startup
# ---------------------------------------------------------------------

# location of .mat output files
MWout = '/home/brian/Desktop/temp/discrete_2_24_19'
# MWout = '/Users/brianlivingston/Documents/discrete_2_25_19'
sim = False # True/False

# ---------------------------------------------------------------------
# Function calls
# ---------------------------------------------------------------------

# generate m-file that aggreagates .mat files and creates a table
gen_mfile_aggregator(MWout,sim)

