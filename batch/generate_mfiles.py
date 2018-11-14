import os
import re

# This script 

# ---------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------

def gen_mfiles(mfile,args):
    # generate mfiles named matlab_nb#_#.m which are copies of master.m
    # with tweaks to run on the server in batch

    # delete any existing batch m-files
    deletepath = os.path.join(args['batchpath'],'mfiles')
    for dirpath, dirs, files in os.walk(deletepath):
        for nfile in files:
            os.remove(os.path.join(dirpath,nfile))
            
 
    # loop over specifications, identified by '_spec<#>' in files
    for nb in ['1','5']:
        count = 0
        for freq in args['frequencies']:
            for name in args['names'][nb]:
                count += 1
                label = '_nb' + nb + '_' + str(count)
                newfpath = os.path.join(args['batchpath'],'mfiles/master'+label+'.m')
                with open(newfpath,'w') as newmfile:
                    for line in mfile.readlines():
                        if line.startswith('runopts.Batch ='):
                            newmfile.write("runopts.Batch = 1;\n")
                        elif line.startswith('runopts.Server ='):
                            newmfile.write('runopts.Server = 1;\n')
                        elif line.startswith('runopts.fast ='):
                            newmfile.write('runopts.fast = '+args['fast']+';\n')
                        elif line.startswith('IncomeProcess'):
                            newmfile.write('IncomeProcess = '+args['Qincvar']+';\n')
                        elif line.startswith('selection.names_to_run'):
                            newmfile.write('selection.names_to_run = '+name+';\n')
                        elif line.startswith('selection.suffix'):
                            newmfile.write("selection.suffix = '"+label+"';\n")
                        elif line.startswith('selection.frequencies'):
                            newmfile.write('selection.frequencies = '+freq+';\n')
                        elif line.startswith('selection.nb'):
                            newmfile.write('selection.nb = '+nb+';\n')
                        else:
                            newmfile.write(line)
                    mfile.seek(0)
        if nb == '1':
            nb1end = count
        elif nb == '5':
            nb5end = count
    return (nb1end,nb5end)
            
def gen_sbatch(mfile,args,nb1end,nb5end):
    # generate sbatch scripts to run the m-files generated in the previous fcn

    # delete existing sbatch scripts
    deletepath = os.path.join(args['batchpath'],'sbatch')
    for dirpath, dirs, files in os.walk(deletepath):
        for nfile in files:
            os.remove(os.path.join(dirpath,nfile))

    for nb in [1,5]: 
        label = '_nb' + str(nb)
        if nb == 5:
            mem = '4000'
            end = nb5end
        elif nb == 1:
            mem = '2000'
            end = nb1end

        with open(os.path.join(args['batchpath'],'sbatch','s'+label+'.sbatch'),'w') as newmfile:
            lines = ['#!/bin/bash',
                     '#SBATCH --job-name=disc' + label,
                     '#SBATCH --output='+'/home/livingstonb/output/matlab'+label+'_%a.out',
                     '#SBATCH --error=' +'/home/livingstonb/output/matlab'+label+'_%a.err',
                     '#SBATCH --partition=broadwl',
                     '#SBATCH --time='+args['time'],
                     '#SBATCH --array=1-'+str(end),
                     '#SBATCH --nodes=1',
                     '#SBATCH --ntasks-per-node=1',
                     '#SBATCH --mem-per-cpu=' + mem,
                     '\nmodule load matlab',
                     '\nmatlab -nodisplay < '+'/home/livingstonb/GitHub/MPCrecode/batch/mfiles/master'
                         +label+'_${SLURM_ARRAY_TASK_ID}.m']
            newmfile.write('\n'.join(lines))
            mfile.seek(0)
                    
def gen_mfile_aggregator(args):
    # create m-file to aggregate .mat files and create a table

    matdir = args['MWout']
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
    
    with open(os.path.join(args['MWout'],'aggregate.m'),'w') as newmfile:
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

args = {}
# location of master.m
args['masterpath'] = '/Users/Brian/Documents/GitHub/MPCrecode/master.m'
# path of 'batch' folder
args['batchpath'] = '/Users/Brian/Documents/GitHub/MPCrecode/batch'
# relative path of income variable inside MPCrecode
args['Qincvar'] = "'IncomeGrids/quarterly_b.mat'"
# location of .mat output files
args['MWout'] = '/Users/Brian/Documents/midway2temp'
# run with small grids for speed (string)
args['fast'] = '0'
# frequencies to run (list of strings)
args['frequencies'] = ['1','4']
# time allocation for each mfile
args['time'] = '01:00:00'

# names for each nb
names = {};
names['1'] = ['{}'] # make matlab use all names
names['5'] = [] # each instance of matlab will use one name
for d in ['NoDeath','Death']:
    for ibw in list(map(str,[0.001,0.005,0.01])):
        names['5'].append("{{'2 FixedBetaHet5 Width{0:s} {1:s}'}}".format(ibw,d))
        for bs in list(map(str,[0.02,0.1])):
            names['5'].append("{{'2 RandomBetaHet5 Width{0:s} SwitchProb{1:s} {2:s}'}}".format(ibw,bs,d))
args['names'] = names

# open file and create new directories if necessary
mfile = open(args['masterpath'])
if not os.path.exists(os.path.join(args['batchpath'],'sbatch')):
    os.mkdir(os.path.join(args['batchpath'],'sbatch'))
if not os.path.exists(os.path.join(args['batchpath'],'mfiles')):
    os.mkdir(os.path.join(args['batchpath'],'mfiles'))

# ---------------------------------------------------------------------
# Function calls
# ---------------------------------------------------------------------
# create batch m-files
#(nb1end,nb5end) = gen_mfiles(mfile,args)

# create batch sbatch scripts
#gen_sbatch(mfile,args,nb1end,nb5end)

# generate m-file that aggreagates .mat files and creates a table
gen_mfile_aggregator(args)
    
mfile.close()
