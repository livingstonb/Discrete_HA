
import os
import re

def gen_mfiles(mfile,args):
    # delete any existing batch m-files
    deletepath = os.path.join(args['batchpath'],'mfiles')
    for dirpath, dirs, files in os.walk(deletepath):
        for nfile in files:
            os.remove(os.path.join(dirpath,nfile))
            
    # loop over specifications, identified by '_spec<#>' in files
    count = 0
    for freq in ['1','4']:
        for nb in ['1','5']:
            for name in args['names'][nb]:
                count += 1
                suffix = '_spec'+str(count)
                with open(os.path.join(args['batchpath'],'mfiles/master')+suffix+'.m','w') as newmfile:
                    for line in mfile.readlines():
                        if line.startswith('runopts.Batch ='):
                            newmfile.write("runopts.Batch = 1;\n")
                        elif line.startswith('runopts.Server ='):
                            newmfile.write('runopts.Server = 1;\n')
                        elif line.startswith('IncomeProcess'):
                            newmfile.write('IncomeProcess = '+args['Qincvar']+';\n')
                        elif line.startswith('selection.names_to_run'):
                            newmfile.write('selection.names_to_run = '+name+';\n')
                        elif line.startswith('selection.suffix'):
                            newmfile.write("selection.suffix = '"+suffix+"';\n")
                        elif line.startswith('selection.frequencies'):
                            newmfile.write('selection.frequencies = '+freq+';\n')
                        elif line.startswith('selection.nb'):
                            newmfile.write('selection.nb = '+nb+';\n')
                        else:
                            newmfile.write(line)
                    mfile.seek(0)
            
def gen_sbatch(mfile,args):
    deletepath = os.path.join(args['batchpath'],'sbatch')
    for dirpath, dirs, files in os.walk(deletepath):
        for nfile in files:
            os.remove(os.path.join(dirpath,nfile))
            
    # loop over specifications, identified by '_spec<#>' in files
    count = 0
    for freq in ['1','4']:
        for nb in ['1','5']:
            for name in args['names'][nb]:
                count += 1
                if nb == '5':
                    time = '01:00:00'
                else:
                    time = '00:25:00'
                    
                if (freq=='4') and (nb=='5'):
                    mem = '4000'
                else:
                    mem = '2000'
                
                suffix = 'spec'+str(count)
                with open(os.path.join(args['batchpath'],'sbatch/')+suffix+'.sbatch','w') as newmfile:
                    lines = ['#!/bin/bash',
                             '#SBATCH --job-name=disc' + '_' +suffix,
                             '#SBATCH --output='+'/home/livingstonb/output/matlab'+'_'+suffix+'.out',
                             '#SBATCH --error=' +'/home/livingstonb/output/matlab'+'_'+suffix+'.err',
                             '#SBATCH --partition=broadwl',
                             '#SBATCH --' + time,
                             '#SBATCH --nodes=1',
                             '#SBATCH --ntasks-per-node=1',
                             '#SBATCH --mem-per-cpu=' + mem,
                             '\nmodule load matlab',
                             '\nmatlab -nodisplay < '+'/home/livingstonb/GitHub/MPCrecode/batch/mfiles/master'+'_'+suffix+'.m']
                    newmfile.write('\n'.join(lines))
                    mfile.seek(0)
                    
def gen_mfile_aggregator(args):
    # create m-file to aggregate .mat files
    matdir = args['MWout']
    matfiles = os.listdir(matdir)
    pattern = re.compile('spec[0-9]+.mat')
    matfiles = list(filter(lambda x: re.fullmatch(pattern,x),matfiles))
    
    matfiles = sorted(matfiles)
    matfiles = list(map(lambda x: os.path.join(matdir,x),matfiles))
    
    with open(os.path.join(args['MWout'],'aggregate.m'),'w') as newmfile:
        newmfile.write("clear\n% Aggregates spec##.mat's into one .mat file\n\n")
        for i,matfile in enumerate(matfiles):
            newmfile.write("matfiles{{{0:d}}} = '{1:s}';\n".format(i+1,matfile))
        
        lines = ['\nspec = 0;',
                 'for im = 1:numel(matfiles)',
                 '    S = load(matfiles{im});',
                 '    for is = 1:numel(S)',
                 '        spec = spec + 1;',
                 '        out(spec) = S(is);',
                 '    end',
                 'end']
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
args['Qincvar'] = "'IncomeVariables/quarterly_b.mat'"
# location of .mat output files
args['MWout'] = '/Users/Brian/Documents/midway2temp'

mfile = open(args['masterpath'])

if not os.path.exists(os.path.join(args['batchpath'],'sbatch')):
    os.mkdir(os.path.join(args['batchpath'],'sbatch'))
if not os.path.exists(os.path.join(args['batchpath'],'mfiles')):
    os.mkdir(os.path.join(args['batchpath'],'mfiles'))

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

# ---------------------------------------------------------------------
# Function calls
# ---------------------------------------------------------------------
gen_mfiles(mfile,args)
gen_sbatch(mfile,args)
gen_mfile_aggregator(args)
    
mfile.close()
