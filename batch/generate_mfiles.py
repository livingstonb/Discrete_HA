
import os

def gen_mfiles(mfile,incvars,filebase,suffix,frequencies):
    deletepath = os.path.join(filebase,'mfiles')
    for dirpath, dirs, files in os.walk(deletepath):
        for nfile in files:
            os.remove(os.path.join(dirpath,nfile))
    for i in range(len(incvars)):
        count = 0
        for nb in [1,5]:
            for name in names[nb]:
                count += 1
                suffix = incsuffix[i]+'_spec'+str(count)
                with open(filebase+'mfiles/master'+suffix+'.m','w') as newmfile:
                    for line in mfile.readlines():
                        if line.startswith('IncomeProcess'):
                            newmfile.write("IncomeProcess = "+incvars[i]+';\n')
                        elif line.startswith('selection.names_to_run'):
                            newmfile.write('selection.names_to_run = '+name+';\n')
                        elif line.startswith('selection.suffix'):
                            newmfile.write("selection.suffix = '"+suffix+"';\n")
                        elif line.startswith('selection.frequencies'):
                            newmfile.write('selection.frequencies = {};\n'.format(frequencies[i]))
                        elif line.startswith('selection.nb'):
                            newmfile.write('selection.nb = {:d};\n'.format(nb))
                        else:
                            newmfile.write(line)
                    mfile.seek(0)
            
def gen_sbatch(mfile,incvars,filebase,incsuffix,frequencies,sbatch):
    for i in range(len(incvars)):
        count = 0
        for nb in [1,5]
            count += 1
            suffix = incsuffix[i]+'_spec'+str(count)
            with open(filebase+'sbatch/batch'+suffix+'.sbatch','w') as newmfile:
                lines = ['#!/bin/bash',
                         '#SBATCH --job-name=' + suffix,
                         '#SBATCH --output=' + sbatch['output'][i],
                         '#SBATCH --error=' + sbatch['error'][i],
                         '#SBATCH --partition=broadwl',
                         '#SBATCH --time=04:00:00',
                         '#SBATCH --nodes=1',
                         '#SBATCH --ntasks-per-node=1',
                         '#SBATCH --mem-per-cpu=4000',
                         '\n module load matlab',
                         '\n matlab -nodisplay < ' + sbatch['mfile'][i]]
                newmfile.write('\n'.join(lines))
                mfile.seek(0)
            

# ---------------------------------------------------------------------
# Startup
# ---------------------------------------------------------------------
mfile = open('/Users/Brian/Documents/GitHub/MPCrecode/master.m')

incvars = ["''",
            "''",
            "'IncomeVariables/quarterly_a.mat'",
            "'IncomeVariables/quarterly_b.mat'",
            "'IncomeVariables/quarterly_c.mat'"]
            
filebase = '/Users/Brian/Documents/GitHub/MPCrecode/batch/'
if not os.path.exists(filebase+'sbatch'):
    os.mkdir(filebase+'sbatch')
if not os.path.exists(filebase+'mfiles'):
    os.mkdir(filebase+'mfiles')
    
incsuffix = ['_A','_Q','_Qa','_Qb','_Qc']
frequencies = [1,4,4,4,4]
# names for each nb
names = {};
names[1] = ['{}']
names[5] = []
for d in ['NoDeath','Death']:
    for ibw in list(map(str,[0.001,0.005,0.01])):
        names[5].append("{{'2 FixedBetaHet5 Width{0:s} {1:s}'}}".format(ibw,d))
        for bs in list(map(str,[0.02,0.1])):
            names[5].append("{{'2 RandomBetaHet5 Width{0:s} SwitchProb{1:s} {2:s}'}}".format(ibw,bs,d))

# sbatch variables
sbatch = {'output':[],'error':[],'mfile':[]}
for i in range(len(incvars)):
    sbatch['output'].append('/home/livingstonb/output/matlab' + incsuffix[i] + '.out')
    sbatch['error'].append('/home/livingstonb/output/matlab' + incsuffix[i] + '.err')
    sbatch['mfile'].append('/home/livingstonb/GitHub/MPCrecode/batch/mfiles/master' + incsuffix[i] + '.m')

# ---------------------------------------------------------------------
# Function calls
# ---------------------------------------------------------------------
gen_mfiles(mfile,incvars,filebase,incsuffix,frequencies)
#gen_sbatch(mfile,incvars,filebase,suffix,frequencies,sbatch)
    
mfile.close()
