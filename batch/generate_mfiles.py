
import os

def gen_mfiles(mfile,incvars,filebase,suffix,frequencies):
    for i in range(len(incvars)):
        with open(filebase+'mfiles/master'+suffix[i]+'.m','w') as newmfile:
            for line in mfile.readlines():
                if line.startswith('IncomeProcess'):
                    newmfile.write("IncomeProcess = "+incvars[i]+';\n')
                elif line.startswith('selection.suffix'):
                    newmfile.write("selection.suffix = '"+suffix[i]+"';\n")
                elif line.startswith('selection.frequencies'):
                    newmfile.write('selection.frequencies = {};\n'.format(frequencies[i]))
                else:
                    newmfile.write(line)
            mfile.seek(0)
            
def gen_sbatch(mfile,incvars,filebase,suffix,frequencies,sbatch):
    for i in range(len(incvars)):
        with open(filebase+'sbatch/batch'+suffix[i]+'.sbatch','w') as newmfile:
            lines = ['#!/bin/bash',
                     '#SBATCH --job-name=' + sbatch['jname'][i],
                     '#SBATCH --output=' + sbatch['output'][i],
                     '#SBATCH --error=' + sbatch['error'][i],
                     '#SBATCH --partition=broadwl',
                     '#SBATCH --time=02:00:00',
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
suffix = ['_A','_Q','_Qa','_Qb','_Qc']
frequencies = [1,4,4,4,4]

# sbatch variables
sbatch = {'jname':[],'output':[],'error':[],'mfile':[]}
for i in range(len(incvars)):
    sbatch['jname'].append('MPCs' + suffix[i])
    sbatch['output'].append('/home/livingstonb/output/matlab' + suffix[i] + '.out')
    sbatch['error'].append('/home/livingstonb/output/matlab' + suffix[i] + '.err')
    sbatch['mfile'].append('/home/livingstonb/GitHub/MPCrecode/batch/mfiles/master' + suffix[i] + '.m')

# ---------------------------------------------------------------------
# Function calls
# ---------------------------------------------------------------------
gen_mfiles(mfile,incvars,filebase,suffix,frequencies)
gen_sbatch(mfile,incvars,filebase,suffix,frequencies,sbatch)
    
mfile.close()
