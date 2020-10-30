batch :
	sbatch "code/batch/server.sbatch"

combine :
	sbatch "code/batch/combine_runs.sbatch"

spath := "$$MW:/home/livingstonb/GitHub/Discrete_HA/output/*.xlsx"
cdate := $(shell date +"%m-%d-%Y-%T")
download :
	-mkdir -p output/server-$(cdate)
	-scp $(spath) output/server-$(cdate)

.PHONY : batch combine download