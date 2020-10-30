batch :
	sbatch "code/batch/server.sbatch"

combine :
	sbatch "code/batch/combine_runs.sbatch"

spath := "$$MW:/home/livingstonb/GitHub/Discrete_HA/output/*.csv"
download :
	-mkdir -p output/server
	-scp $(spath) output/server/

.PHONY : batch combine download