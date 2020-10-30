
serverpath="$MW:/home/livingstonb/GitHub/Discrete_HA/output/*.csv"
outpath="/home/brian/Documents/GitHub/Discrete_HA/output/server/"
scp $serverpath $outpath


spath := "$$MW:/home/livingstonb/GitHub/Discrete_HA/output/*.csv"
download :
	-mkdir -p output/server
	-scp $(spath) output/server/