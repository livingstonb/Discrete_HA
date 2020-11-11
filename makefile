batch :
	sbatch "code/batch/server.sbatch"

combine :
	sbatch "code/batch/combine_runs.sbatch"

spath := "$$MW:/home/livingstonb/GitHub/Discrete_HA/output/tables/*"
cdate := $(shell date +"%m-%d-%Y-%T")
download : tables
	-mkdir -p output/server-$(cdate)
	-scp $(spath) output/server-$(cdate)

tables :
	python code/+tables/create_tex_tables.py "output/tables"
	cp code/+tables/tables.tex output/tables/tables.tex
	cd output/tables && pdflatex tables
	cd output/tables && pdflatex tables

clean_tables :
	rm output/tables/*

.PHONY : batch combine download tables clean_tables