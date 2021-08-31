#!/bin/sh
# properties = {"type": "single", "rule": "SEACR", "local": false, "input": ["alignment/bed/Sample_2_IgGControl_1.bedgraph", "alignment/bed/Sample_2_IgGControl_1.bedgraph"], "output": ["peaks/Sample_2_IgGControl_1.peaks"], "wildcards": {"sample": "Sample_2_IgGControl_1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 10, "cluster": {"Job_name": "SEACR", "Mem": 60000, "Partition": "cpu_short", "Cores": 20, "Time": 240, "Error": "logs/slurm_reports/slurm-%j.out"}}
 cd /gpfs/home/gildem01/workflows/CUT-RUN && \
PATH='/gpfs/data/fisherlab/conda_envs/CUT-RUN/bin':$PATH /gpfs/data/fisherlab/conda_envs/CUT-RUN/bin/python3.9 \
-m snakemake peaks/Sample_2_IgGControl_1.peaks --snakefile /gpfs/home/gildem01/workflows/CUT-RUN/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.s7a7l_ff alignment/bed/Sample_2_IgGControl_1.bedgraph alignment/bed/Sample_2_IgGControl_1.bedgraph --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules SEACR --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.s7a7l_ff/10.jobfinished || (touch /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.s7a7l_ff/10.jobfailed; exit 1)

