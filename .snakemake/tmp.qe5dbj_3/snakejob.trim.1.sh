#!/bin/sh
# properties = {"type": "single", "rule": "trim", "local": false, "input": ["fastq/Sample_1_ExpControl_1_R1.fastq.gz", "fastq/Sample_1_ExpControl_1_R2.fastq.gz"], "output": ["trim/Sample_1_ExpControl_1_trimmed_R1.fastq.gz", "trim/Sample_1_ExpControl_1_trimmed_R2.fastq.gz", "logs/trim_reports/Sample_1_ExpControl_1.html", "logs/trim_reports/Sample_1_ExpControl_1.json"], "wildcards": {"sample": "Sample_1_ExpControl_1"}, "params": {}, "log": ["logs/trim_reports/Sample_1_ExpControl_1.log"], "threads": 4, "resources": {}, "jobid": 1, "cluster": {}}
 cd /gpfs/home/gildem01/workflows/CUT-RUN && \
PATH='/gpfs/data/fisherlab/conda_envs/CUT-RUN/bin':$PATH /gpfs/data/fisherlab/conda_envs/CUT-RUN/bin/python3.9 \
-m snakemake trim/Sample_1_ExpControl_1_trimmed_R1.fastq.gz --snakefile /gpfs/home/gildem01/workflows/CUT-RUN/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.qe5dbj_3 fastq/Sample_1_ExpControl_1_R1.fastq.gz fastq/Sample_1_ExpControl_1_R2.fastq.gz --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules trim --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.qe5dbj_3/1.jobfinished || (touch /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.qe5dbj_3/1.jobfailed; exit 1)

