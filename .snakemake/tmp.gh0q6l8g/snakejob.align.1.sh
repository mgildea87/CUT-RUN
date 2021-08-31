#!/bin/sh
# properties = {"type": "single", "rule": "align", "local": false, "input": ["trim/Sample_1_ExpControl_1_trimmed_R1.fastq.gz", "trim/Sample_1_ExpControl_1_trimmed_R2.fastq.gz"], "output": ["alignment/Sample_1_ExpControl_1.bam", "alignment/Sample_1_ExpControl_1_ecoli.bam"], "wildcards": {"sample": "Sample_1_ExpControl_1"}, "params": {"p1": "--end-to-end --very-sensitive --no-mixed --no-unal --no-discordant --phred33 -I 10 -X 700", "p2": "--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-unal --no-discordant --phred33 -I 10 -X 700"}, "log": ["logs/alignment_reports/Sample_1_ExpControl_1.log"], "threads": 20, "resources": {}, "jobid": 1, "cluster": {"Job_name": "align", "Mem": 60000, "Partition": "cpu_short", "Cores": 20, "Time": 240, "Error": "logs/slurm_reports/slurm-%j.out"}}
 cd /gpfs/home/gildem01/workflows/CUT-RUN && \
PATH='/gpfs/data/fisherlab/conda_envs/CUT-RUN/bin':$PATH /gpfs/data/fisherlab/conda_envs/CUT-RUN/bin/python3.9 \
-m snakemake alignment/Sample_1_ExpControl_1.bam --snakefile /gpfs/home/gildem01/workflows/CUT-RUN/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.gh0q6l8g trim/Sample_1_ExpControl_1_trimmed_R1.fastq.gz trim/Sample_1_ExpControl_1_trimmed_R2.fastq.gz --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules align --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.gh0q6l8g/1.jobfinished || (touch /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.gh0q6l8g/1.jobfailed; exit 1)

