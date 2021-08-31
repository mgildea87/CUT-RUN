#!/bin/sh
# properties = {"type": "single", "rule": "spike_in_norm", "local": false, "input": ["alignment/Sample_1_ExpControl_1.bam", "alignment/Sample_1_ExpControl_1_ecoli.bam"], "output": ["alignment/bed/Sample_1_ExpControl_1.bed", "alignment/bed/Sample_1_ExpControl_1.bedgraph"], "wildcards": {"sample": "Sample_1_ExpControl_1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 1, "cluster": {}}
 cd /gpfs/home/gildem01/workflows/CUT-RUN && \
PATH='/gpfs/data/fisherlab/conda_envs/CUT-RUN/bin':$PATH /gpfs/data/fisherlab/conda_envs/CUT-RUN/bin/python3.9 \
-m snakemake alignment/bed/Sample_1_ExpControl_1.bed --snakefile /gpfs/home/gildem01/workflows/CUT-RUN/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.ghwq2kcq alignment/Sample_1_ExpControl_1.bam alignment/Sample_1_ExpControl_1_ecoli.bam --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules spike_in_norm --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.ghwq2kcq/1.jobfinished || (touch /gpfs/home/gildem01/workflows/CUT-RUN/.snakemake/tmp.ghwq2kcq/1.jobfailed; exit 1)

