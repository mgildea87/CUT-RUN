#!/bin/bash -l
module load miniconda3/cpu/4.9.2
conda activate /gpfs/data/fisherlab/conda_envs/CUT-RUN
mkdir fastq
if ! python cat_rename.py /gpfs/data/sequence/results/moorelab/2021-07-27/fastq/ ${1}; then
    echo "Exiting..."
    exit
 fi

snakemake --cluster "sbatch -J {cluster.Job_name} --mem {cluster.Mem} -c {cluster.Cores} -p {cluster.Partition} -t {cluster.Time} --output {cluster.Error}" --cluster-config cluster_config.yml --jobs 6
snakemake --report snake_make_report.html
multiqc .