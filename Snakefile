import pandas as pd
import os

for directory in ['fastqc', 'fastqc_post_trim', 'trim', 'logs', 'logs/slurm_reports', 'logs/trim_reports', 'alignment','alignment/bed', 'alignment/frag_len', 'logs/alignment_reports', 'peaks']:
	if not os.path.isdir(directory):
		os.mkdir(directory)

configfile: "config.yaml"
sample_file = config["sample_file"]
sample = pd.read_table(sample_file)['Sample']
replicate = pd.read_table(sample_file)['Replicate']
condition = pd.read_table(sample_file)['Condition']
Antibody = pd.read_table(sample_file)['Antibody']
File_R1 = pd.read_table(sample_file)['File_Name_R1']
File_R2 = pd.read_table(sample_file)['File_Name_R2']
File_names = File_R1.append(File_R2)
genome = config["genome"]
spike_genome = config["spike_genome"]

sample_ids = []
for i in range(len(sample)):
	sample_ids.append('%s_%s_%s' % (sample[i], condition[i], replicate[i]))
sample_ids = pd.unique(sample_ids).tolist()

sample_ids_file = []
for i in range(len(sample)):
	sample_ids_file.append('%s_%s_%s_%s' % (sample[i], condition[i], replicate[i], Antibody[i]))

read = ['_R1', '_R2']

rule all:
	input:
		expand('fastqc/{sample_file}{read}_fastqc.html', sample_file = sample_ids_file, read = read),
		expand('fastqc_post_trim/{sample_file}_trimmed{read}_fastqc.html', sample_file = sample_ids_file, read = read),
		expand('peaks/{sample}.stringent.bed', sample = sample_ids),
		'FRP.txt',
		expand('alignment/frag_len/{sample}.txt', sample = sample_ids_file)

rule fastqc:
	input: 
		fastq = "fastq/{sample}{read}.fastq.gz"
	output:  
		"fastqc/{sample}{read}_fastqc.html",
	params:
		'fastqc/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule fastqc_post_trim:
	input: 
		fastq = "trim/{sample}{read}.fastq.gz"
	output:  
		"fastqc_post_trim/{sample}{read}_fastqc.html",
	params:
		'fastqc_post_trim/'
	shell: 
		'fastqc {input.fastq} -o {params}'

rule trim:
	input:
		R1='fastq/{sample}_R1.fastq.gz',
		R2='fastq/{sample}_R2.fastq.gz'
	output:
		R1='trim/{sample}_trimmed_R1.fastq.gz',
		R2='trim/{sample}_trimmed_R2.fastq.gz',
		html='logs/trim_reports/{sample}.html',
		json='logs/trim_reports/{sample}.json'
	threads: 4
	log:
		'logs/trim_reports/{sample}.log'
	params:
		'--detect_adapter_for_pe'
	shell:
		'fastp -w {threads} {params} -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --html {output.html} --json {output.json} 2> {log}'

rule align:
	input:
		R1='trim/{sample}_trimmed_R1.fastq.gz',
		R2='trim/{sample}_trimmed_R2.fastq.gz'
	output:
		'alignment/{sample}.bam'
	threads: 20
	log:
		'logs/alignment_reports/{sample}.log'
	params:
		'--end-to-end --very-sensitive --no-mixed --no-unal --no-discordant --phred33 -I 10 -X 700'
	shell:
		'bowtie2 {params} -x %s --threads {threads} -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -bh -q 3 > alignment/{wildcards.sample}.bam' % (genome)

rule align_spike:
	input:
		R1='trim/{sample}_trimmed_R1.fastq.gz',
		R2='trim/{sample}_trimmed_R2.fastq.gz'
	output:
		'alignment/{sample}_ecoli.bam'
	threads: 20
	log:
		'logs/alignment_reports/{sample}_ecoli.log'
	params:
		'--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-unal --no-discordant --phred33 -I 10 -X 700'
	shell:
		'bowtie2 {params} -x %s --threads {threads} -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -bh -q 3 > alignment/{wildcards.sample}_ecoli.bam' % (spike_genome)

rule spike_in_norm:
	input:
		sample_bam='alignment/{sample}.bam',
		spike_bam='alignment/{sample}_ecoli.bam'
	output:
		'alignment/bed/{sample}.bed',
		'alignment/bed/{sample}.bedgraph'
	shell:
		"""
		depth=`samtools view alignment/{wildcards.sample}_ecoli.bam | wc -l`
		depth=$((depth/2))
		echo $depth
		scale_fac=`echo "10000 / $depth" | bc -l`
		echo $scale_fac
		bedtools bamtobed -bedpe -i alignment/{wildcards.sample}.bam | cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n > alignment/bed/{wildcards.sample}.bed
		bedtools genomecov -bg -i alignment/bed/{wildcards.sample}.bed -scale $scale_fac -g /gpfs/data/fisherlab/genomes/mm10/STAR_75_2.7.7a/chrNameLength.txt > alignment/bed/{wildcards.sample}.bedgraph
		"""

rule SEACR:
	input:
		exp='alignment/bed/{sample}_Antibody.bedgraph',
		con='alignment/bed/{sample}_Control.bedgraph'
	output:
		'peaks/{sample}.stringent.bed'
	params:
		'non stringent'
	shell:
		'bash SEACR_1.3.sh {input.exp} {input.con} {params} peaks/{wildcards.sample}'

rule fragment_size:
	input:
		'alignment/{sample}.bam'
	output:
		'alignment/frag_len/{sample}.txt'
	shell:
		"""
		samtools view {input} | awk -F'\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{print abs($9)}}' | sort | uniq -c | awk -v OFS="\t" '{{print $2, $1/2}}' > {output}
		"""

rule FRP:
	input:
		expand('peaks/{sample}.stringent.bed', sample = sample_ids)
	output:
		'FRP.txt'
	script:
		'FRP.py'