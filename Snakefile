# Snakefile for Iceman analysis
# (c) Fredrik Boulund 2017



rule all:
	input:
		expand('interleaved_background_reads/1018.interleaved.fq.gz'),
		expand('filtered_human/{sample}_R{readnum}.human_filtered.fq.gz', sample=config["samples"], readnum=["1", "2"]),
		expand('filtered_background/{sample}_R{readnum}.background_filtered.fq.gz', sample=config["samples"], readnum=["1", "2"]),
		expand('merged_reads/{sample}.merged.fq.gz', sample=config["samples"]),
		expand('mapDamage/{sample}/Fragmisincorporation_plot.pdf', sample=config["samples"]),


rule interleave_background_reads:
	"""Interleave reads from two paired end fastq into a single file."""
	input:
		in1 = config['background_reads_dir']+'/1018_Iceman_R1.fastq.gz',
		in2 = config['background_reads_dir']+'/1018_Iceman_R2.fastq.gz',
	output:
		'interleaved_background_reads/1018.interleaved.fq.gz'
	threads: 20
	shell:
		"""
		reformat.sh \
			in1={input.in1} \
			in2={input.in2} \
			out={output}
			threads={threads}
		"""


rule filter_human:
	"""Filter (remove) human contamination."""
	input:
		read1 = config['samples_dir']+'/{sample}_Iceman_R1.fastq.gz',
		read2 = config['samples_dir']+'/{sample}_Iceman_R2.fastq.gz',
	output:
 		out_read1 = 'filtered_human/{sample}_R1.human_filtered.fq.gz',
		out_read2 = 'filtered_human/{sample}_R2.human_filtered.fq.gz',
		out_matched = 'filtered_human/{sample}_matched.human_filtered.sam.gz',
	threads: 20
	shell:
		"""
		bbmap.sh \
			in1={input.read1} \
			in2={input.read2} \
			out1={output.out_read1} \
			out2={output.out_read2} \
			outm={output.out_matched} \
			path={config[bbduk_human_db_path]} \
			statsfile=logs/filtered_human/{wildcards.sample}.statsfile.txt \
			minid=0.95 \
			maxindel=3 \
			minhits=2 \
			bandwidthratio=0.16 \
			bandwidth=12 \
			quickmatch \
			fast \
			qtrim=rl \
			trimq=10 \
			untrim \
			threads={threads}
		"""


rule filter_background:
	"""Filter (remove) background reads."""
	input:
		read1 = 'filtered_human/{sample}_R1.human_filtered.fq.gz',
		read2 = 'filtered_human/{sample}_R2.human_filtered.fq.gz',
		interleaved_lung_reads = rules.interleave_background_reads.output,
	output:
		out_read1 = 'filtered_background/{sample}_R1.background_filtered.fq.gz',
		out_read2 = 'filtered_background/{sample}_R2.background_filtered.fq.gz',
	threads: 20
	shell:
		"""
		bbmap.sh \
			in1={input.read1} \
			in2={input.read2} \
			out1={output.out_read1} \
			out2={output.out_read2} \
			ref={input.interleaved_lung_reads} \
			statsfile=logs/filtered_background/{wildcards.sample}.statsfile.txt \
			minid=0.95 \
			maxindel=3 \
			minhits=2 \
			bandwidthratio=0.16 \
			bandwidth=12 \
			quickmatch \
			fast \
			qtrim=rl \
			trimq=10 \
			untrim \
			threads={threads}
		"""


rule merge_reads:
	"""Merge overlapping paired end reads."""
	input: 
		read1 = 'filtered_background/{sample}_R1.background_filtered.fq.gz',
		read2 = 'filtered_background/{sample}_R2.background_filtered.fq.gz'
	output:
		'merged_reads/{sample}.merged.fq.gz'
	threads: 20
	shell:
		"""
		bbmerge.sh \
			in1={input.read1} \
			in2={input.read2} \
			out={output}
			threads={threads}
		"""


rule mapDamage:
	"""Run mapDamage to assess ancient DNA content."""
	input:
		rules.filter_human.output.out_matched
	output:
		'mapDamage/{sample}/Fragmisincorporation_plot.pdf'
	shell:
		"""
		source activate mapDamage
		mapDamage \
			--input {input} \
			--reference {config[hg19_fasta]} \
			--folder mapDamage/{wildcards.sample} \
		"""
		
		
		
