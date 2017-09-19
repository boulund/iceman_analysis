# Snakefile for Iceman analysis
# (c) Fredrik Boulund 2017
# Edit config.yaml to conform to your system, 
# then run with snakemake flags:
# 	--snakefile path/to/Snakefile 
# 	--configfile path/to/config.yaml
#   --cores <to_your_liking>
#



rule all:
	input:
		expand('qc_reads/{sample}_R{readnum}.fq.gz', sample=config["samples"], readnum=["1", "2"]),
		expand('qc_reads/{background_sample}_R{readnum}.fq.gz', background_sample=config["background_sample"], readnum=["1", "2"]),
		expand('interleaved_background_reads/{background_sample}.interleaved.fq.gz', background_sample=config["background_sample"]),
		expand('filtered_human/{sample}_R{readnum}.human_filtered.fq.gz', sample=config["samples"], readnum=["1", "2"]),
		expand('filtered_background/{sample}_R{readnum}.background_filtered.fq.gz', sample=config["samples"], readnum=["1", "2"]),
		expand('merged_subtracted_reads/{sample}.merged_subtracted.fq.gz', sample=config["samples"]),
		expand('mapDamage/{sample}/Fragmisincorporation_plot.pdf', sample=config["samples"]),


rule qc_reads:
	"""Run reads through BBDuk for QC."""
	input:
		in1 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R1.fastq.gz',
		in2 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R2.fastq.gz',
	output:
		out1 = 'qc_reads/{sample}_R1.fq.gz',
		out2 = 'qc_reads/{sample}_R2.fq.gz',
		stats = 'qc_reads/{sample}/{sample}.stats.txt',
		bhist = 'qc_reads/{sample}/{sample}.bhist.txt',
		qhist = 'qc_reads/{sample}/{sample}.qhist.txt',
		qchist = 'qc_reads/{sample}/{sample}.qchist.txt',
		aqhist = 'qc_reads/{sample}/{sample}.aqhist.txt',
		bqhist = 'qc_reads/{sample}/{sample}.bqhist.txt',
		lhist = 'qc_reads/{sample}/{sample}.lhist.txt',
		gchist = 'qc_reads/{sample}/{sample}.gchist.txt',
	threads: 40
	shell:
		"""
		bbduk.sh \
			in1={input.in1} \
			in2={input.in2} \
			ref={config[bbduk_ref]} \
			out1={output.out1} \
			out2={output.out2} \
			stats={output.stats} \
			bhist={output.bhist} \
			qhist={output.qhist} \
			qchist={output.qchist} \
			aqhist={output.aqhist} \
			bqhist={output.bqhist} \
			lhist={output.lhist} \
			gchist={output.gchist} \
			qtrim=rl \
			trimq=20 \
			ktrim=r \
			k=25 \
			mink=11 \
			hdist=1 
		"""


rule interleave_background_reads:
	"""Interleave reads from two paired end fastq into a single file."""
	input:
		in1 = 'qc_reads/{sample}_R1.fq.gz',
		in2 = 'qc_reads/{sample}_R2.fq.gz',
	output:
		'interleaved_background_reads/{sample}.interleaved.fq.gz'
	threads: 40
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
		read1 = 'qc_reads/{sample}_R1.fq.gz',
		read2 = 'qc_reads/{sample}_R2.fq.gz',
	output:
 		out_read1 = 'filtered_human/{sample}_R1.human_filtered.fq.gz',
		out_read2 = 'filtered_human/{sample}_R2.human_filtered.fq.gz',
		out_matched = 'filtered_human/{sample}_matched.human_filtered.sam.gz',
	threads: 40
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
	threads: 40
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
		'merged_subtracted_reads/{sample}.merged_subtracted.fq.gz'
	threads: 40
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
		
		
rule metaphlan2:
	"""Run MetaPhlAn2 on the merged subtracted reads."""
	input:
		reads = 'merged_subtracted_reads/{sample}.merged_subtracted.fq.gz'
	output:
		
	shell:
		"""
		"""
