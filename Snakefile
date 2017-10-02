# Snakefile for Iceman analysis
# vim: syntax=python expandtab
# (c) Fredrik Boulund 2017
# Edit config.yaml to conform to your system, 
# then run with snakemake flags:
#     --snakefile path/to/Snakefile 
#     --configfile path/to/config.yaml
#     --cores <to_your_liking>
#



rule all:
    input:
        expand('qc_reads/{sample}_R{readnum}.fq.gz', sample=config["samples"], readnum=["1", "2"]),
        expand('qc_reads/{background_sample}_R{readnum}.fq.gz', background_sample=config["background_sample"], readnum=["1", "2"]),
        expand('merged_background_reads/{background_sample}.merged.fa.gz', background_sample=config["background_sample"]),
        expand('assembled_background_reads/{background_sample}.assembled.fa.gz', background_sample=config["background_sample"]),
        expand('filtered_human/{sample}_R{readnum}.human_filtered.fq.gz', sample=config["samples"], readnum=["1", "2"]),
        expand('subtracted_background/{sample}_R{readnum}.background_subtracted.fq.gz', sample=config["samples"], readnum=["1", "2"]),
        expand('merged_subtracted_reads/{sample}.merged_subtracted.fq.gz', sample=config["samples"]),
        expand('mapDamage/{sample}/Fragmisincorporation_plot.pdf', sample=config["samples"]),
        expand('metaphlan2/metaphlan2_outputs/{sample}.metaphlan2.txt', sample=config["samples"]),
        expand('metaphlan2/metaphlan2_outputs/pre/{sample}.metaphlan2.txt', sample=config["samples"]),
        expand('metaphlan2/metaphlan2_outputs/pre/{background_sample}.metaphlan2.txt', background_sample=config["background_sample"]),


rule assess_input_data_quality:
    """Run reads through BBDuk for quality stats."""
    input:
        in1 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R1.fastq.gz',
        in2 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R2.fastq.gz',
    output:
        stats =  'logs/input_read_quality/{sample}.stats.txt',
        bhist =  'logs/input_read_quality/{sample}.bhist.txt',
        qhist =  'logs/input_read_quality/{sample}.qhist.txt',
        qchist = 'logs/input_read_quality/{sample}.qchist.txt',
        aqhist = 'logs/input_read_quality/{sample}.aqhist.txt',
        bqhist = 'logs/input_read_quality/{sample}.bqhist.txt',
        lhist =  'logs/input_read_quality/{sample}.lhist.txt',
        gchist = 'logs/input_read_quality/{sample}.gchist.txt',
    threads: 40
    shell:
        """
        bbduk.sh \
            in1={input.in1} \
            in2={input.in2} \
            stats={output.stats} \
            bhist={output.bhist} \
            qhist={output.qhist} \
            qchist={output.qchist} \
            aqhist={output.aqhist} \
            bqhist={output.bqhist} \
            lhist={output.lhist} \
            gchist={output.gchist} \
        """


rule qc_reads:
    """Run reads through BBDuk for QC and adapter trimming."""
    input:
        in1 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R1.fastq.gz',
        in2 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R2.fastq.gz',
    output:
        out1 = 'qc_reads/{sample}_R1.fq.gz',
        out2 = 'qc_reads/{sample}_R2.fq.gz',
        stats =  'logs/qc_reads/{sample}.stats.txt',
        bhist =  'logs/qc_reads/{sample}.bhist.txt',
        qhist =  'logs/qc_reads/{sample}.qhist.txt',
        qchist = 'logs/qc_reads/{sample}.qchist.txt',
        aqhist = 'logs/qc_reads/{sample}.aqhist.txt',
        bqhist = 'logs/qc_reads/{sample}.bqhist.txt',
        lhist =  'logs/qc_reads/{sample}.lhist.txt',
        gchist = 'logs/qc_reads/{sample}.gchist.txt',
    threads: 40
    shell:
        """
        bbduk.sh \
            in1={input.in1} \
            in2={input.in2} \
            ref={config[bbduk_ref]} \
            build=1 \
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


rule merge_background_reads:
    """Merge background reads from two paired end fastq into a single file."""
    input:
        in1 = 'qc_reads/1018_R1.fq.gz',
        in2 = 'qc_reads/1018_R2.fq.gz',
    output:
        'merged_background_reads/1018.merged.fa.gz'
    threads: 40
    shell:
        """
        bbmerge.sh \
            in1={input.in1} \
            in2={input.in2} \
            out={output} \
            threads={threads}
        """


rule assemble_background_reads:
    """Assemble background reads."""
    input:
        'merged_background_reads/1018.merged.fa.gz'
    output:
        'assembled_background_reads/1018.assembled.fa.gz'
    log:
        'assembled_background_reads/1018.assembly_stats.txt'
    threads: 40
    shell:
        """
        tadpole.sh \
            in={input} \
            out={output} \
            -Xmx60g \
            > {log}
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
        statsfile = 'logs/filtered_human/{sample}.human_filtered.statsfile.txt'
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
            statsfile={output.statsfile} \
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


rule subtract_background:
    """Subtract/filter (remove) background reads from sample."""
    input:
        read1 = 'filtered_human/{sample}_R1.human_filtered.fq.gz',
        read2 = 'filtered_human/{sample}_R2.human_filtered.fq.gz',
        background_reads = 'assembled_background_reads/1018.assembled.fa.gz'
    output:
        out_read1 = 'subtracted_background/{sample}_R1.background_subtracted.fq.gz',
        out_read2 = 'subtracted_background/{sample}_R2.background_subtracted.fq.gz',
        out_matched = 'subtracted_background/{sample}.background_matched.fq.gz',
        statsfile = 'logs/subtracted_background/{sample}.background_subtracted.statsfile.txt'
    threads: 40
    shell:
        """
        bbmap.sh \
            in1={input.read1} \
            in2={input.read2} \
            outu1={output.out_read1} \
            outu2={output.out_read2} \
            outm={output.out_matched} \
            statsfile={output.statsfile} \
            ref={input.background_reads} \
            build=2 \
            minid=0.95 \
            maxindel=1 \
            minhits=2 \
            quickmatch \
            fast \
            qtrim=rl \
            untrim \
            threads={threads} \
            -Xmx60g
        """


rule merge_subtracted_reads:
    """Merge overlapping paired end reads."""
    input: 
        read1 = 'subtracted_background/{sample}_R1.background_subtracted.fq.gz',
        read2 = 'subtracted_background/{sample}_R2.background_subtracted.fq.gz'
    output:
        'merged_subtracted_reads/{sample}.merged_subtracted.fq.gz'
    threads: 20
    shell:
        """
        bbmerge.sh \
            in1={input.read1} \
            in2={input.read2} \
            out={output} \
            threads={threads}
        """


rule mapDamage:
    """Run mapDamage to assess ancient DNA content."""
    input:
        rules.filter_human.output.out_matched
    output:
        'mapDamage/{sample}/Fragmisincorporation_plot.pdf'
    threads: 1
    shell:
        """
        mapDamage \
            --input {input} \
            --reference {config[hg19_fasta]} \
            --folder mapDamage/{wildcards.sample} 
        """
        

rule metaphlan2_pre:
    """Run MetaPhlAn2 on the QC reads."""
    input:
        reads1 = 'qc_reads/{sample}_R1.fq.gz',
        reads2 = 'qc_reads/{sample}_R2.fq.gz'
    output:
        bowtie2out = 'metaphlan2/bowtie2_outputs/pre/{sample}.bowtie2.bz2',
        mpl_out = 'metaphlan2/metaphlan2_outputs/pre/{sample}.metaphlan2.txt' 
    threads: 20
    shell:
        """
        metaphlan2.py \
            --nproc {threads} \
            --bowtie2out {output.bowtie2out} \
            --input_type fastq \
            --mpa_pkl {config[mpa_pkl]} \
            --bowtie2db {config[mpa_bowtie2db]} \
            --sample_id pre_{wildcards.sample} \
            {input.reads1},{input.reads2} \
            {output.mpl_out}
        """
        

rule metaphlan2:
    """Run MetaPhlAn2 on the merged subtracted reads."""
    input:
        reads = 'merged_subtracted_reads/{sample}.merged_subtracted.fq.gz'
    output:
        bowtie2out = 'metaphlan2/bowtie2_outputs/{sample}.bowtie2.bz2',
        mpl_out = 'metaphlan2/metaphlan2_outputs/{sample}.metaphlan2.txt' 
    threads: 20
    shell:
        """
        metaphlan2.py \
            --nproc {threads} \
            --bowtie2out {output.bowtie2out} \
            --input_type fastq \
            --mpa_pkl {config[mpa_pkl]} \
            --bowtie2db {config[mpa_bowtie2db]} \
            --sample_id {wildcards.sample} \
            {input.reads} \
            {output.mpl_out}
        """
