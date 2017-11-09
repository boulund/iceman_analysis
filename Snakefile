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
        expand('microbecensus/{sample}.microbecensus.txt', sample=config["samples"]),
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
        expand('logs/input_read_quality/{sample}.input_read_quality.stats.txt', sample=config["samples"]),
        expand('logs/input_read_quality/{background_sample}.input_read_quality.stats.txt', background_sample=config["background_sample"]),
        expand('megares/{sample}.megares.covstats.txt', sample=config["samples"]),
        expand('bbcountunique/{sample}.bbcountunique.histogram.txt', sample=config["samples"]),
        expand('bbcountunique/{sample}.bbcountunique.histogram.pdf', sample=config["samples"]),
        expand('PMDtools/{sample}.deamination.pdf', sample=config["samples"]),


rule assess_input_data_quality:
    """Run reads through BBDuk for quality stats."""
    input:
        in1 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R1.fastq.gz',
        in2 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R2.fastq.gz',
    output:
        stats =  'logs/input_read_quality/{sample}.input_read_quality.stats.txt',
        bhist =  'logs/input_read_quality/{sample}.input_read_quality.bhist.txt',
        qhist =  'logs/input_read_quality/{sample}.input_read_quality.qhist.txt',
        qchist = 'logs/input_read_quality/{sample}.input_read_quality.qchist.txt',
        aqhist = 'logs/input_read_quality/{sample}.input_read_quality.aqhist.txt',
        bqhist = 'logs/input_read_quality/{sample}.input_read_quality.bqhist.txt',
        lhist =  'logs/input_read_quality/{sample}.input_read_quality.lhist.txt',
        gchist = 'logs/input_read_quality/{sample}.input_read_quality.gchist.txt',
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


rule average_genome_size:
    """Run MicrobeCensus to estimate Average Genome Size and number of 
    genome equivalents for use in normalization."""
    input:
        in1 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R1.fastq.gz',
        in2 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R2.fastq.gz',
    output:
        'microbecensus/{sample}.microbecensus.txt'
    log:
        'logs/microbecensus/{sample}.microbecensus_log.txt'
    threads: 8
    shell:
        """
        run_microbe_census.py \
            -t {threads} \
            -v \
            {input.in1},{input.in2} \
            {output} \
            > {log}
        """


rule qc_reads:
    """Run reads through BBDuk for QC and adapter trimming."""
    input:
        in1 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R1.fastq.gz',
        in2 = config['samples_dir']+'{sample}_Iceman/{sample}_Iceman_R2.fastq.gz',
    output:
        out1 = 'qc_reads/{sample}_R1.fq.gz',
        out2 = 'qc_reads/{sample}_R2.fq.gz',
        stats =  'logs/qc_reads/{sample}.qc_reads.stats.txt',
        bhist =  'logs/qc_reads/{sample}.qc_reads.bhist.txt',
        qhist =  'logs/qc_reads/{sample}.qc_reads.qhist.txt',
        qchist = 'logs/qc_reads/{sample}.qc_reads.qchist.txt',
        aqhist = 'logs/qc_reads/{sample}.qc_reads.aqhist.txt',
        bqhist = 'logs/qc_reads/{sample}.qc_reads.bqhist.txt',
        lhist =  'logs/qc_reads/{sample}.qc_reads.lhist.txt',
        gchist = 'logs/qc_reads/{sample}.qc_reads.gchist.txt',
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


rule assess_saturation:
    """Count unique kmers to assess sequencing depth/saturation."""
    input:
        in1 = 'qc_reads/{sample}_R1.fq.gz',
        in2 = 'qc_reads/{sample}_R2.fq.gz'
    output:
        histogram_data = 'bbcountunique/{sample}.bbcountunique.histogram.txt',
        histogram_plot = 'bbcountunique/{sample}.bbcountunique.histogram.pdf',
    log:
        'logs/bbcountunique/{sample}.bbcountunique.stdout'
    shell:
        """
        bbcountunique.sh \
            in1={input.in1} \
            in2={input.in2} \
            out={output.histogram_data} \
            > {log}
        plot_bbcountunique.py \
            {output.histogram_data} \
            {output.histogram_plot} \
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
        'logs/assembled_background_reads/1018.assembly_stats.txt'
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
    log:
        'logs/filtered_human/{sample}.human_filtered.statsfile.txt'
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
            statsfile={log} \
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
    log:
        'logs/subtracted_background/{sample}.background_subtracted.statsfile.txt'
    threads: 40
    shell:
        """
        bbmap.sh \
            in1={input.read1} \
            in2={input.read2} \
            outu1={output.out_read1} \
            outu2={output.out_read2} \
            outm={output.out_matched} \
            statsfile={log} \
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
    log:
        'logs/merged_subtracted_reads/{sample}.merged_subtracted.stdout'
    threads: 20
    shell:
        """
        bbmerge.sh \
            in1={input.read1} \
            in2={input.read2} \
            out={output} \
            threads={threads} \
            > {log}
        """


rule mapDamage:
    """Run mapDamage to assess ancient DNA content."""
    input:
        rules.filter_human.output.out_matched
    output:
        'mapDamage/{sample}/Fragmisincorporation_plot.pdf'
    log:
        'logs/mapDamage/{sample}.mapDamage.stdout'
    threads: 1
    shell:
        """
        mapDamage \
            --input {input} \
            --reference {config[hg19_fasta]} \
            --folder mapDamage/{wildcards.sample} \
            > {log}
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


rule map_antibiotic_resistance:
    """Map quality filtered and subtracted reads to MEGARes."""
    input:
        reads = 'merged_subtracted_reads/{sample}.merged_subtracted.fq.gz'
    output:
       samfile = 'megares/{sample}.megares.sam',
       scafstats = 'megares/{sample}.megares.scafstats.txt',
       covstats = 'megares/{sample}.megares.covstats.txt',
       rpkm = 'megares/{sample}.megares.rpkm.txt',
       basecov = 'megares/{sample}.megares.basecov.txt',
       dataframe = 'megares/{sample}.megares.covstats.dataframe.csv',
    log:
       statsfile = 'logs/megares/{sample}.megares.statsfile.txt',
       count_matrix = 'logs/megares/{sample}.megares.dataframe.txt',
    threads: 40
    shell:
        """
        bbmap.sh \
            ref={config[megares_fasta]} \
            build=3 \
            in={input.reads} \
            minid=0.90 \
            threads={threads} \
            outm={output.samfile} \
            scafstats={output.scafstats} \
            covstats={output.covstats} \
            rpkm={output.rpkm} \
            basecov={output.basecov} \
            statsfile={log.statsfile} \
        create_count_matrix.py \
            {output.covstats} \
            -o {output.dataframe} \
            -a {config[megares_annotations]} \
            > {log.count_matrix} 
        """

rule pmd_tools_ar:
    """Use PMD tools to remove modern contamination from AR-mappings."""
    input:
        samfile = 'megares/{sample}.megares.sam'
    output:
        fillmd = 'PMDtools/{sample}.megares.fillmd.sam',
        bamfile = 'PMDtools/{sample}.pmds3filter.bam',
        deamination = 'PMDtools/{sample}.deamination.txt',
        plot = 'PMDtools/{sample}.deamination.pdf',
    shell:
        """
        samtools fillmd {input.samfile} {config[megares_fasta]} > {output.fillmd}
        
        samtools view -h {output.fillmd} \
            | {config[pmdtools_path]} --threshold 3 --header \
            | samtools view -Sb - \
            > {output.bamfile}

        samtools view {output.bamfile} \
            | {config[pmdtools_path]} --deamination --range 30 \
            > {output.deamination}

        plot_pmd.py {output.deamination} {output.plot}
        """
