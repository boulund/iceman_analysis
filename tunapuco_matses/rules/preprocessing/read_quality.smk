# Snakemake rules for read quality assessment

fastqc_output = expand("{outdir}/fastqc/{sample}_{readpair}_fastqc.{ext}",
        outdir=OUTDIR,
        sample=SAMPLES,
        readpair=[1,2],
        ext=["zip", "html"])
trimmed_qa = expand("{outdir}/qc_reads/{sample}_{readpair}.fq.gz",
        outdir=OUTDIR,
        sample=SAMPLES,
        readpair=[1,2])
kmer_saturation = expand(joinpath(OUTDIR, 'bbcountunique', '{sample}.bbcountunique.histogram.pdf'),
        sample=SAMPLES)
all_outputs.extend(fastqc_output)
all_outputs.extend(trimmed_qa)
all_outputs.extend(kmer_saturation)


rule fastqc:
    input:
        joinpath(config["inputdir"], config["input_fn_pattern"])
    output:
        html=joinpath(OUTDIR, "fastqc", "{sample}_{readpair}_fastqc.html"),
        zip=joinpath(OUTDIR, "fastqc", "{sample}_{readpair}_fastqc.zip"),
    log: joinpath(OUTDIR, "logs", "fastqc", "{sample}_{readpair}.log")
    shadow: 
        "shallow"
    wrapper:
        "0.22.0/bio/fastqc"


rule qc_reads:
    """Run reads through BBDuk for QC and adapter trimming."""
    input:
        read1=joinpath(INDIR, config["input_fn_pattern"]).format(sample="{sample}", readpair="1"),
        read2=joinpath(INDIR, config["input_fn_pattern"]).format(sample="{sample}", readpair="2")
    output:
        out1 = joinpath(OUTDIR, 'qc_reads', '{sample}_1.fq.gz'),
        out2 = joinpath(OUTDIR, 'qc_reads', '{sample}_2.fq.gz'),
        stats =  'logs/qc_reads/{sample}.qc_reads.stats.txt',
        bhist =  'logs/qc_reads/{sample}.qc_reads.bhist.txt',
        qhist =  'logs/qc_reads/{sample}.qc_reads.qhist.txt',
        qchist = 'logs/qc_reads/{sample}.qc_reads.qchist.txt',
        aqhist = 'logs/qc_reads/{sample}.qc_reads.aqhist.txt',
        bqhist = 'logs/qc_reads/{sample}.qc_reads.bqhist.txt',
        lhist =  'logs/qc_reads/{sample}.qc_reads.lhist.txt',
        gchist = 'logs/qc_reads/{sample}.qc_reads.gchist.txt',
    log: joinpath(OUTDIR, "logs", "qc_reads", "{sample}.bbduk.log")
    threads: 20
    conda: "../../envs/iceman.yaml"
    shell:
        """
        bbduk.sh \
            in1={input.read1} \
            in2={input.read2} \
            ref=adapters \
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
            hdist=1 \
            2> {log}
        """

rule assess_saturation:
    """Count unique kmers to assess sequencing depth/saturation."""
    input:
        in1 = joinpath(OUTDIR, 'qc_reads', '{sample}_1.fq.gz'),
        in2 = joinpath(OUTDIR, 'qc_reads', '{sample}_2.fq.gz'),
    output:
        histogram_data = joinpath(OUTDIR, 'bbcountunique', '{sample}.bbcountunique.histogram.txt'),
        histogram_plot = joinpath(OUTDIR, 'bbcountunique', '{sample}.bbcountunique.histogram.pdf'),
    log:
        joinpath(OUTDIR, 'logs/bbcountunique/{sample}.bbcountunique.stdout')
    conda: "../../envs/iceman.yaml"
    threads: 8
    shell:
        """
        bbcountunique.sh \
            in1={input.in1} \
            in2={input.in2} \
            out={output.histogram_data} \
            > {log}
        scripts/plot_bbcountunique.py \
            {output.histogram_data} \
            {output.histogram_plot} \
        """

