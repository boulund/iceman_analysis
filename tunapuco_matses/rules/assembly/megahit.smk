# Snakemake rules for metagenome assembly

megahit_output = expand(joinpath(OUTDIR, "megahit", "{sample}", "{sample}.contigs.fa"),
        sample=SAMPLES)
all_outputs.extend(megahit_output)

rule assemble_reads:
    """Assemble reads using MEGAHIT"""
    input:
        read1=joinpath(OUTDIR, "filtered_human", "{sample}_1.filtered_human.fq.gz"), 
        read2=joinpath(OUTDIR, "filtered_human", "{sample}_2.filtered_human.fq.gz"), 
    output:
        contigs = joinpath(OUTDIR, "megahit", "{sample}", "{sample}.contigs.fa"),
    log: joinpath(OUTDIR, "logs", "megahit", "{sample}.megahit.log")
    threads: 20
    conda: "../../envs/iceman.yaml"
    params:
        outdir=joinpath(OUTDIR, "megahit", "{wildcards.sample}")
    shell:
        """
        megahit \
            -1 {input.read1} \
            -2 {input.read2} \
            --out-dir {params.outdir} \
            --out-prefix {wildcards.sample} \
            --num-cpu-threads {threads} \
            2>&1 > {log}
        """
