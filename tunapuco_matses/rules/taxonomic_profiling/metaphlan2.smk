# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using MetaPhlAn2
from snakemake.exceptions import WorkflowError
import os.path

mpa_config = config["metaphlan2"]
bt2_db_ext = ".1.bt2"
if not any([os.path.isfile(mpa_config["mpa_pkl"]),
            os.path.isfile(mpa_config["bt2_db_prefix"]+bt2_db_ext)]):
    err_message = "No MetaPhlAn2 pickle or database found at: '{}', '{}'!\n".format(mpa_config["mpa_pkl"], mpa_config["bt2_db_prefix"])
    err_message += "Specify relevant paths in the metaphlan2 section of config.yaml.\n"
    err_message += "Run 'snakemake build_metaphlan2_index' to download and build the default mpa_v20_m200 database in '{dbdir}/metaphlan2'\n".format(dbdir=config["dbdir"])
    err_message += "If you do not want to run MetaPhlAn2 for taxonomic profiling, set metaphlan2: False in config.yaml"
    raise WorkflowError(err_message)

# Add MetaPhlAn2 output files to 'all_outputs' from the main Snakefile scope.
# SAMPLES is also from the main Snakefile scope.
mpa_outputs = expand(joinpath(OUTDIR, "metaphlan2", "{sample}.{output_type}"),
        sample=SAMPLES,
        output_type=("bowtie2.bz2", "metaphlan2.txt", "metaphlan2.krona"))
mpa_combined = expand(joinpath(OUTDIR, "metaphlan2", "all_samples.metaphlan2.{ext}"),
        ext=("txt", "pdf", "krona.html"))
all_outputs.extend(mpa_outputs)
all_outputs.extend(mpa_combined)

rule download_metaphlan2_database:
    """Download MetaPhlAn2 db_v20_m200"""
    output:
        joinpath(OUTDIR, "db", "metaphlan2", "mpa_v20_m200.fna"),
        joinpath(OUTDIR, "db", "metaphlan2", "mpa_v20_m200.pkl"),
    shadow:
        "shallow"
    conda:
        "../../envs/metaphlan2.yaml"
    params:
        dbdir=joinpath(OUTDIR, "db", "metaphlan2")
    shell:
        """
        cd {params.dbdir}
        wget https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200.tar \
        && \
        tar -xf mpa_v20_m200.tar \
        && \
        bunzip2 mpa_v20_m200.fna.bz2 \
        && \
        rm -v mpa_v20_m200.tar
        """

rule build_metaphlan2_index:
    """Build MetaPhlAn2 bowtie2 index."""
    input:
        joinpath(OUTDIR, "db", "metaphlan2", "mpa_v20_m200.fna")
    output:
        [joinpath(OUTDIR, "db", "metaphlan2", "mpa_v20_m200.{n}.bt2").format(n=num) for num in (1,2,3,4)],
        [joinpath(OUTDIR, "db", "metaphlan2", "mpa_v20_m200.rev.{n}.bt2").format(n=num) for num in (1,2)],
    log: joinpath(OUTDIR, "logs", "metaphlan2", "bowtie2-build.log")
    shadow: "shallow"
    conda: "../../envs/metaphlan2.yaml"
    threads: 4
    params:
        dbdir=joinpath(OUTDIR, "db", "metaphlan2")
    shell:
        """
        cd {params.dbdir}
        bowtie2-build \
            mpa_v20_m200.fna \
            mpa_v20_m200 \
            --threads {threads} \
            2>&1 > {log}
        """

rule metaphlan2:
    """Taxonomic profiling using MetaPhlAn2."""
    input:
        read1=joinpath(OUTDIR, "filtered_human", "{sample}_1.filtered_human.fq.gz"),
        read2=joinpath(OUTDIR, "filtered_human", "{sample}_2.filtered_human.fq.gz"),
    output:
        bt2_out=joinpath(OUTDIR, "metaphlan2", "{sample}.bowtie2.bz2"),
        mpa_out=joinpath(OUTDIR, "metaphlan2", "{sample}.metaphlan2.txt"),
        krona=joinpath(OUTDIR, "metaphlan2", "{sample}.metaphlan2.krona"),
    log:
        stdout=joinpath(OUTDIR, "logs", "metaphlan2", "{sample}.metaphlan2.stdout.log"),
        stderr=joinpath(OUTDIR, "logs", "metaphlan2", "{sample}.metaphlan2.stderr.log"),
    shadow: "shallow"
    conda: "../../envs/metaphlan2.yaml"
    threads: 4
    params:
        mpa_pkl=mpa_config["mpa_pkl"],
        bt2_db_prefix=mpa_config["bt2_db_prefix"],
        extra=mpa_config["extra"],
    shell:
        """
        metaphlan2.py \
            --input_type fastq \
            --nproc {threads} \
            --sample_id {wildcards.sample} \
            --bowtie2out {output.bt2_out} \
            --mpa_pkl {params.mpa_pkl} \
            --bowtie2db {params.bt2_db_prefix} \
            {input.read1},{input.read2} \
            {output.mpa_out} \
            {params.extra} \
            > {log.stdout} \
            2> {log.stderr} \
        && \
        metaphlan2krona.py \
            --profile {output.mpa_out} \
            --krona {output.krona}
        """


rule combine_metaphlan2_outputs:
    """Combine metaphlan2 outputs into a large table and plot heatmap."""
    input:
        expand(joinpath(OUTDIR, "metaphlan2", "{sample}.metaphlan2.txt"), sample=SAMPLES)
    output:
        txt=joinpath(OUTDIR, "metaphlan2", "all_samples.metaphlan2.txt"),
        pdf=joinpath(OUTDIR, "metaphlan2", "all_samples.metaphlan2.pdf"),
    shadow: "shallow"
    conda: "../../envs/metaphlan2.yaml"
    threads: 1
    params:
        tax_lev=mpa_config["tax_lev"],
        minv=mpa_config["minv"],
        s=mpa_config["s"],
        top=mpa_config["top"],
        m=mpa_config["m"],
        d=mpa_config["d"],
        f=mpa_config["f"],
        c=mpa_config["c"],
    shell:
        """
        merge_metaphlan_tables.py {input} \
            | sed 's/\.metaphlan2//g' > {output.txt} \
        && \
        metaphlan_hclust_heatmap.py \
            --in {output.txt} \
            --out {output.pdf} \
            --tax_lev {params.tax_lev} \
            --minv {params.minv} \
            -s {params.s} \
            --top {params.top} \
            -m {params.m} \
            -d {params.d} \
            -f {params.f} \
            -c {params.c}
        """

rule create_metaphlan2_krona_plots:
    input:
        expand(joinpath(OUTDIR, "metaphlan2", "{sample}.metaphlan2.krona"), sample=SAMPLES)
    output:
        html=joinpath(OUTDIR, "metaphlan2", "all_samples.metaphlan2.krona.html"),
    shadow: "shallow"
    conda: "../../envs/metaphlan2.yaml"
    threads: 1
    shell:
        """
        ktImportText \
            -o {output.html} \
            {input}
        """
