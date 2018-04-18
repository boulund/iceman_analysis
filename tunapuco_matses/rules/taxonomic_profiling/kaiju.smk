# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kaiju
from snakemake.exceptions import WorkflowError
import os.path

kaiju_config = config["kaiju"]
if not all([os.path.isfile(kaiju_config["db"]), 
            os.path.isfile(kaiju_config["nodes"]),
            os.path.isfile(kaiju_config["names"])]):
    err_message = "No Kaiju database files at: '{}', '{}', '{}'!\n".format(kaiju_config["db"], kaiju_config["nodes"], kaiju_config["names"])
    err_message += "Specify relevant paths in the kaiju section of config.yaml.\n"
    err_message += "Run 'snakemake download_kaiju_database' to download a copy into '{dbdir}/kaiju'\n".format(dbdir=joinpath(OUTDIR, "db"))
    err_message += "If you do not want to run Kaiju for taxonomic profiling, set 'kaiju: False' in config.yaml"
    raise WorkflowError(err_message)

# Add Kaiju output files to 'all_outputs' from the main Snakefile scope.
# SAMPLES is also from the main Snakefile scope.
kaiju = expand("{outdir}/kaiju/{sample}.kaiju", outdir=OUTDIR, sample=SAMPLES)
kaiju_reports = expand("{outdir}/kaiju/{sample}.kaiju.summary.species", outdir=OUTDIR, sample=SAMPLES)
kaiju_krona = expand("{outdir}/kaiju/all_samples.kaiju.krona.html", outdir=OUTDIR)
all_outputs.extend(kaiju)
all_outputs.extend(kaiju_reports)
all_outputs.extend(kaiju_krona)

rule download_kaiju_database:
    output:
        db=joinpath(OUTDIR, "kaiju", "kaiju_db.fmi"),
        names=joinpath(OUTDIR, "kaiju", "names.dmp"),
        nodes=joinpath(OUTDIR, "kaiju", "nodes.dmp"),
    shadow:
        "shallow"
    params:
        dbdir=joinpath(OUTDIR, "db", "kaiju")
    shell:
        """
        wget http://kaiju.binf.ku.dk/database/kaiju_index_pg.tgz \
        && \
        tar -xf kaiju_index_pg.tgz \
        && \
        mv kaiju_db.fmi names.dmp nodes.dmp {params.dbdir}
        """

rule kaiju:
    input:
        read1=joinpath(OUTDIR, "filtered_human", "{sample}_1.filtered_human.fq.gz"),
        read2=joinpath(OUTDIR, "filtered_human", "{sample}_2.filtered_human.fq.gz"),
    output:
        kaiju=joinpath(OUTDIR, "kaiju", "{sample}.kaiju"),
    shadow: "shallow"
    threads: 4
    conda: "../../envs/iceman.yaml"
    params:
        db=kaiju_config["db"],
        nodes=kaiju_config["nodes"],
    shell:
        """
        if [[ "{input.read1}" == *.gz ]]
        then
            kaiju \
                -z {threads} \
                -t {params.nodes} \
                -f {params.db} \
                -i <(gunzip -c {input.read1}) \
                -j <(gunzip -c {input.read2}) \
                -o {output.kaiju} 
        else 
            kaiju \
                -z {threads} \
                -t {params.nodes} \
                -f {params.db} \
                -i {input.read1} \
                -j {input.read2} \
                -o {output.kaiju}
        fi
        """


rule kaiju_report:
    input:
        kaiju=joinpath(OUTDIR, "kaiju", "{sample}.kaiju"),
    output:
        krona=joinpath(OUTDIR, "kaiju", "{sample}.krona"),
        family=joinpath(OUTDIR, "kaiju", "{sample}.summary.family"),
        genus=joinpath(OUTDIR, "kaiju", "{sample}.kaiju.summary.genus"),
        species=joinpath(OUTDIR, "kaiju", "{sample}.kaiju.summary.species"),
    shadow: "shallow"
    conda: "../../envs/iceman.yaml"
    params:
        nodes=kaiju_config["nodes"],
        names=kaiju_config["names"],
    shell:
        """
		kaiju2krona \
			-t {params.nodes} \
			-n {params.names} \
			-i {input.kaiju} \
			-o {output.krona} \
            -u
        kaijuReport \
            -t {params.nodes} \
            -n {params.names} \
            -i {input.kaiju} \
            -r species \
            -p \
            -o {output.species}
        kaijuReport \
            -t {params.nodes} \
            -n {params.names} \
            -i {input.kaiju} \
            -r genus \
            -p \
            -o {output.genus}
        kaijuReport \
            -t {params.nodes} \
            -n {params.names} \
            -i {input.kaiju} \
            -r family \
            -p \
            -o {output.family}
        """

rule create_kaiju_krona_plot:
    input:
        expand(joinpath(OUTDIR, "kaiju", "{sample}.krona"), sample=SAMPLES)
    output:
        krona_html=joinpath(OUTDIR, "kaiju", "all_samples.kaiju.krona.html"),
    shadow: "shallow"
    conda: "../../envs/iceman.yaml"
    shell:
        """
		ktImportText \
			-o {output.krona_html} \
			{input}
        """
