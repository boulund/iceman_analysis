
mapdamage_plots = expand(joinpath(OUTDIR, "mapdamage", "{sample}", "Fragmisincorporation_plot.pdf"),
        sample=SAMPLES)
all_outputs.extend(mapdamage_plots)

rule mapDamage:
    """Run mapDamage to assess ancient DNA content."""
    input:
        joinpath(OUTDIR, "filtered_human", "{sample}_human.fq.gz")
    output:
        joinpath(OUTDIR, "mapdamage", "{sample}", "Fragmisincorporation_plot.pdf")
    log:
        joinpath(OUTDIR, "logs", "mapdamage", "{sample}.mapDamage.stdout")
    conda: "../../envs/iceman.yaml"
    shadow: "shallow"
    threads: 1
    shell:
        """
        mapDamage \
            --input {input} \
            --reference {config[hg19_fasta]} \
            --folder mapDamage/{wildcards.sample} \
            > {log}
        """
