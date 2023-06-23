rule salmon_combine:
    input:
        quant = lambda w: expand("Analyses/{{analysis}}/{{run}}/salmon/{sample}/quant.genes.sf", sample=RUNS[w.run]["samples"])
    output:
        table = "Analyses/{analysis}/{run}/salmon/{run}_eukaryota_count_table.txt"
    params:
        nb_ech = lambda w: len(RUNS[w.run]["samples"]),
        samples_names = lambda w: RUNS[w.run]["samples"]
    conda:
        "environment.yaml"
    shell:
        "Rscript {PATH}/scripts/combine_quant_files.R {params.nb_ech} {output} {params.samples_names} {input} "
