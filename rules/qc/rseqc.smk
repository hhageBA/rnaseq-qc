rule rseqc_gb_coverage:
    input:
        bam = "Analyses/{analysis}/{run}/mapping/{name}.Aligned.sortedByCoord.out.bam",
        index = "Analyses/{analysis}/{run}/mapping/{name}.Aligned.sortedByCoord.out.bam.bai"
    output:
        "Analyses/{analysis}/{run}/mapping/{name}.geneBodyCoverage.curves.jpeg"
    params:
        bed = config["bed"]
    log:
        "Analyses/{analysis}/{run}/logs/rseqc/{name}.geneBodyCoverage.log"
    conda:
        "environment.yaml"
    shell:
        "geneBody_coverage.py "\
            "-i {input.bam} "\
            "-r {params.bed} "\
            "-o Analyses/{wildcards.analysis}/{wildcards.run}/mapping/{wildcards.name} "\
            "-f jpeg"


rule list_bamfiles_for_gb_coverage:
    input:
        bam = lambda w: expand("Analyses/{{analysis}}/{{run}}/mapping/{name}.Aligned.sortedByCoord.out.bam", name=RUNS[w.run]["samples"]),
        index = lambda w: expand("Analyses/{{analysis}}/{{run}}/mapping/{name}.Aligned.sortedByCoord.out.bam.bai", name=RUNS[w.run]["samples"])
    output:
        "Analyses/{analysis}/{run}/mapping/BAMfiles_list_EUK.txt"
    shell:
        """
        echo {input.bam} | tr " " "\n" > {output}
        """

rule rseqc_read_distribution:
    input:
        bam = "Analyses/{analysis}/{run}/mapping/{name}.Aligned.sortedByCoord.out.bam"
    output:
        "Analyses/{analysis}/{run}/mapping/{name}_read_distrib.txt"
    params:
        bed = config["bed"]
    log:
        "Analyses/{analysis}/{run}/logs/rseqc/{name}.readDistrib.log"
    conda:
        "environment.yaml"
    shell:
        "read_distribution.py -i {input} -r {params.bed} &> {output}"


rule rseqc_read_duplication:
    input:
        "Analyses/{analysis}/{run}/mapping/{name}.Aligned.sortedByCoord.out.bam"
    output:
        "Analyses/{analysis}/{run}/mapping/{name}.read_duplication.pos.DupRate.xls"
    params:
        prefix= "Analyses/{analysis}/{run}/mapping/{name}.read_duplication"
    log:
        "Analyses/{analysis}/{run}/logs/rseqc/{name}.readDuplication.log"
    conda:
        "environment.yaml"
    shell :
        "read_duplication.py -i {input} -o {params.prefix}"


rule rseqc_inner_distance:
    input: "Analyses/{analysis}/{run}/mapping/{name}.Aligned.sortedByCoord.out.bam"
    output: "Analyses/{analysis}/{run}/mapping/{name}.inner_distance.txt"
    params:
        bedfile=config["bed"],
        prefix = "Analyses/{analysis}/{run}/mapping/{name}",
    log:
        "Analyses/{analysis}/{run}/logs/rseqc/{name}.innerDistance.log"
    conda:
        "environment.yaml"
    shell :
        "inner_distance.py -i {input} -r {params.bedfile} -o {params.prefix}"


rule infer_exp:
    input: "Analyses/{analysis}/{run}/mapping/{name}.Aligned.sortedByCoord.out.bam"
    output: "Analyses/{analysis}/{run}/mapping/{name}.infer_exp.txt"
    params:
         bedfile=config["bed"]
    log:
        "Analyses/{analysis}/{run}/logs/rseqc/{name}.inferExperiment.log"
    conda:
        "environment.yaml"
    shell :
        "infer_experiment.py -i {input} -r {params.bedfile} &> {output}"
