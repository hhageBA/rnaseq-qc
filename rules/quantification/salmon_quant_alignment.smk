rule salmon_quant_alignment:
    input:
        bam = "Analyses/{analysis}/{run}/mapping/{sample}.Aligned.toTranscriptome.out.bam",
    output:
        quant = "Analyses/{analysis}/{run}/salmon/{sample}/quant.genes.sf",
    log:
        "Analyses/{analysis}/{run}/logs/salmon/{sample}.log"
    params:
        # optional parameters
        transcriptome = config["transcripts"],
        annot = config["annotation"],
        libtype = "A",
    threads: 2
    conda:
        "environment.yaml"
    shell:
        "salmon quant "\
            "-t {params.transcriptome} "\
            "-l {params.libtype} "\
            "-a {input.bam} "\
            "-g {params.annot} "\
            "-o Analyses/{wildcards.analysis}/{wildcards.run}/salmon/{wildcards.sample}"
