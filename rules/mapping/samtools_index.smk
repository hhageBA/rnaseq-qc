rule samtools_index:
    input:
        "Analyses/{analysis}/{run}/mapping/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "Analyses/{analysis}/{run}/mapping/{sample}.Aligned.sortedByCoord.out.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.35.2/bio/samtools/index"
