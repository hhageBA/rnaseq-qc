rule star_pe:
    input:
        fq = expand(
            "Analyses/{{analysis}}/{{run}}/preprocessing/{{sample}}{r}_preprocessed.fastq.gz",
            r=config['read_type']),
    output:
        bam="Analyses/{analysis}/{run}/mapping/{sample}.Aligned.sortedByCoord.out.bam",
        trs_bam=temp("Analyses/{analysis}/{run}/mapping/{sample}.Aligned.toTranscriptome.out.bam"),
        log="Analyses/{analysis}/{run}/mapping/{sample}.Log.final.out",
        dir_tmp=temp(directory("Analyses/{analysis}/{run}/mapping/{sample}._STARgenome")),
    log:
        "Analyses/{analysis}/{run}/logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index = config["star_index"],
        annotation = config["annotation"],
        # optional parameters
        extra = "--outSAMtype BAM SortedByCoordinate"
    threads: 8
    conda:
        "environment.yaml"
    shell:
        "mkdir -p Analyses/{wildcards.analysis}/{wildcards.run}/mapping; "\
            "STAR --runThreadN {threads} --genomeDir {params.index} "\
            "--readFilesIn {input.fq}  "\
            "--sjdbGTFfile {params.annotation} "\
            "--outSAMtype BAM SortedByCoordinate "\
            "--runDirPerm All_RWX "\
            "--quantMode TranscriptomeSAM "\
            "--readFilesCommand zcat "\
            "--outFileNamePrefix Analyses/{wildcards.analysis}/{wildcards.run}/mapping/{wildcards.sample}."

