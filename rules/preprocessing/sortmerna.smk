rule seqtk:
    input:
        fq1 = expand(
            "Analyses/{{analysis}}/{{run}}/preprocessing/{{sample}}{r1}_preprocessed.fastq.gz",
            r1=config['read_type'][0]),
    output:
        fastq = temp(expand(
            "{tmp}/{{analysis}}/{{run}}/preprocessing/{{sample}}_sub.fastq",
            tmp=config["tmp_dir"]))
    threads: 1
    conda:
        "environment_sortmerna.yaml"
    shell:
        "seqtk sample {input.fq1} 1000000 > {output.fastq}; "


rule sortmerna:
    input:
        fastq = expand(
            "{tmp}/{{analysis}}/{{run}}/preprocessing/{{sample}}_sub.fastq",
            tmp=config["tmp_dir"])
    output:
        log = "Analyses/{analysis}/{run}/preprocessing/{sample}.rrna.out.log"
    conda:
        "environment_sortmerna.yaml"
    params:
        sortrna = config["tmp_sortmerna"],
        db_rRNA = config["db_rRNA"],
        prefix = "Analyses/{analysis}/{run}/preprocessing/{sample}.rrna.out"
    threads: 20
    shell:
        "rm -r -f {params.sortrna}; "
        "sortmerna " \
            "-a {threads} " \
            "--ref {params.db_rRNA} " \
            "--reads {input.fastq} " \
            "--aligned {params.prefix} " 




