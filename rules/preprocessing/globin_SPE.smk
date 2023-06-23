if len(config["read_type"])==2:
    rule globin_quantification_PE:
        input:
            fq1 = expand(
                "Analyses/{{analysis}}/{{run}}/preprocessing/{{sample}}{r1}_preprocessed.fastq.gz",
                r1=config['read_type'][0]),
            fq2 = expand(
                "Analyses/{{analysis}}/{{run}}/preprocessing/{{sample}}{r2}_preprocessed.fastq.gz",
                r2=config['read_type'][1])
        output:
            bam = temp(expand("{tmp}/{{analysis}}/{{run}}/preprocessing/{{sample}}.globin.bam", tmp=config["tmp_dir"])),
            txt = "Analyses/{analysis}/{run}/preprocessing/{sample}.globin.log"
        log:
            "Analyses/{analysis}/{run}/logs/bowtie2/{sample}.log"
        params:
            index = config['globin_index'],
            max_frag_length = 800
        threads: 4
        conda:
            "environment_bowtie.yaml"
        shell:
            """
            bowtie2 -p {threads} --very-sensitive-local -x {params.index} \
              -1 {input.fq1} -2 {input.fq2} -X {params.max_frag_length} \
              2> {output.txt} | samtools  view -bS - |  samtools sort -o {output.bam}
            """
else:
    rule globin_quantification_SE:
        input:
            fq1 = expand(
                "Analyses/{{analysis}}/{{run}}/preprocessing/{{sample}}{r1}_preprocessed.fastq.gz",
                r1=config['read_type'][0]),
        output:
            bam = temp(expand("{tmp}/{{analysis}}/{{run}}/preprocessing/{{sample}}.globin.bam", tmp=config["tmp_dir"])),
            txt = "Analyses/{analysis}/{run}/preprocessing/{sample}.globin.log"
        log:
            "Analyses/{analysis}/{run}/logs/bowtie2/{sample}.log"
        params:
            index = config['globin_index'],
            max_frag_length = 800
        threads: 1
        conda:
            "environment_bowtie.yaml"
        shell:
            """
            bowtie2 -p {threads} --very-sensitive-local -x {params.index} -U {input.fq1}  \
            -X {params.max_frag_length} \
            2> {output.txt} | samtools  view -bS - |  samtools sort -o {output.bam}
            """
