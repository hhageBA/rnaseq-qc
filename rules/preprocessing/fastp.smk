def get_fastqR1(wildcards):
    return([SAMPLES[wildcards.sample]['R1']])

def get_fastqR2(wildcards):
    return([SAMPLES[wildcards.sample]['R2']])

rule fastp_pe:
    input:
        sampleR1= get_fastqR1,
        sampleR2= get_fastqR2
    output:
        trimmedR1 = [expand(
                    "Analyses/{{analysis}}/{{run}}/preprocessing/{{sample}}{r}_preprocessed.fastq.gz",
                    r=config["read_type"][0])],
        trimmedR2 = [expand(
                    "Analyses/{{analysis}}/{{run}}/preprocessing/{{sample}}{r}_preprocessed.fastq.gz",
                    r=config["read_type"][1])],
        html = "Analyses/{analysis}/{run}/preprocessing/{sample}_fastp.html",
        json = "Analyses/{analysis}/{run}/preprocessing/{sample}_fastp.json"
    log:
        "Analyses/{analysis}/{run}/logs/fastp/{sample}.log"
    conda:
        "environment_fastp.yaml"
    params:
        extra = "-z 9",
        adapter_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    threads: 4
    shell:
        """
        fastp -z 9 -i {input.sampleR1} -I {input.sampleR2} \
        --json {output.json} --html {output.html} \
        --adapter_sequence={params.adapter_R1} --adapter_sequence_r2={params.adapter_R2} \
        -o {output.trimmedR1} -O {output.trimmedR2}
        """

