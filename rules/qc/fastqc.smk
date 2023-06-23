def get_rawdata(wildcards):
    listfiles=[]
    for spl in SAMPLES.keys():
        #print(wildcards.spl_id)
        if SAMPLES[spl]["demultiplexed"]== wildcards.spl_id :
            return([SAMPLES[spl]["R1"], SAMPLES[spl]["R2"]])

rule fastqc:
    input:
        fq = get_rawdata,
   
    output:
        html = [expand(
                "Analyses/{{analysis}}/{{run}}/qc/{{spl_id}}{r1}_fastqc.html",
                r1=config["read_type"])],
        zip = [expand(
                "Analyses/{{analysis}}/{{run}}/qc/{{spl_id}}{r1}_fastqc.zip",
                r1=config["read_type"])],
    
    params:
        outdir= "Analyses/{analysis}/{run}/qc/"
    conda:"environment.yaml"
    shell:
        """
        fastqc --outdir {params.outdir} {input.fq}
        """
    
rule fastqc_filtered:
    input:
        fq = "Analyses/{analysis}/{run}/preprocessing/{sample}_preprocessed.fastq.gz",
    output:
        html = "Analyses/{analysis}/{run}/preprocessing/FastQC/{sample}_preprocessed_fastqc.html",
        zip = "Analyses/{analysis}/{run}/preprocessing/FastQC/{sample}_preprocessed_fastqc.zip",
    params:
        outdir=  "Analyses/{analysis}/{run}/preprocessing/FastQC/"
    conda:"environment.yaml"    
    shell:
        """
        fastqc --outdir {params.outdir} {input.fq}
        """