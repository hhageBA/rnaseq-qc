import os

rule generate_sample_correspondances:
    input:
        meta = config["metadata"]
    output:
        corresp= os.path.join(LAUNCHED_SMK, "sample_correspondances.txt")
    params:
        headers=["SampleID", "SampleName"]
    run:
        with open(output.corresp, "w") as fileout:
            fileout.write("{s}\t{b}\n".format(s=params.headers[0], b=params.headers[1]))
            for spl in SAMPLES.keys():
                fileout.write("{s}\t{b}\n".format(s=SAMPLES[spl]["demultiplexed"], b=spl))



def get_rawdata_fastqc(wildcards):
    listfiles=[]
    for spl in SAMPLES.keys():
        if SAMPLES[spl]["RunID"]== wildcards.run :
            spl_id = SAMPLES[spl]["demultiplexed"]
            pR1="Analyses/{analysis}/{run}/qc/{spl_id}{r1}_fastqc.zip".format(analysis=wildcards.analysis, run=wildcards.run, spl_id=spl_id, r1=config['read_type'][0])
            pR2="Analyses/{analysis}/{run}/qc/{spl_id}{r2}_fastqc.zip".format(analysis=wildcards.analysis, run=wildcards.run, spl_id=spl_id, r2=config['read_type'][1])
            listfiles.append(pR1)
            listfiles.append(pR2)
    return(listfiles)

rule multiqc:
    input:
        raw = get_rawdata_fastqc,
        filt= lambda w: expand("Analyses/{{analysis}}/{{run}}/preprocessing/FastQC/{sample}{r}_preprocessed_fastqc.zip",
            r=config["read_type"], sample=RUNS[w.run]["samples"]),
        inner_distance = lambda w: expand("Analyses/{{analysis}}/{{run}}/mapping/{sample}.inner_distance.txt", sample=RUNS[w.run]["samples"]),
        infer_exp = lambda w: expand("Analyses/{{analysis}}/{{run}}/mapping/{sample}.infer_exp.txt", sample=RUNS[w.run]["samples"]),
        read_distrib = lambda w: expand("Analyses/{{analysis}}/{{run}}/mapping/{sample}_read_distrib.txt", sample=RUNS[w.run]["samples"]),
        quantif = lambda w: expand("Analyses/{{analysis}}/{{run}}/salmon/{sample}/quant.genes.sf", sample=RUNS[w.run]["samples"]),
        sortme = lambda w: expand("Analyses/{{analysis}}/{{run}}/preprocessing/{sample}.rrna.out.log", sample=RUNS[w.run]["samples"]),
        corresp_names = os.path.join(LAUNCHED_SMK, "sample_correspondances.txt"),
        additionalqc = "Analyses/{analysis}/{run}/additionalQC/Number_of_detected_genes_mqc.png"
    output:
        html = "Analyses/{analysis}/{run}/Reporting/{run}_eukaryota_multiqc.html",
    conda:
        "environment.yaml"
    params:
        multiqc_conf = os.path.join(PATH, "multiqc_rnaseq_config.yaml"),
        out_dir = lambda w: "Analyses/{}/{}".format(w.analysis, w.run),
        desc_sentence = config["descriptive_sentence"],
       
    log:
        "Analyses/{analysis}/{run}/multiqc.log"
    shell:
        """
        multiqc -c {params.multiqc_conf} \
        -i "{wildcards.run} : MultiQC analysis" \
        --replace-names {input.corresp_names} \
        -n {output.html} \
        {params.out_dir}
        """
