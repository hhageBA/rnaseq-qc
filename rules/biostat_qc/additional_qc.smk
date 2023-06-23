rule additional_QC:
    input:
        quant = "Analyses/{analysis}/{run}/salmon/{run}_eukaryota_count_table.txt",
        stat = "Analyses/{analysis}/{run}/qc/all_eukaryota_stats_per_samples.txt"
    output:
        fig = "Analyses/{analysis}/{run}/additionalQC/Number_of_detected_genes_mqc.png"
    params:
        metadata = config["metadata"],
        typeFile = config["gene_type_corresp"],
        groupColor = config["condition"],
        Filter_threshold = config["Filter_threshold"]
    conda:
        "environment.yaml"
    shell:
        "Rscript {PATH}/scripts/Stats_scripts/additional_qc.r "\
        "--inputData {input.quant} " \
        "--metadata {params.metadata} " \
        "--statFile {input.stat} " \
        "--typeFile {params.typeFile} " \
        "--groupColor {params.groupColor} " \
        "--filterThreshold {params.Filter_threshold} " \
        "--output {output.fig} " \
        "--scripts {PATH} "
        
   
