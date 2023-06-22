RNASEQ-QCS is a bioinformatics pipeline for the quality control analysis of RNA sequencing data obtained from organism with an available reference genome and annotation.

STAGES:
    - preprocessing:
        - fastp : Read filtering
        - FASTQC : Read QC
        - SortMeRNA : ribosomal RNA quantification
        - bowtie2 : globin quantification
