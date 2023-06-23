RNASEQ-QCS is a bioinformatics pipeline that can be used for the quality control analysis of RNA sequencing data obtained from organism with an available reference genome and annotation. 

STAGES:
    - preprocessing: 
        - fastp : Read filtering 
        - FASTQC : Read QC
        - SortMeRNA : ribosomal RNA quantification
        - bowtie2 : globin quantification

    - processing:
        - STAR : alignment to the reference genome 
        - salmon: quantification
        - RSeQC : quality control 

    - postprocessing:
        - R : additional qc
        - multiqc : regroup all qc


QUICK START:

1. Prepare input files:
    Metadata Requirements
        The metadata file should include the following mandatory column:
    **SampleID:** should contain the name of the fastq file for each sample
    **SampleName:** refer to the name of the sample as you want it to appear in the figures 
    **Group:** needed to color samples belonging to same group
    **Run_Name:** a new output directory will be created based on run name  
    **Experiment_Name:** an additional new output subdirectory will be created based on experiment name
    **Include:** Specify 'Y' to include the sample in the analysis or 'N' to exclude it.
    **Path:** relative path that specifies the location of the fastq raw file relative to the current working directory  

   Reference Genome and annotation files:
    - Download from your favorite database: gtf, bed transcript files
    - Generate mapping indexes with STAR (STAR   --runMode genomeGenerate   --runThreadN 25   --genomeDir [OUTPUT_DIRECTORY]   --genomeFastaFiles [FASTA_FILE] --sjdbGTFfile [GTF8FILE]   --sjdbOverhang 74)
    - Generate correspondance file between gene id and gene biotype from gtf file  (use get_corresp_geneID_geneBiotype_from_gtf.py provided in script directory)
    - download ribosomal RNA database from https://github.com/sortmerna/sortmerna/tree/master/data/rRNA_databases
    - create index for hemoglobin transcript if your samples are from blood origin

    Config file 
    Prepare a config file similar to the one provided in the test dataset. The config file  mostly contain different paths to the necessary input files. 


2. Create conda environement using provided yaml file 
    > conda env create -n rnaseq-qc -f rnaseq-qc.yaml

    Activate the conda environement
    > conda activate rnaseq-qc

    Run the pipeline with command line: 
    > snakemake --directory test/ --snakefile Snakefile --use-conda --conda-prefix path/to/conda --configfile test/config.json -p 


