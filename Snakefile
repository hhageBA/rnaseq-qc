from snakemake.utils import min_version
min_version("5.5.4")
import pprint
include: "rules/common.smk"

# get the path where the main code is stored = var have different name depending of snakemake version
try:
    LAUNCHED_SMK = os.path.dirname(workflow.overwrite_configfile[0])
except:
    LAUNCHED_SMK = os.path.dirname(workflow.overwrite_configfiles[0])


WORK = os.getcwd()
print(WORK)
PATH = os.path.dirname(os.path.abspath(workflow.snakefile))
SAMPLES=read_XLS_metadata(config, WORK)
pprint.pprint(SAMPLES)
RUNS=create_RUNS_dict(SAMPLES, config)
config["runs"]=RUNS
print(RUNS)


rule all:
    input:
        ##### FASTQC
        expand(getfiles(SAMPLES, RUNS, "qc", "{RT}_fastqc.zip", name_type="sampleID"), RT=config["read_type"]),

        ##### Fastp
        expand(getfiles(SAMPLES, RUNS, "preprocessing", "{RT}_preprocessed.fastq.gz"), RT=config["read_type"]),

        ##### count globin proportion
        getfiles(SAMPLES, RUNS, "preprocessing", ".globin.log"),

        ##### Salmon quantification
        getfiles(SAMPLES, RUNS, "salmon", "/quant.genes.sf"),
        getfiles(SAMPLES, RUNS, "salmon", "_eukaryota_count_table.txt", name_type="run"),

        ##### RSeQCs
        getfiles(SAMPLES, RUNS, "mapping", "_read_distrib.txt"),
        getfiles(SAMPLES, RUNS, "mapping", ".inner_distance.txt"),
        getfiles(SAMPLES, RUNS, "mapping", ".infer_exp.txt"),
        ##getfiles(SAMPLES, RUNS, "mapping", "BAMfiles_list_EUK.txt", name_type="noRun2"),
        getfiles(SAMPLES, RUNS, "qc", ".reads.stats.txt"),
        getfiles(SAMPLES, RUNS, "Reporting", "_eukaryota_stats_per_samples.txt", name_type="run"),

        ###### additional QCs
        getfiles(SAMPLES, RUNS, "additionalQC", "Number_of_detected_genes_mqc.png", name_type="noRun2"),

        ###### MultiQC
        getfiles(SAMPLES, RUNS, "Reporting", "_eukaryota_multiqc.html", name_type="run"),


include: "rules/qc/fastqc.smk"
include: "rules/qc/rseqc.smk"
include: "rules/qc/multiqc.smk"
include: "rules/preprocessing/fastp.smk"
include: "rules/preprocessing/sortmerna.smk"
include: "rules/preprocessing/globin_SPE.smk"
include: "rules/mapping/star.smk"
include: "rules/mapping/samtools_index.smk"
include: "rules/quantification/salmon_quant_alignment.smk"
include: "rules/quantification/salmon_combine.smk"
include: "rules/qc/stats_per_sample.smk"
include: "rules/biostat_qc/additional_qc.smk"


