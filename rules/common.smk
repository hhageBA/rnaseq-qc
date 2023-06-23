from snakemake.utils import validate
import re
import os
import os.path
import glob
from subprocess import Popen, PIPE
import pandas as pd
import pprint


suffix_spl = {"Nextseq2000": "_S*", "Nextseq500": "_S*",  "Novaseq": "", "SP_NextSeq":"", "Miseq": "_S*"}
def read_XLS_metadata(config, WORK): 
    
    try:
        if (os.path.splitext(config["metadata"])[1]==".xlsx" or os.path.splitext(config["metadata"])[1]==".xlsm"):
            df = pd.read_excel(config["metadata"], config["analysis_type"])  # read sheet
        elif os.path.splitext(config["metadata"])[1]==".csv":
            df = pd.read_csv(config["metadata"], sep=";", encoding='latin-1')
    except:
        print("Field 'metadata' of the config file is missing or file could not be read. Exiting.")
        raise
    df=df.astype({'SampleID': 'str'})
    d = df.to_dict('index')
    SAMPLES = {}
    for k in d:
        samp_id = str(d[k]["SampleID"]) # name of the sample in raw file demultiplexed (in the sample sheet)
        samp_name= str(d[k]["SampleName"]) # informative name of the sample

        #check for special characters in sample names
        special_characters = "!@#$%^&*()+?=,<>/ "
        if any(c in special_characters for c in samp_id):
            sys.exit(f"ERROR: Whitespace and/or special character detected in sample ID \"{samp_id}\". Please correct the metadata file. Exiting.")
        if any(c in special_characters for c in samp_name):
            sys.exit(f"ERROR: Whitespace and/or special character detected in Sample name \"{samp_name}\". Please correct the metadata file. Exiting.")

        if d[k]["Include"] == "Y" :
            datapath = os.path.join(WORK, "Data", d[k]["Path"])

            try:
                try:
                    fastq_R1 = glob.glob(os.path.join(datapath, f"{samp_id}_*R1*.fastq.gz"))[0]
                    fastq_R2 = glob.glob(os.path.join(datapath, f"{samp_id}_*R2*.fastq.gz"))[0]
                except:
                    fastq_R1 = glob.glob(os.path.join(datapath, f"{samp_name}_*R1*.fastq.gz"))[0]
                    fastq_R2 = glob.glob(os.path.join(datapath, f"{samp_name}_*R2*.fastq.gz"))[0]
            except:
                sys.exit(f"ERROR: Sample raw data fastq not found for sample \"{samp_name}\" in directory \"{datapath}\". Please check what happend. Exiting.")
            # get demultiplexed names
            suffix_to_replace=config["read_type"][0]+".fastq.gz"
            #print(suffix_to_replace)
            demultiplexed_name = os.path.basename(fastq_R1).replace(suffix_to_replace,"")
            # create sample dict
            SAMPLES[samp_name] = {"datapath":datapath, "R1":fastq_R1, "R2":fastq_R2, "sampleID":samp_id,
                                            "demultiplexed":demultiplexed_name,
                                            "Path": d[k]["Path"], 
                                            "analysis": os.path.dirname(d[k]["Path"]),
                                            "RunID": d[k]["Run_Name"]}
    return(SAMPLES)


# function to create RUNs dictionary
suffix_spl = {"Nextseq500": "_S*",  "Novaseq": "", "SP_NextSeq":"", "Miseq": "_S*"}
def create_RUNS_dict(SAMPLES, config):
    RUNS={}
    for spl in SAMPLES.keys():
        if  SAMPLES[spl]['RunID'] not in RUNS.keys():
            RUNS[SAMPLES[spl]['RunID']] = {"samples": [spl], "Path":SAMPLES[spl]['Path'], "analysis": SAMPLES[spl]['analysis']}
        else:
            # add sample
            RUNS[SAMPLES[spl]['RunID']]['samples'].append(spl)
    return(RUNS)


def getfiles(SAMPLES, RUNS, method, suffix, name_type="sample"):
    listfiles = []
    if "sample" in name_type:
        for spl in SAMPLES.keys():
            if name_type=="sample":
                spl_id = spl
            if name_type=="sampleID":
                spl_id = SAMPLES[spl]["demultiplexed"]
            analysis = SAMPLES[spl]["Path"]
            #if name_type=="run":
            #    spl_id = spl
            listfiles.append(
                f"Analyses/{analysis}/{method}/{spl_id}{suffix}")
    if "run" in name_type:
        for run in RUNS.keys():
            analysis=RUNS[run]["Path"]
            listfiles.append(
                    f"Analyses/{analysis}/{method}/{run}{suffix}")
    if "noRun2" in name_type:
        for run in RUNS.keys():
            analysis=RUNS[run]["Path"]
            listfiles.append(
                    f"Analyses/{analysis}/{method}/{suffix}")
    return listfiles


