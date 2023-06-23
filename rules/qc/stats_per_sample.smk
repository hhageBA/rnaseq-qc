headers_of_sequencers = { "Nextseq2000": "@VH00881", "Nextseq500": "@NB501625",  "Novaseq": "@A00729", "SP_NextSeq":"@NS500165"}

rule stats_per_sample: 
	input:
		fastp_json = "Analyses/{analysis}/{run}/preprocessing/{sample}_fastp.json",
		mapped = "Analyses/{analysis}/{run}/mapping/{sample}.Log.final.out",
		rrna = "Analyses/{analysis}/{run}/preprocessing/{sample}.rrna.out.log",
		globin= "Analyses/{analysis}/{run}/preprocessing/{sample}.globin.log"
	output: "Analyses/{analysis}/{run}/qc/{sample}.reads.stats.txt"
	benchmark: "Analyses/{analysis}/{run}/benchmarks/stats_sample/{sample}.txt"
	params:
		to_grep=headers_of_sequencers[config["sequencer"]], # @NB501625 for NextSeq #
	shell:
		"""
		nb_raw=`grep '"before_filtering":' -A 1 {input.fastp_json} | grep -oP "(?<=\\"total_reads\\":).*(?=,)"` ; \
		nb_cleaned=`grep '"after_filtering":' -A 1 {input.fastp_json} | grep -oP "(?<=\\"total_reads\\":).*(?=,)"` ; \
		nb_unique_mapped=`grep "Uniquely mapped reads number" {input.mapped} |  cut -d "|" -f2 | sed "s/ //g" | sed "s/\t//g" | sed s/%//g` ; \
		perc_unique_mapped=`grep "Unique" {input.mapped} | grep % | cut -d "|" -f2 | sed "s/ //g" | sed "s/\t//g" | sed s/%//g` ; \
		nb_multi_mapped=`grep "Number of reads mapped to multiple loci" {input.mapped} |  cut -d "|" -f2 | sed "s/ //g" | sed "s/\t//g" | sed s/%//g` ; \
		perc_multi_mapped=`grep "% of reads mapped to multiple loci" {input.mapped} |  cut -d "|" -f2 | sed "s/ //g" | sed "s/\t//g" | sed s/%//g` ; \
		perc_rrna=`grep "passing E-value threshold" {input.rrna} | cut -d"(" -f2 | cut -d")" -f1 ` ; \
		perc_glob=`grep "overall" {input.globin} | grep -oP '[0-9\.]*(?=%)'` ; \
		echo -e "{wildcards.sample}\t${{nb_raw}}\t${{nb_cleaned}}\t${{nb_unique_mapped}}\t${{perc_unique_mapped}}\t${{nb_multi_mapped}}\t${{perc_multi_mapped}}\t${{perc_rrna}}\t${{perc_glob}}" > {output}
		"""


rule compute_coding_reads:
	input:
		count_tb="Analyses/{analysis}/{run}/salmon/{run}_eukaryota_count_table.txt",
		stats="Analyses/{analysis}/{run}/qc/all_eukaryota_stats_per_samples.txt"
	output:
		"Analyses/{analysis}/{run}/Reporting/{run}_eukaryota_stats_per_samples.txt"
	params:
		corresp_gtf=config["gene_type_corresp"],
		run_nb=lambda w: w.run
	log:
		"Analyses/{analysis}/{run}/logs/stats/{run}.eukaryota_codingReads.log"
	conda:
		"environment.yaml"
	shell:
		"Rscript {PATH}/scripts/compute_coding_reads.R {input.count_tb} {input.stats}  {params.corresp_gtf} {params.run_nb}"




rule concat_stats: 
	input: lambda w: expand("Analyses/{{analysis}}/{{run}}/qc/{sample}.reads.stats.txt", sample=RUNS[w.run]["samples"])
	output: "Analyses/{analysis}/{run}/qc/all_eukaryota_stats_per_samples.txt"
	shell:
		"""
		echo -e "Sample\tRaw_reads\tCleaned_reads\tUniquely_mapped\tPerc_uniquely_mapped\tMulti_mapped\tPerc_multi_mapped\tPerc_rRNA\tPerc_globin" > {output} ; \
		cat {input} >> {output}
		"""
