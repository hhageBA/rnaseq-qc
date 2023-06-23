rule combine_counts_tables:
	input: expand("Analyses/{analysis}/{run}/salmon/{run}_eukaryota_count_table.txt", analysis=config["analysis"], run=config["runs"].keys())
	output: "Analyses/{analysis}/all/All_eukaryota_runs_count_table.txt"
	params:
		runs=list(config["runs"].keys()),
		#spls_names=[y for x in [config["runs"][o] for o in config["runs"].keys()] for y in x]
	conda:
		"environment.yaml"
	shell:
		"Rscript {PATH}/scripts/combine_counts_tables_runs.R {output} {params.runs} "
