# Rscript to combine count tables from different runs in one project 

args = commandArgs(TRUE)
output=args[1]
runs=args[2:length(args)]


prefix_path <- paste0(dirname(dirname(output)),"/")
for (i in seq(1:length(runs))){
	run=runs[i]
	# premier run
	if (i==1){
		combined_tb <- read.table(paste0(prefix_path, run, "/salmon/", run, "_eukaryota_count_table.txt"), header=TRUE, stringsAsFactor=FALSE)
		combined_tb$Gene <- rownames(combined_tb)
	}else{ # autres runs 
		to_add_tb <- read.table(paste0(prefix_path, run, "/salmon/", run, "_eukaryota_count_table.txt"), header=TRUE, stringsAsFactor=FALSE)
		to_add_tb$Gene <- rownames(to_add_tb)
		combined_tb <- merge(combined_tb, to_add_tb, by="Gene")
	}
}

write.table(combined_tb,file=output, quote=FALSE, row.names=FALSE, sep='\t')
