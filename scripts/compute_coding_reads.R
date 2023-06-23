## Count number of reads mapping on transcriptom and number of coding reads. 

# arguments
args = commandArgs(TRUE)
table_path = args[1]
pre_stats= args[2]
path_type= args[3]
run_nb= args[4]
output_path_tmp <- dirname(table_path)
output_path_final <- dirname(dirname(table_path))

# 1. compute transcriptome and coding reads
t <- read.table(table_path, header=TRUE, row.names=1, check.names=FALSE)
type <- read.table(path_type, header=FALSE, check.names=FALSE)
transcriptome <- colSums(t)
write.table(transcriptome, paste0(output_path_tmp,"/colsums_quantified_reads_all.txt"), sep="\t", quote=FALSE)
coding <- colSums(t[rownames(t)%in%type[type$V2=="protein_coding",]$V1,])
write.table(coding, paste0(output_path_tmp,"/colsums_prot_cod.txt"), sep="\t", quote=FALSE)


# 2. get table
df_T <- as.data.frame(transcriptome)
df_T$Sample <- gsub("X", "", rownames(df_T))
df_C <- as.data.frame(coding)
df_C$Sample <- gsub("X", "", rownames(df_C))
###df_C$Sample <- gsub(".", "-", rownames(df_C))
mrg1 <- merge(df_T, df_C, by="Sample")
colnames(mrg1) <- c("Sample", "Transcriptome_reads", "Coding_reads")

# merge both
pre_stats_tb <- read.table(pre_stats, header=TRUE)
mrg2 <- merge(pre_stats_tb, mrg1, by="Sample")

write.table(mrg2, paste0(output_path_final, "/Reporting/", run_nb, "_eukaryota_stats_per_samples.txt"), quote=FALSE, row.names=FALSE, sep='\t')


