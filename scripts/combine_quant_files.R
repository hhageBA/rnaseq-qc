#!/usr/bin/env Rscript
library(tximport)
library(readr)

# Recup des arguments
args = commandArgs(TRUE)
nbr_echs = args[1]
output=args[2]
end=2+as.numeric(nbr_echs)
spls_names=args[3:end]
paths=args[as.numeric(end+1):as.numeric(length(args))]

# Lecture des tables de quantif
full_table = tximport(paths, type="salmon", txIn=FALSE, dropInfReps=TRUE, geneIdCol="Name")


# pour les noms des echantillons
name_splitted <- strsplit(spls_names, "_")
recup_names <- c()
for (i in 1:length(name_splitted)){
  recup_names <- c(recup_names, name_splitted[[i]][1])
}

# creation de la table
quantification_table <- as.data.frame(full_table$counts, rownames=NULL)
#dim(quantification_table)
#length(recup_names)

colnames(quantification_table) <- c(spls_names)

write.table(quantification_table, file=paste(output), quote=FALSE, sep="\t")
