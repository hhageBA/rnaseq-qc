import sys

path_to_gtf_file=sys.argv[1]

geneid="gene_id"
biotype="gene_type"  # "gene_type" #"gene_biotype"


with open(path_to_gtf_file, "r") as filein:
	for line in filein:
		if line[0]!="#":
			l=dict()
			list_gtf = line.split('\t')
			if list_gtf[2]=="gene" or list_gtf[2]=="pseudogene":
				l=dict(a.strip('"').split(' "') for a in list_gtf[-1].strip("; \n").split("; ") if '"' in a)
				print(str(l[geneid])+"\t"+str(l[biotype]))

