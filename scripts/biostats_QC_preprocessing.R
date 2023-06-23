######### Arguments #########
args = commandArgs(TRUE)
biotracs_path=args[1]
inputDataPath=args[2]
outputDataPath=args[3]
run=args[4]
org=args[5]
metadata_path =args[6] 
corresp_geneID=args[7]
groupChosen=args[8]


#########  Packages and functions #############
## Loading Dependencies
pathAutoload= paste0(biotracs_path, "/biotracs-r/tests")
source(paste0(pathAutoload, "/autoloadLib.R"))

pathBiotracs= paste0(biotracs_path)
autoloadLib(dirBiotracs = pathBiotracs, dependencies = c('qc','metagenomics','atlas',"transcription"))


## input data
if (groupChosen==""){
    groupChosen="SampleID"
}
dataName <- paste0(run, "_", org, "_count_table")


# upload dataset 
featureSet <- read.table(paste0(inputDataPath,dataName,".txt"),sep="\t",header=TRUE, row.names=1)
featureSet <- RemoveAddedX(featureSet)
colnames(featureSet) <- gsub("_S.*", "", colnames(featureSet))
dim(featureSet)
featureSet <- round(featureSet, 0)
colnames(featureSet)

# upload metadata 
metaData <- read.table(paste0(metadata_path),sep = ";",header=TRUE)
metaData <- metaData[which(as.character(metaData$SampleID) %in% colnames(featureSet)),]

featureSet <- featureSet[,as.character(metaData$SampleID)]
dim(featureSet)

# get gene annotation 
gtf = read.table(paste0(corresp_geneID), sep ="\t", stringsAsFactors=F)
row.names(gtf) <- gtf$V1
dim(gtf)

biotypePerSpl <- getPercBiotypePerSample(featureSet, gtf)
biotypePerSpl$Gene_Type <- rownames(biotypePerSpl)
biotypePerSpl <- biotypePerSpl[order(biotypePerSpl[,2], decreasing = TRUE),]
biotypePerSpl2 <- biotypePerSpl[seq(1,10),]

mt <- melt(biotypePerSpl2,id.vars = "Gene_Type")

plotBiotype <- ggplot(mt, aes(x=variable, y=value, fill=Gene_Type))+geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 60))
ggsave(plotBiotype, file=paste0(outputDataPath, "geneBiotypeDistribution.png"))

write.csv(biotypePerSpl, file=paste0(outputDataPath, "ReadpercentagesPrurigoRun1.csv"), sep=";")


mRNA <- gtf$V1[which(gtf$V2 == "protein_coding")]
featureSet <- featureSet[which(row.names(featureSet) %in% mRNA),]
dim(featureSet)



###PreProcess 

ro <- rowSums(featureSet)
featureSetFilter <- featureSet[which(ro > 10),]
dim(featureSetFilter)


#Normalization 
#metaData[,which(colnames(metaData)==groupDESeq)] <- as.character(metaData[,which(colnames(metaData)==groupDESeq)])

countsPrePross <- preProcessingTranscriptionFeatureSet(featureSet = featureSet, minCount = 10,minSamples = 0,
                                                       groups = 0,
                                                       priorCountZero = 1,
                                                       normalizeMethod ="RLE",outputPath = outputDataPath)


countsPreProssFilterLog <- log(countsPrePross$countsFilter$counts+1)

countsPreProssFilterNormalizeLog <- countsPrePross$logCpmFilterNormalize
dim(countsPreProssFilterNormalizeLog)



colorByMouse <- ColorPalettes(ncol(featureSet), "Set2")
names(colorByMouse) <- unique(metaData[,which(colnames(metaData)==groupChosen)])

#1- boxplot des expresstion 
metaData[,which(colnames(metaData)==groupChosen)] <- as.character(metaData[,which(colnames(metaData)==groupChosen)])
boxPlotSamplesByGroups(data = countsPreProssFilterNormalizeLog, 
                       sampleNames = metaData$SampleID,
                       groups = metaData[,which(colnames(metaData)==groupChosen)],colors = colorByMouse,
                       ylab = "counts", mainTitle = "Counts RLE",
                       outputPath = outputDataPath)



plotDensityDistribution(logRawData = log(featureSet +1),
                        logFilterNormData = countsPreProssFilterNormalizeLog,
                        metaData =metaData , ylimFN = c(0,0.3),
                        colorBy = groupChosen ,colors = colorByMouse,
                        outputPath = outputDataPath)

colorBySample <- c()
PCAModelRLE <- ordinationModelClinicalData(featureSetSampleRows=
                                             t(countsPreProssFilterNormalizeLog),
                                           metaData=metaData,
                                           model = "PCA",
                                           ntop=0,
                                           clinicalVariables =colnames(metaData),
                                           distance = "bray", 
                                           pc=3, 
                                           scale=TRUE, 
                                           colors = colorBySample,
                                           sampleNames = metaData$SampleID,
                                           normalizartion = "RLE", paletteName= "Set 2",
                                           ellipse =FALSE,
                                           groupeShape =0,
                                           nameGroupShape =0,
                                           outputPath= outputDataPath)



metaData$SampleName <- paste0("Samp",metaData$SampleID)
MeanNumberOffeaturesBySample(featureSet=featureSet , metaData=metaData,
                             colors = colorByMouse,
                             Group=groupChosen,sampleNames = "SampleName", countsThres =10,ylim = c(5000,12000),
                             main = "",
                             outputPath=outputDataPath)
nbGenesZymoMouseallRunq=colSums(featureSet>10)




####rarefaction curves
#colnames(metaData)[1] <- "SampleID"
phylo <- PhyloseqObject(countsData= featureSet,
                        metaData=metaData,
                        sampleNames= metaData$SampleID, taxTable = NULL)
#courbe de rerafraction 
depths = rep(c(1, 1:4*10^4, 1:9*10^5 ,1:10* 10^6,1:3*10^7), each = 10)
phylObject <-otu_table(phylo,taxa_are_rows = TRUE)
countsThres <- 10
metaData[,which(colnames(metaData)==groupChosen)] <- as.character(metaData[,which(colnames(metaData)==groupChosen)])
RarefractionCurves <- plotRarefactionCurves(phylObject = phylObject, depths = depths, 
                                                         metaData=metaData, countsThres, 
                                                         colorBy = groupChosen,colors = colorByMouse,
                                                         ncores = 15, main = "", 
                                                         lineSize = 0.2,outputPath = outputDataPath)
save(RarefractionCurves,file=paste0(outputDataPath,"RarefractionCurves.Rdata"))









