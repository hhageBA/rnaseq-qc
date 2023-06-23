if(!require("BiocManager")) {install.packages("BiocManager")}
library("BiocManager")

##A function of used libaraies in general 
requiredLibrariesFunction <- function(requiredLibraries){
  install.lib <- requiredLibraries[!requiredLibraries %in% installed.packages()]
  for(lib in install.lib) BiocManager::install(lib, ask = FALSE)
  sapply(requiredLibraries,require,character=TRUE)  
} 
update.packages("rlang", "1.1.0")
librariesGeneral <- c("ggplot2","grid","gridExtra", "lattice","reshape", "stringr", "parallel","doParallel","optparse",
                      "plyr","gplots", "RColorBrewer","ggcorrplot","openxlsx","edgeR","phyloseq")


requiredLibrariesFunction(librariesGeneral)

# arguments

option_list <- list(
  make_option(c("--inputData"), type ="character"),
  make_option(c("--metadata"), type ="character"),
  make_option(c("--statFile"), type ="character"),
  make_option(c("--typeFile"), type ="character"),
  make_option(c("--output"), type ="character"),
  make_option(c("--groupColor"), type ="character"),
  make_option(c("--filterThreshold"), type ="integer"), 
  make_option(c("--scripts"), type ="character")

)

opt <- parse_args(OptionParser(option_list=option_list))

inputData <- opt$inputData
metadata <- opt$metadata
statFile <- opt$statFile
typeFile <- opt$typeFile
output <- opt$output
groupColor <- opt$groupColor
filterThreshold <- opt$filterThreshold
pathToScripts <- opt$scripts

sourcingFiles <-function(pathway){
  
  listFiles = list.files(path=pathway,pattern=".R", include.dirs = TRUE, recursive=TRUE)
  for(i in 1:length(listFiles)){
    source(paste0(pathway, listFiles[i]))
  }
  
}

sourcingFiles(paste0(pathToScripts,"/scripts/Stats_scripts/"))

fileName <- ""
outputDataPath <- paste0(dirname(output),"/")
minCount <- filterThreshold

##upload featureSet and metaData 
featureSet <- read.table(inputData,sep="\t",header=TRUE, row.names=1)
featureSet <- round(featureSet, 0)

# Read the metadata file  (if nme sheet ==> sheetName ="")
metaDataAll <- read.xlsx(metadata,skipEmptyCols = TRUE,  sheet = 1) 

metaDataAll[,"SampleName"] <- as.character(metaDataAll[,"SampleName"])

# MetaData$SampleID must correspond to the featureSet columns
metaData <- metaDataAll[which(metaDataAll[,"SampleName"] %in% 
                             colnames(featureSet)),]

#colors for groupColor
ColorPalettes <- function(colorNumber,paletteName){
  
  colors <- hcl.colors(n= colorNumber,palette = paletteName)
  return(colors)
}
colorBy <- ColorPalettes(length(unique(metaData[,groupColor])), "Set2")
names(colorBy) <- unique(metaData[,groupColor])

#order featureSet as metaData
featureSet <- featureSet[,as.character(metaData$SampleName)]

#Name columns featureSet for graph 
colnames(featureSet) <- metaData[,"SampleName"]


stat <- read.table(statFile, sep="\t",header=TRUE)
stat[,"Perc_rRNA"] <- as.numeric(gsub("%","",stat[,"Perc_rRNA"]))
stat[,"Perc_globin"] <- as.numeric(gsub("%","",stat[,"Perc_globin"]))
row.names(stat) <- stat$Sample
stat <- stat[as.character(metaData[,"SampleName"]),]
stat <- cbind(Names = metaData[,"SampleName"],stat)

PlotPercentOfrRNAandGlobin(stat,sampleNames = "Sample",rRNAcol="Perc_rRNA",
                           GlobinCol="Perc_globin",sizex=16,sizey=16, sizeaxistitle=18,
                           GlobinPresent = TRUE,
                           outputPath=outputDataPath)
  

### LOAD ENSG ANNOTATION
gtf = read.table(typeFile, sep ="\t", stringsAsFactors=F)
biotypePerSpl <- getPercBiotypePerSample(featureSet = featureSet,SampleNames =metaData[,"SampleName"] ,gtf =  gtf,NbTypeToPlot = 1:5,
                                         FileName = fileName,
                                          outputPath = outputDataPath)

#Only protein coding genes
mRNA <- as.character(gtf$V1[which(gtf$V2 == "protein_coding")])
featureSet <- featureSet[which(row.names(featureSet) %in% mRNA),]

## Preprocessing 
countsPrePross <- preProcessingTranscriptionFeatureSet(featureSet = featureSet,
                                                       minCount = minCount,
                                                       minSamples = 0,
                                                       priorCountZero = 1,
                                                       normalizeMethod ="RLE",fileName = fileName,
                                                       outputPath = outputDataPath)


countsPreProssFilterNormalizeLog <- countsPrePross$logCpmFilterNormalize


MeanNumberPlotAndSummary <- MeanNumberOffeaturesBySample(featureSet=featureSet , metaData=metaData,
                             colors = colorBy,
                             Group=groupColor,sampleNames ="SampleName",
                             countsThres =10,sizex = 16,sizeaxistitle = 18,sizetext = 6,sizey = 16,
                             outputPath=outputDataPath)


## Rarefaction Curves 

rarefactionPlot <- paste0(outputDataPath,"RarefactionCurvescoloredByHuman_additional_qc.png")
phylo <- PhyloseqObject(countsData= featureSet,
                        metaData=metaData,
                        sampleNames= metaData[,"SampleName"], taxTable = NULL)
depths = rep(c(1, 1:4*10^4, 1:6*10^5 ,1:5* 10^6,1:3*10^7), each = 10)
depths = rep(c(1, 1:4*10^4), each = 2)


phylObject <-otu_table(phylo,taxa_are_rows = TRUE)
countsThres <- 10
RarefractionCurves <- plotRarefactionCurves(phylObject = phylObject, depths = depths,
                                            metaData=metaData, countsThres, 
                                            colorBy =groupColor ,colors = colorBy,
                                            metaDataColMerge = "SampleName",
                                            lineSize = 0.3,
                                            outputPath = outputDataPath)



boxPlotSamplesByGroups(data = countsPreProssFilterNormalizeLog, 
                       sampleNames = metaData[,"SampleName"],
                       groups = metaData[,groupColor],colors = colorBy,angle=90,
                       xlab = "Normalized counts", mainTitle = " Normalized counts ",
                       outputPath = outputDataPath) 

plotDensityDistribution(logRawData = as.matrix(log(featureSet +1)),
                        logFilterNormData = countsPreProssFilterNormalizeLog,
                        metaData =metaData, 
                        colorBy = groupColor ,colors = colorBy,
                        outputPath = outputDataPath)

