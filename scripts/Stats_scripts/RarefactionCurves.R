#' @description a function to compute the rarefying richness for a given depth (number of reads)
#' 
#' 
#' @param  phylObject a phyloseq object class output of the function PhyloseqObject (containing at least counts_table 
#' and sample Data. 
#' @param depth: one value representing the total number of reads  
#' @param countsThres: minimum number of counts for each. equal to 10 by default

estimateRarifiedRichness <- function(phylObject,depth,countsThres = 10, measures = 0){
   
  if(max(sample_sums(phylObject)) < depth) return()
  
  data <- prune_samples(sample_sums(phylObject) >= depth, phylObject)
  rarified_phylObject <- rarefy_even_depth(data, depth, verbose = FALSE,replace=FALSE)
 if(measures != 0)
  {
   alpha_diversity <- estimate_richness(rarified_phylObject, measures = measures)
   RarefiedCounts <- data.frame(Sample=  row.names(alpha_diversity) , alpha_diversity)
   colnames(RarefiedCounts) <- c("Sample", "Measures")
 }
  else
  {
  countsHigherThres <- colSums(rarified_phylObject >= countsThres)
  RarefiedCounts <- data.frame(Sample=  names(countsHigherThres) , countsHigherThres)
  colnames(RarefiedCounts) <- c("Sample", "Measures")
  }
  
  RarefiedCounts
}




#' @description  a function that create rarefaction curves using a phyloseq object, a vector of depths and 
#'  the function estimateRarifiedRichness
#' 
#' @param  phylObject a phyloseq object class output of the function PhyloseqObject (containing at least counts_table 
#' and sample Data. 
#' @param depths: a vector of depth : each value represent the total number of reads
#' @param countsThres: minimum number of counts for each. equal to 10 by default
 

calculateRarefactionCurves <- function(phylObject, depths, countsThres=10,  measures = 0){
  
names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
RarefactionCurvesData <- ldply(depths, estimateRarifiedRichness, phylObject =phylObject ,
                                countsThres = countsThres, measures = measures , .id = 'Depth', 
                                .progress = ifelse(interactive(), 'text', 'none'), .parallel=TRUE,
                               .paropts = list(.packages = c( "phyloseq", "dplyr")))


RarefactionCurvesData$Depth <- as.numeric(levels(RarefactionCurvesData$Depth))[RarefactionCurvesData$Depth]

return(RarefactionCurvesData)
}


############ 
#  rarefaction curves
############ 
#' @description  a function to compute, create and plot rarefaction curves (based on both 
#' estimateRarifiedRichness& calculateRarefactionCurves functions)
#' 
#' @param  phylObject a phyloseq object class output of the function PhyloseqObject (containing at least counts_table 
#' and sample Data. 
#' @param depths: a vector of depth : each value represent the total number of reads
#' @param countsThres: minimum number of counts for each. equal to 10 by default   
#' @param colorBy the variable to be used for colors, 
#' @param ncores number of cores to be used for parallel 
#' @param lineSize size of the curves. 
#' @param legendposition: position of the legend, default to "none"
#' @param sizeaxis: size of the labels and the title of the axes
#' @param sizelegend: size of the labels and the title of the legend

PhyloseqObject <- function(countsData, metaData, sampleNames, taxTable = NULL){
  colnames(countsData) <- sampleNames
  countsData <- otu_table(countsData, taxa_are_rows = TRUE)
  
  row.names(metaData) <- sampleNames
  metaData <- sample_data(as.data.frame(metaData))
  
  if(!is.null(taxTable))
  {taxTable <- tax_table(taxTable)
     tax = TRUE}
  else 
  {tax =FALSE}
  phylObject <- phyloseq(countsData, metaData, taxTable)
  
  
  return(phylObject = phylObject)
}



plotRarefactionCurves <- function(phylObject, depths, metaData, countsThres, measures = 0, legendTitle =  0 ,
                                  colorBy,colors, metaDataColMerge = "SampleID",legendposition="none",sizeaxis=20,
                                  sizelegend=30,
                                  ncores,  lineSize = 0.2,outputPath){
  
  
  
  registerDoParallel(ncores)
  
  RarefactionCurvesData <-  calculateRarefactionCurves(phylObject = phylObject, depths= depths, 
                                                       countsThres=countsThres, measures = measures)
  
  save(RarefactionCurvesData,file=paste0(outputPath,"RarefractionData",".Rdata"))
  
  RarefactionCurvesDataSummary <- ddply(RarefactionCurvesData, c('Depth', 'Sample'), summarise, 
                                        countsHigherThresMean = mean(Measures))
  

  RarefactionCurvesDataSummary=merge(RarefactionCurvesDataSummary,metaData,by.x="Sample",by.y=metaDataColMerge)
  if(max(RarefactionCurvesDataSummary$Depth) > 10^6)
  {RarefactionCurvesDataSummary$Depth=RarefactionCurvesDataSummary$Depth*10^(-6)
  xlab = "Depth (million coding reads)"} 
  else
  {xlab = "Depth (reads)"}
  
  RarefactionCurvesDataSummary[,colorBy] <- as.factor(RarefactionCurvesDataSummary[,colorBy])
  colnames(RarefactionCurvesDataSummary)[which(colnames(RarefactionCurvesDataSummary) == colorBy)] <- "colVar"

  if(measures == 0)
  {
    ylab = paste0("Number of features with counts higher than",countsThres)
  }
  else
  {
    ylab = paste0("alpha diversity: ", measures)
  }
  
  if(legendTitle == 0)
  {
    legendTitle = colorBy
  }
  
  RarefractionCurves = list(RarefactionCurvesData = RarefactionCurvesData, 
                RarefactionCurvesDataSummary = RarefactionCurvesDataSummary)
 
  save(RarefactionCurvesDataSummary,file=paste0(outputPath,"RarefactionCurvesDataSummary",".Rdata"))
  
  graph <-  
    ggplot(data = RarefactionCurvesDataSummary,  mapping = aes(x = Depth,y =countsHigherThresMean,
                                                               colour = colVar,
                                                               group = Sample)) + 
    ylim(range(RarefactionCurvesDataSummary[,"countsHigherThresMean"])) +
    scale_color_manual(limits = as.character(unique(RarefactionCurvesDataSummary[,"colVar"])),values=colors) +
    ylab(ylab) +
    xlab(xlab) +
    theme(strip.text.x = element_text(size=sizeaxis),axis.title =element_text(size=sizeaxis),
          axis.text.y=element_text(size=sizeaxis),axis.text.x=element_text(size=sizeaxis),
          legend.text = element_text(size=sizelegend),legend.position = legendposition,
          legend.title = element_text(size=sizelegend),legend.key.size =  unit(0.40, "in")) + geom_line()+ 
    geom_line(linewidth=lineSize) 
    #labs(color = legendTitle) + 
    #geom_dl(aes(label = Sample), method = list(dl.combine("last.points")), cex = 0.8) 
 
  grid.arrange(graph, ncol=1)
  
  png(filename = paste0(outputPath,"Rarefaction Curves","_mqc.png"),
      width = 1600, height = 953, units = "px", pointsize = 12,
      bg = "white")
 grid.arrange(graph, ncol=1)
  dev.off()

  return(list(RarefactionCurvesData = RarefactionCurvesData, 
              RarefactionCurvesDataSummary = RarefactionCurvesDataSummary))
  
}
