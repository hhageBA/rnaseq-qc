#' @description  this function is used to get a table of percentage of gene type per sample for RNAseq data

#' @param featureSet: count table with samples as columns and genes as row (no column for gene name) 
#' @param gtf: correspondance table of gene id with gene type (2 columns) with no header names
#' @param NbTypeToplot:  the number of genetypes to be plotted
#' @param SampleNames character, colnames of metaData used as labels for samples (x axes)
#' @param GroupToFacet: character, colnames of metaData used for facet_grid (group BY)
#' @param angle: angle of the text of the x-axis labels 
#' @param widthpng: width of the png file to be saved 
#' @param heightpng: height of the png file to be saved
#' @param sizetextx: size of the x axis labels
#' @param sizetexty: size of the y axis labels
#' @param sizetitleaxis: size of the title of the axes
#' @param sizetitle: size of the main title
#' @param sizelegend: size of the legend elements
#' @param showxaxislabels: whether to show x axis labels or not (boolean)


getPercBiotypePerSample <- function(featureSet, gtf, NbTypeToPlot = 1:10,SampleNames,FileName, outputPath,size=10, groupToFacet = NULL,
                                    angle=90,widthpng=2000,heightpng=800,sizetitleaxis=15,sizetexty=15,sizetextx=20,
                                    pointsizepng=12,sizelegend=15,sizetitle=20, showxaxislabels=TRUE,legend.position = "bottom"
                                    ){
  percReadsBiotypes <- data.frame(matrix(NA, nrow = length(unique(gtf$V2)), ncol = ncol(featureSet)+1))
  colnames(percReadsBiotypes) <- c(SampleNames, "Mean")
  rownames(percReadsBiotypes) <- c(as.vector(unique(gtf$V2)))
  
  for (type in unique(gtf$V2)){
    genes_selected<- gtf[gtf$V2==type,]$V1
    selectedPerc = colSums(featureSet[rownames(featureSet)%in%genes_selected,]) / colSums(featureSet)
    mean_val <- mean(selectedPerc)
    percReadsBiotypes[type,] <- as.numeric(c(selectedPerc, mean_val))
  }
  percReadsBiotypes <- percReadsBiotypes[order(percReadsBiotypes[,2], decreasing = TRUE),]
  
  
  #Plot and csv
  percReadsBiotypes$Gene_Type <- rownames(percReadsBiotypes)
  percReadsBiotypes <- percReadsBiotypes[order(percReadsBiotypes[,2], decreasing = TRUE),]
  
  percReadsBiotypesPlot <- percReadsBiotypes[NbTypeToPlot, -which(colnames(percReadsBiotypes) == "Mean")]
  co <- colSums(percReadsBiotypesPlot[,1:ncol(featureSet)])
  percReadsBiotypesPlot <- rbind(percReadsBiotypesPlot, Others = c(1- co, "Others"))
  
  
  
  mt <- melt(percReadsBiotypesPlot,id.vars = "Gene_Type")
  mt$value <- as.numeric(as.character(mt$value))
  mt$Gene_Type <- factor(mt$Gene_Type, levels = row.names(percReadsBiotypesPlot))
  

  plotBiotype <- ggplot(mt, aes(y=variable, x=value, fill=Gene_Type))+geom_bar(stat="identity") 
  if (showxaxislabels==TRUE)
  plotBiotype <- plotBiotype + theme( axis.text.x = element_text(size = sizetextx), axis.text.y = element_text(size = sizetexty),
           axis.title.x = element_text(size = sizetitleaxis),axis.title.y = element_text(size = sizetitleaxis),
           legend.title=element_text(size=sizelegend), legend.text=element_text(size=sizelegend),legend.position = legend.position,
           plot.title = element_text(size = sizetitle)) + 
    ylab("Samples")+ xlab("Gene Type Frequency")
  else 
    plotBiotype <- plotBiotype + theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size = sizetextx),
                                        axis.title.y = element_blank(),axis.title.x = element_text(size = sizetitleaxis),
                                        legend.title=element_text(size=sizelegend), legend.text=element_text(size=sizelegend),legend.position = legend.position,
                                        plot.title = element_text(size = sizetitle)) + 
                                        ylab("Gene Type Frequency")

  if(length(groupToFacet) > 0)
  {
    mt$Group <- rep(groupToFacet,each = max(NbTypeToPlot)+1)
    plotBiotype <- plotBiotype + facet_grid(. ~Group,scales='free_x', space ="free_x") + 
      theme(strip.text.y = element_text(size = 15, angle = 0))
  }
  png(filename = paste0(outputPath, "Gene biotype distribution",FileName,"_mqc.png"),
      width = widthpng, height = heightpng, units = "px", pointsize = pointsizepng,
      bg = "white")
  grid.arrange(plotBiotype, ncol=1)
  dev.off()
  
  # save the dataframe to a csv file
  write.csv2(percReadsBiotypes,
            file=paste0(outputPath,
                        "Readpercentages",FileName,".csv"))
  
  return(list(Table = percReadsBiotypes,Plot = plotBiotype))
}


