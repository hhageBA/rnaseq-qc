
#' @description Mean number of featuresBy sample and by conditions 
#' 
#' @param featureSet:  count matrix with samples on columns 
#' @param metaData: metaData with rows as samples
#' @param Group: character, colnames of metaData used for facet_grid (group BY)
#' @param sampleNames character, colnames of metaData used as labels for samples (x axes)
#' @param countsThres threshold of counts to be considered. By default is equal to 10 
#' @param ylim by default is equal to c(5000,15000)
#' @param ylim: position of additional text in the y-axis
#' @param colors: vector with same length of levels(Group)
#' @param order: the levels of the metaData groups, set to 0 by default
#' @param angle: angle of the text of the x-axis labels 
#' @param angletext: angle of the text added to the plot
#' @param size: size of the text of the axes
#' @param pointSize: size of points of the plots
#' @param sizetext: size of the text added to the plot
#' @param anglefacet: angle of the facet labels
#' @param sizefacet: size of the facet labels
#' @param legendposition: position of the legend, default to "none"
#' @param widthpng: width of the png file to be saved 
#' @param heightpng: height of the png file to be saved 
#' @param sizex: size of the x axis labels
#' @param sizey: size of the y axis labels
#' @param sizetitleaxis: size of the  axes titles
#' @param showxaxislabels: whether to show x axis labels or not (boolean)


MeanNumberOffeaturesBySample <- function(featureSet, metaData,Group,sampleNames, colors, countsThres =10 ,xlim = c(5000,15000),
                                         order = 0, sizey = 15, sizex=15, pointSize =3,sizetext = 3,sizefacet=15,
                                         xtext = 0,angletext=0, showxaxislabels=TRUE,sizeaxistitle=11,
                                         legendposition="none", plotfacet=FALSE,
                                         widthpng=2000,heightpng=800,pointsizepng=12, outputPath)
                                         
{
  # this statement is to make sure to display the mean number of genes in the plot
  if(xtext == 0)
  {
    xtext = min(xlim)+1000
  }
  
  if(length(order)==1)
  {
    order <- unique(metaData[,Group])
  }
  metaData[,Group] = factor(metaData[,Group], levels = order)

  # The following dataframe contains values of mean number of genes for each sample
  dfToPlot=data.frame(ID =metaData[,sampleNames], nbMeanGenes=colSums(featureSet>countsThres),
                      Group=metaData[,Group])
  # the  following function plot the mean number of genes as fonction of the sample names and partitions the plots
  # into k facets defined par the levels of the metaData grouping variable
  MeanNumberPlot=ggplot(dfToPlot,aes(x=nbMeanGenes,y=ID,colour=Group))+geom_point(size = pointSize) + 
    expand_limits(x=xlim)
  if (plotfacet==TRUE)
  MeanNumberPlot <- MeanNumberPlot + facet_grid(.~ Group, scales="free",space="free") + theme(strip.text.y = element_text(size=sizefacet,face="bold"))
 
  MeanNumberPlot <- MeanNumberPlot + scale_color_manual(values=colors) +
    labs(colour = colors, y= "Samples", x= "Number of features expressed")
  if (showxaxislabels==TRUE) # we choose to display the x axis labels if there are not many samples
  MeanNumberPlot <- MeanNumberPlot + theme(axis.title =element_text(size=sizeaxistitle), 
          axis.text.x=element_text(size=sizex),axis.text.y=element_text(size=sizey,face = "bold"),
          legend.text = element_blank(),legend.title = element_blank(),legend.position = legendposition) + 
    geom_text(aes(label = nbMeanGenes,x=xtext),size = sizetext, angle=angletext)
  else 
    MeanNumberPlot <- MeanNumberPlot + theme(axis.title.x =element_text(size=sizeaxistitle), 
                                             axis.text.x=element_text(size=sizex),axis.text.y=element_blank(),axis.title.y=element_blank(),
                                             legend.text = element_blank(),legend.title = element_blank(),legend.position = legendposition) + 
    geom_text(aes(label = nbMeanGenes,x=xtext),size = sizetext, angle=angletext) 

  title <- "Number_of_detected_genes"
 
  #save the summary table as an txt file
  write.table(dfToPlot, file=paste0(outputPath,"SummaryTable",title,".txt"))
  grid.arrange(MeanNumberPlot)
  
  # the plot is saved as a png file
  png(filename = paste0(outputPath,title,"_mqc.png"),
      width = widthpng, height = heightpng, units = "px", pointsize = pointsizepng,
      bg = "white")
  grid.arrange(MeanNumberPlot)
  dev.off()
  return(list(Plot=MeanNumberPlot,Summary=summary(dfToPlot$nbMeanGenes)))
}