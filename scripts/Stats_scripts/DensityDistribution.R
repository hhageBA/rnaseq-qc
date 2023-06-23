#' @description Density plot of raw, filtered and normalized data 
#' 
#' @param logRawData: log counts data with all genes 
#' @param logFilterData: log counts data filtered (without low counts)
#' @param logFilterNormData: log counts data filtered and normalized 
#' @param metaData: metaData variables with a least the column to be used for colors (colorBy)
#' @param colorBy: chracter (a colnames of metaData) of the variable to be used for colors
#' @param colors: a vector with the colors to use. names(colors) are essensials and should be equal to levels(colorBy)
#' @param outputPath
plotDensityDistribution <- function(logRawData, logFilterNormData, metaData, colorBy, colors, 
                                    ylimRaw = c(0,1.5), ylimFN= c(0,1),main = 0,
                                    outputPath){
  
  if(length(names(colors)) == 0){ # creation du vecteur de couleurs si ce n'est pas deja fait 
    print("names(colors) must be not NULL")
  }else{
    nbcol <- colors
    Groups <- names(colors)
  }
  
  nsamples=dim(logRawData)[2]
  colnames(metaData)[which(colnames(metaData)==colorBy)] <- "Group_to_plot"
  metaData$Group_to_plot <- as.character(metaData$Group_to_plot)

  title <- "DensityPlot"
 
  
  par(mfrow=c(1,2))
  plot(stats::density(logRawData[,1]), col=nbcol[metaData$Group_to_plot[1]], lwd=2, ylim=ylimRaw, las=1,
       main=paste0("A. Raw data (Nb features = ",dim(logRawData)[1],")"), xlab="Log-counts")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- stats::density(logRawData[,i])
    lines(den$x, den$y, col=nbcol[metaData$Group_to_plot[i]], lwd=0.3)
  }
  legend("topright",Groups,col=nbcol,lty=1,lwd=1)
  
  
  plot(stats::density(logFilterNormData[,1]), col=nbcol[metaData$Group_to_plot[1]], lwd=2, ylim=ylimFN, las=1, 
       main=paste0("C. Filtered and normalized data (Nb features = ",dim(logFilterNormData)[1],")"), 
       xlab="Log-counts")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- stats::density(logFilterNormData[,i])
    lines(den$x, den$y, col=nbcol[metaData$Group_to_plot[i]], lwd=0.3)
  }
  legend("topright",Groups,col=nbcol,lty=1,lwd=1)

  pdf(file= paste0(outputPath,title,".pdf"),
      width = 12, height = 5,  pointsize = 12, 
      bg = "white")
  par(mfrow=c(1,2))
  plot(stats::density(logRawData[,1]), col=nbcol[metaData$Group_to_plot[1]], lwd=0.5, ylim=ylimRaw, las=1,
       main=paste0("A. Raw data (Nb features = ",dim(logRawData)[1],")"), xlab="Log-counts")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- stats::density(logRawData[,i])
    lines(den$x, den$y, col=nbcol[metaData$Group_to_plot[i]], lwd=0.3)
  }
  legend("topright",Groups,col=nbcol,lty=1,lwd=1)
  
  
  plot(stats::density(logFilterNormData[,1]), col=nbcol[metaData$Group_to_plot[1]], lwd=0.5, ylim=ylimFN, las=1, 
       main=paste0("C. Filtered and normalized data (Nb features = ",dim(logFilterNormData)[1],")"), xlab="Log-counts")
  abline(v=0, lty=3)
  for (i in 2:nsamples){
    den <- stats::density(logFilterNormData[,i])
    lines(den$x, den$y, col=nbcol[metaData$Group_to_plot[i]], lwd=0.3)
  }
  legend("topright",Groups,col=nbcol,lty=1,lwd=1)
  dev.off()
}
