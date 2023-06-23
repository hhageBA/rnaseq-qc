#' @description Mean number of features By sample and by conditions 
#' 
#' @param stat: matrix of stats for all samples 
#' @param sampleNames character, colnames of metaData used as labels for samples (x axes)
#' @param ylim: position of additional text in the y-axis
#' @param main if needed, else is equal to 0 
#' @param sizex: size of x axis labels
#' @param sizey: size of y axis labels
#' @param pointSize: size of points of the plots
#' @param sizetitleaxis: size of the title of the axes
#' @param sizetitle: size of the main title
#' @param sizelegend: size of the legend elements
#' @param legendposition: position of the legend, default to "top"
#' @param rRNAcol: column in the stat table containing the rRNA percentage
#' @param GlobinCol: column in the stat table containing the globin percentage
#' @param GlobinPresent: TRUE if the globin percentage is present
#' @param widthpng: width of the png file to be saved 
#' @param heightpng: height of the png file to be saved 

PlotPercentOfrRNAandGlobin <- function(stat,sampleNames, xlim = c(0,80),rRNAcol="Perc_rRNA",GlobinCol="Perc_globin",
                                          sizey = 15, sizex=15, pointSize =3,main=" percentage_mqc.png",
                                          sizeaxistitle=10,position.dodge=0.8,sizelegend=15,
                                         legend.position="bottom", sizetitle = 20, GlobinPresent=TRUE,
                                         widthpng=2000,heightpng=800,pointsizepng=12, outputPath, angle=90)
  
{ 
  
  if(GlobinPresent==FALSE)
  {
    GlobinCol <- NULL
    Percent <- as.numeric(c(stat[,rRNAcol]))
    Type <- rep("rRNA",nrow(stat))
    Samples <- as.character(stat[,sampleNames])
  }
  else {
    Percent <- as.numeric(c(stat[,rRNAcol],stat[,GlobinCol]))
    Type <- rep(c("rRNA","Globin"), each= nrow(stat))
    Samples <- rep(as.character(stat[,sampleNames]), 2)
  }
  
  # The following dataframe contains values 
  df <- as.data.frame(cbind(Samples = Samples,Percent = Percent, Type = Type))
  df$Percent <- as.numeric(as.character(df$Percent))
 
  p <- ggplot(df, aes(x = Percent , y = Samples)) +
    scale_y_discrete(limits=stat[,sampleNames])  +
    geom_col(aes(color = Type, fill = Type), position = position_dodge(position.dodge), width = 0.7) +
    scale_color_manual(values = c("rRNA" = "#00BFC4", "Globin" = "#CC79A7")) +
    scale_fill_manual(values = c("rRNA" = "#00BFC4", "Globin" = "#CC79A7"))+ 
    theme( axis.text.x = element_text(size = sizex), axis.text.y = element_text(size = sizey),
           axis.title.x = element_text(size = sizeaxistitle),axis.title.y = element_text(size = sizeaxistitle),
           legend.title=element_text(size=sizelegend), legend.text=element_text(size=sizelegend),legend.position = legend.position,
           plot.title = element_text(size = sizetitle)) + expand_limits(x=xlim)
if(GlobinPresent==TRUE)
  main=paste0(paste(c(unique(df$Type)),collapse =" and "),main)
else 
    main=paste0(unique(df$Type),main)
  png(filename = paste0(outputDataPath,main),
      width = widthpng, height = heightpng, units = "px", pointsize = pointsizepng,
      bg = "white")
  grid.arrange(p, ncol=1)
  dev.off()
  p 
  

}

    
