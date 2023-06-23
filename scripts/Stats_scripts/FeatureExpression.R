 

#' @param data: data frame with sample in columns and features in rows 
#' @param groups: a vector of length number of samples ( = dim(data)[2]) giving the group of each sample 
#' @param sampleNames: vector with names attributed to each colomnes of data (name shonw on yhe x-axis)
#' @param colors:  a vector of length groups giving a color for each group 
#' @param ylab: names on the y axis
#' @param sizex: size of the x axis labels
#' @param sizey: size of the y axis labels
#' @param sizetitleaxis: size of the  axes titles
#' @param widthpng: width of the png file to be saved 
#' @param heightpng: height of the png file to be saved

boxPlotSamplesByGroups <- function(data, groups, sampleNames,colors, mainTitle, xlab,outputPath="0", sizey = 15,
                                   legTitle= 0,angle = 60,sizex=20,widthpng=2000,heightpng=800,sizetitleaxis=20,
                                   pointsizepng=12,sizelegend=20,hjust = 0.5){
  gr <- rep(groups, each = dim(data)[1])
  data.plot <- melt(data)
  if (ncol(data.plot)==3)
    colnames(data.plot)[2] <- "Sample"
  else 
    colnames(data.plot)[1]<-"Sample"
  data.plot <- cbind(data.plot, group = gr)
  data.plot$Sample <- as.character(data.plot$Sample)
  
  # if no legend title provide, choose "Groups"
  if(legTitle == 0)
  {
    legTitle = "Groups"
  }
  graph <- ggplot(aes(x=value, y= Sample, fill = group), data=data.plot)+
    geom_boxplot()+
    scale_y_discrete(limits = unique(data.plot$Sample),labels = sampleNames)+
    scale_fill_manual(limits = levels(data.plot$group),values=colors)+
    theme(axis.text.x= element_text(size= sizex),
          #axis.text.y = element_text(fsize=sizey),
          axis.line = element_line(color = "black", linewidth = 1,linetype = "solid"),
          legend.text = element_text(size=sizelegend, 
                                     face="bold"),
          legend.title = element_text(size=sizelegend, 
                                      face="bold"),
          legend.position = "bottom",
          axis.title.x = element_text(size=sizetitleaxis),
          plot.title = element_text(face="bold",hjust = hjust)) + 
    ylab(NULL)+
    xlab(xlab)+
    labs(fill = legTitle) 
    ggtitle(mainTitle)
    
  
    
  if(outputPath != "0")
  {
    png(filename = paste0(outputPath,mainTitle, "boxplot_mqc.png"),
        width = 1700, height = 900, units = "px", pointsize = 12,
        bg = "white")
    gr <- grid.arrange(graph,ncol=1,nrow = 1)
  dev.off()}
  else
  {
    grid.arrange(graph,ncol=1,nrow = 1)
  }

  return(boxPlot = graph) 
}





