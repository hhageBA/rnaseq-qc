#' @description  A function that takes a counts featureSet with samples as columns, 
#' transform counts to DEG, filter the data with minimum Counts and samples, compute normalization factors and return 
#' log cpm counts and an output logger 

#' @param featureSet: data set of gene expressions, with genes (features) on rows and samples in columns 
#' @param minCount = 10 nombre minimale de reads accept�, par d�fault 10
#' @param minSamples should be set to the minimum % of sample that have more then minCount  
#' @param normalizeMethod = "TMM" default method of normalization 
#' @param priorCountZero number to add to 0 before log 
#' @param fileName for the output name 
#' @param groups = 0 if no groups needed, else list of at least 2 groups (same length as samples (from metaData)) 
                   ##that needed to be tested for minimum count 
preProcessingTranscriptionFeatureSet <- function(featureSet, 
                                                 minCount = 10,
                                                 minSamples=0,
                                                 groups = 0, 
                                                 priorCountZero = 0.5, 
                                                 normalizeMethod = "TMM",fileName,
                                                 outputPath){ 
  
  #create a  DGEList object using a featureSet with gene expression 
    countsNoFilter <- DGEList(counts=featureSet)
  
  #filtering the genes by total and sample count 
    if(is.list(groups) == FALSE)
    {#countsTokeep <- rowSums(featureSet > minCount) >= 5}
      countsTokeep <- which(rowSums(featureSet) > minCount)}
    else  ##if their are two groups, filtring by nb of minimal  samples in both groups combined
    {
      if(length(groups) == 1)
      {
        countsTokeep <- rowSums(countsNoFilter$counts > minCount) >= 
                                  min(table(as.character(groups[[1]])))/2
      }
      
      if(length(groups) == 2)  
     {
       countsTokeep <-rowSums(countsNoFilter$counts > minCount) >= 
         min(table(as.character(groups[[1]]),as.character(groups[[2]])))/2
     }
       
       #countsTokeep <- unique(unlist(lapply(unique(groups[[1]]),
                  #    function(t){lapply(unique(groups[[2]]),
 # function(g){which(rowSums(countsNoFilter$counts[,(groups[[2]]==g)&(groups[[1]]==t)] > minCount)>=
                    #  minSamples*sum((groups[[2]]==g)&(groups[[1]]==t)))})})))}
    }
    
     
  ##Samples to consider in the study 
    countsFilter <- countsNoFilter[countsTokeep, , keep.lib.sizes=FALSE]

  #Normalization
    countsNoFilterNormalize <- calcNormFactors(countsNoFilter,method=normalizeMethod)
    countsFilterNormalize <- calcNormFactors(countsFilter,method=normalizeMethod)
   
    #log counts per million  
    
   CpmFilterNormalize <- cpm(countsFilterNormalize, log=FALSE,prior.count = priorCountZero)
   logCpmNoFilter  <- cpm(countsNoFilter, log=TRUE,prior.count = priorCountZero)
   logCpmNoFilterNormalize <-  cpm(countsNoFilterNormalize, log=TRUE,prior.count = priorCountZero)
   logCpmFilter <- cpm(countsFilter, log=TRUE,prior.count = priorCountZero)
   logCpmFilterNormalize <- cpm(countsFilterNormalize, log=TRUE,prior.count = priorCountZero)
  
  

  return(list(countFeatureSetFilterNormalize = countsFilterNormalize,
              countsNoFilter = countsNoFilter, countsNoFilterNormalize = countsNoFilterNormalize,
              countsFilter = countsFilter, countsFilterNormalize = CpmFilterNormalize,
              logCpmNoFilter = logCpmNoFilter, logCpmNoFilterNormalize = logCpmNoFilterNormalize,
              logCpmFilter = logCpmFilter, logCpmFilterNormalize = logCpmFilterNormalize))
}

