sourcingFiles <-function(pathway){
  
  listFiles = list.files(path=pathway,pattern=".R", include.dirs = TRUE, recursive=TRUE)
  for(i in 1:length(listFiles)){
    source(paste0(pathway, listFiles[i]))
  }
  
}


autoloadLib <- function(dirBiotracs,dependencies){
  corePath = paste0(dirBiotracs, "/biotracs-r/","core/")
  source(paste0(dirBiotracs, "/biotracs-r/","core/","installLibrariesFunction.R"))
  sourcingFiles(corePath)
  
  for (i in 1:length(dependencies))
  {
    dependenciesPath= paste0(dirBiotracs, "/biotracs-r-", dependencies[i], "/", dependencies[i], "/" )
    sourcingFiles(dependenciesPath)
  }
  
}

