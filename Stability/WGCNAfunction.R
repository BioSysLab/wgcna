WGCNAfunction<-function(multiExpr) { 
  
  #Remove the name of the multiExpr
  
  library(rlist)
  Name<-multiExpr[[5]]$name
  multiExpr<-multiExpr[-5]  
  
 
  signednet = blockwiseConsensusModules(multiExpr, power = 20,
                                        minModuleSize = 30,maxBlockSize = 5000,
                                        networkCalibration = "single quantile",
                                        corType = "bicor",maxPOutliers = 0.98,
                                        deepSplit = 2,networkType = "signed",
                                        pamRespectsDendro = TRUE, pamStage = TRUE,
                                        mergeCutHeight = 0.2, numericLabels = TRUE,
                                        minKMEtoStay = 0.2,
                                        saveTOMs = TRUE, verbose = 5)
  list.append(signednet,name=Name)
  

    return(signednet)
}