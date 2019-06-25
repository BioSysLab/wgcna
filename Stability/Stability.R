######################Authors : Michael Neidlin          ##################### 
#                               Smaragda Dimitrakopoulou                     #
#                               25-06-2019                                   #
##############################################################################
#THIS SCRIPT PERFOMRS STABILITY ANALYSIS FOR 50x resampling and 10% removal analyses

# INPUT NEEDED
# 1. multiExpr - Expression data (multiExpr.RData)
# 2. signednet - Outcome of WGCNA on original data


#OUTPUT: Dendrogram Plots
#OUTPUT: SignedNet and multiExpr lists for 100 different conditions (50+50)
####################################################################


# Load packages 
library(WGCNA)
library(rlist)


# Load R-Objects
source("WGCNAfunction.R")


##################################################


#########################################################

# Make the random Re-sampling:Rsamples and the 10% out:Rsamples1
nSets=4
NPermutations=50
size<-vector(length=nSets)
limits<-vector(length=nSets)

for (i in 1:4) {

#This is for the removal of Samples  
    size[i]<-round(10/100*length(multiExpr[[i]]$data[,1]))
  
#This is for the Re-Sampling 
    limits[i]<-length(multiExpr[[i]]$data[,1])  
}

#For Bootstap reasmpling you want in the Sample command, replace=TRUE
set.seed(41)

Rsample<-vector(mode="list",length=NPermutations)
Removal<-vector(mode="list",length=NPermutations)
for (j in 1:NPermutations) {
  Rsample[[j]]<-list(data=sample(1:limits[1],limits[1],replace = TRUE),data=sample(1:limits[2],limits[2],replace = TRUE)
                     ,data=sample(1:limits[3],limits[3],replace = TRUE),data=sample(1:limits[4],limits[4],replace = TRUE))
  
  Removal[[j]]<-list(data=sample(1:limits[1],size[1],replace = TRUE),data=sample(1:limits[2],size[2],replace = TRUE)
                     ,data=sample(1:limits[3],size[3],replace = TRUE),data=sample(1:limits[4],size[4],replace = TRUE))
}

# Make the MegaMatrix

Mega_multiExpr<-vector(mode="list",length=NPermutations)
setLabels = c("Meniscus", "Subchondral Bone","Synovium","Cartilage","name")
for (j in 1:(2*NPermutations)) {
#For Removal of Samples you have -Rsample
#For Re-Sampling you have Rsample 
  if (j<=50) {
    
  
  multiExp1<-vector(mode="list",length = nSets)
  multiExp1[[1]]<-list(data=multiExpr[[1]]$data[Rsample[[j]][1]$data,])
  multiExp1[[2]]<-list(data=multiExpr[[2]]$data[Rsample[[j]][2]$data,])
  multiExp1[[3]]<-list(data=multiExpr[[3]]$data[Rsample[[j]][3]$data,])
  multiExp1[[4]]<-list(data=multiExpr[[4]]$data[Rsample[[j]][4]$data,])
  multiExp1[[5]]<-list(name=as.character(j))
  names(multiExp1)<-setLabels
  Mega_multiExpr[[j]]<-multiExp1
  } else {
   
    
    multiExp1<-vector(mode="list",length = nSets)
    multiExp1[[1]]<-list(data=multiExpr[[1]]$data[-Removal[[j-50]][1]$data,])
    multiExp1[[2]]<-list(data=multiExpr[[2]]$data[-Removal[[j-50]][2]$data,])
    multiExp1[[3]]<-list(data=multiExpr[[3]]$data[-Removal[[j-50]][3]$data,])
    multiExp1[[4]]<-list(data=multiExpr[[4]]$data[-Removal[[j-50]][4]$data,])
    multiExp1[[5]]<-list(name=as.character(j))
    names(multiExp1)<-setLabels
    Mega_multiExpr[[j]]<-multiExp1
  }
}
names(Mega_multiExpr)<-rep("multiExpr",NPermutations)


#Run WGCNA
Mega_signednet<-mclapply(Mega_multiExpr,WGCNAfunction,mc.cores = 1) 


#Save Mega_OAmultiExp and Mega_OAsignednet
saveRDS(Mega_multiExpr, file = "Mega_multiExpr.rds")
saveRDS(Mega_signednet, file = "Mega_signednet.rds")


#After Analysis for the plots
## The code for Removal and Resampling  is the same, the index "j" differs

#LOAD signednet from the detection w/o noise (WGCNA function)
moduleLabels = signednet$colors
moduleColors = labels2colors(moduleLabels)
detected_genes<-vector(length = NPermutations)

# Define a matrix of labels for the original and all resampling runs

nGenes<-length(Mega_multiExpr[[1]]$Meniscus$data[1,])
labels = matrix(0, nGenes, NPermutations + 1)
labels[,1] = moduleColors
genes =sum(moduleColors!="grey")


#Finding percentage of detected genes and Resampling plot

for (j in 1:50) {
  Noise_moduleLabels<-Mega_signednet[[j]]$colors
  Noise_moduleColors<-labels2colors(Noise_moduleLabels)
  orig.is.module<-moduleColors!="grey"
  noise.is.module<-Noise_moduleColors!="grey"
  overlap<-orig.is.module+noise.is.module
  
  detected_genes[j]<- (sum(overlap==2)/genes)*100
  labels[,j+1]<-matchLabels(Noise_moduleColors, labels[, 1])
}
Resampling_detected_genes<-detected_genes[1:50]



nBlocks = length(signednet$dendrograms)
pdf(file = "Resampling of Sampes_SamplesdendrogramAndSampledColors.pdf", wi=20, h=25)

for (block in 1:nBlocks) {
  label_string<-labels[signednet$blockGenes[[block]],1]
  for (i in 2:NPermutations+1) {
    label_string<-cbind(label_string,labels[signednet$blockGenes[[block]],i])
  }
  plotDendroAndColors(signednet$dendrograms[[block]],label_string,
                      c("Full data set", paste("Resampling", c(1:50))),
                      main = "Gene dendrogram and module labels from resampled data sets",
                      autoColorHeight = FALSE, colorHeight = 0.9,
                      dendroLabels = FALSE, hang = 0.03, guideHang = 0.05,
                      addGuide = TRUE,
                      guideAll = FALSE,
                      cex.main = 1.2, cex.lab = 1.2, cex.colorLabels = 0.7, marAll = c(0, 0, 3, 0))
  
}
dev.off()




#Finding percentage of detected genes and Removal of 10% plot

labels = matrix(0, nGenes, NPermutations + 1)
labels[,1] = moduleColors

for (j in 51:100) {
  Noise_moduleLabels<-Mega_signednet[[j]]$colors
  Noise_moduleColors<-labels2colors(Noise_moduleLabels)
  orig.is.module<-moduleColors!="grey"
  noise.is.module<-Noise_moduleColors!="grey"
  overlap<-orig.is.module+noise.is.module
  
  detected_genes[j]<- (sum(overlap==2)/genes)*100
  
  labels[,j+1-50]<-matchLabels(Noise_moduleColors, labels[, 1])
}
Removal_detected_genes<-detected_genes[51:100]

 pdf(file = "Removal of Sampes_SamplesdendrogramAndSampledColors.pdf", wi=20, h=15)
 for (block in 1:nBlocks) {
   label_string<-labels[signednet$blockGenes[[block]],1]
   for (i in 2:NPermutations+1) {
     label_string<-cbind(label_string,labels[signednet$blockGenes[[block]],i])
   }
     plotDendroAndColors(signednet$dendrograms[[block]],label_string,
                         c("Full data set", paste("Resampling", c(1:50))),
                         main = "Gene dendrogram and module labels from removal data sets",
                         autoColorHeight = FALSE, colorHeight = 0.9,
                         dendroLabels = FALSE, hang = 0.03, guideHang = 0.05,
                         addGuide = TRUE,
                         guideAll = FALSE,
                         cex.main = 1.2, cex.lab = 1.2, cex.colorLabels = 0.7, marAll = c(0, 0, 3, 0))

}
dev.off()


 

#Boxplot of the percentage of the detected genes with each method
boxplot(cbind(Resampling_detected_genes,Removal_detected_genes),main="Stability",ylab="Percentage %",names=c('Resampling','Removal'))
