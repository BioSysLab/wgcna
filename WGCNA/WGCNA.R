######################Authors : Michael Neidlin          ##################### 
#                               Smaragda Dimitrakopoulou                     #
#                               11-06-2019                                   #
##############################################################################

#THIS SCRIPT PERFORMS ALL STEPS OF THE WGCNA ANALYSIS AS PRESENTED IN THE MANUSCRIPT

#INPUT NEEDED: multiExpr variable taken from the preprocessing folder (multiExpr.RData)

#OUTPUT: GRAPHS AS PRESENTED IN THE MANUSCRIPT
#OUTPUT: DISEASE SIGNATURE AS A LIST OF GENES (DiseaseSig.txt)
#OUTPUT: GENE NAMES OF ALL META MODULES (MMi.txt)

####################################################################

library(WGCNA)
library(tidyverse)
library(gplots)


##################################
##### VARIABLE DECLARATION #######
##################################

#Set Labels
nSets=4
setLabels = c("Meniscus", "Subchondral Bone","Synovium","Cartilage")


#Set Trait matrix
Traits<-vector(mode = "list", length = nSets)
Traits[[1]]<-list(data=as.data.frame(c(rep(0,10),rep(1,12))))
Traits[[2]]<-list(data=as.data.frame(c(rep(0,5),rep(1,20))))
Traits[[3]]<-list(data=as.data.frame(c(rep(0,10),rep(1,10))))
Traits[[4]]<-list(data=as.data.frame(c(rep(0,11),rep(1,10))))

for (set in 1:nSets){
 
  names(Traits[[set]]$data)<- "Disease"
}
names(Traits)<-setLabels


##################################
##### VARIABLE DECLARATION END ###
##################################



#################################
###NETWORK CONSTRUCTION (WGCNA) #
#################################

signednet = blockwiseConsensusModules(multiExpr, power = 20,
                                      minModuleSize = 30,maxBlockSize = 5000,
                                      networkCalibration = "single quantile",
                                      corType = "bicor",maxPOutliers = 0.98,
                                      deepSplit = 2,networkType = "signed",
                                      pamRespectsDendro = TRUE, pamStage = TRUE,
                                      mergeCutHeight = 0.2, numericLabels = TRUE,
                                      minKMEtoStay = 0.2,
                                      saveTOMs = TRUE, verbose = 5)
consMEs = signednet$multiMEs

moduleLabels = signednet$colors
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)

table(moduleColors)


#############################################
#######NETWORK CONSTRUCTION (WGCNA) END######
#############################################



#################################
##########DENDROGRAM PLOT########
#################################

sizeGrWindow(12,6)
pdf(file = "PlotsBlockwiseGeneDendrosAndColors.pdf", wi = 12, he = 6);
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:2), 2, 1), heights = c(0.8, 0.2), widths = c(1,1))
#layout.show(4);
nBlocks = length(signednet$dendrograms)
# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nBlocks) {
  plotDendroAndColors(signednet$dendrograms[[block]], moduleColors[signednet$blockGenes[[block]]],
                      "Module colors",
                      main = paste("Gene dendrogram and module colors" ),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)
}
dev.off()
#################################
##########DENDROGRAM PLOT END####
#################################


#################################
######META-MODULE DETECTION######
#####################################################
### THREE PARTS EXIST:                            ### 
### 1. EIGENGENE DISSIMILARITY CONSESUS MATRIX    ### 
### 2. MDS and k-means CLUSTERING                 ### 
### 3. CREATE METAMODULES BASED ON CLUSTERING END ### 
#####################################################

###############################################
######EIGENGENE DISSIMILARITY CONSESUS MATRIX #
###############################################

#Drop grey module (unassigned genes)
for (set in 1:nSets){
  consMEs[[set]]$data<-consMEs[[set]]$data[,!names(consMEs[[set]]$data) == 'ME0']
}

#Compute eigengene network (AdjME) and Eigengene DissimilarityMatrix
AdjME<-vector(mode="list",length = nSets)
names(AdjME)<-setLabels

for (set in 1:nSets){
  AdjME[[set]] = (1+cor(consMEs[[set]]$data,consMEs[[set]]$data))/2
}


#Consesus of AdjMEs
size<-dim(AdjME[[1]])[1]
ConsME<-matrix(0,size,size)
rownames(ConsME)<-rownames(AdjME[[1]])
colnames(ConsME)<-colnames(AdjME[[1]])

for (i in 1:size){
  for (j in 1:size){
    ConsME[i,j]<-min(AdjME[[1]][i,j],AdjME[[2]][i,j],AdjME[[3]][i,j],AdjME[[4]][i,j])
  }
}
DissConsME<-1-ConsME
heatmap.2(DissConsME,trace="none",dendrogram = "none",col="redblue",
          main="Dissimilarity of Consesus Eigengene Network")

###############################################
######EIGENGENE DISSIMILARITY CONSESUS END#####
###############################################


###################################################
######MDS ON DISSCONSME AND k-MEANS CLUSTERING#####
###################################################

d <- dist(DissConsME) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim

x <- fit$points[,1]
y <- fit$points[,2]
comp<-cbind(x,y)

###K MEANS CLUSTERING#####
# Determine number of clusters (Look at the "elbow" in the scree plot)
wss <- (nrow(comp)-1)*sum(apply(comp,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(comp,nstart=25, iter.max=1000,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 

# From scree plot elbow occurs at k = x
# Apply k-means with k=x
k <- kmeans(comp, 6, nstart=25, iter.max=1000)
plot(comp[,1:2], col=labels2colors(k$cluster), pch=16,ylim=c(-2,2),xlim=c(-1.5,1.5),xlab="Coordinate 1",ylab="Coordinate 2",cex=2)+
  title("MDS and k-means for metamodule identification")
text(x+0.1, y-0.1, labels = row.names(DissConsME), cex=1) 
dev.off()
###############################################
######MDS ON ME AND k-MEANS CLUSTERING END#####
###############################################


###################################################
######CREATE METAMODULES BASED ON CLUSTERING#######
###################################################


# Get the original ME labels
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
changed.moduleColors<-moduleColors

#Merges the ME by declaring the same color for set of MEs
for (clus in 1:length(table(k$cluster))){
  group<-k$cluster==clus
  colors.to.change<-MEColors[which(group==TRUE)]
  changed.moduleColors[changed.moduleColors %in% colors.to.change] <- clus
}

check<-changed.moduleColors=="grey"
changed.moduleColors[check]<- NA

changed.moduleColors<-labels2colors(changed.moduleColors, zeroIsGrey = F,naColor = "grey")

metaMEs<-vector(mode="list",length = nSets)

metaMEs = multiSetMEs(multiExpr, universalColors = changed.moduleColors,excludeGrey = T)

#Rename group names of metaMEs data
for (set in 1:nSets){
names(metaMEs[[set]]$data)<-c("MMblue","MMbrown","MMgreen","MMred","MMturquoise","MMyellow" )
}

###################################################
###CREATE METAMODULES BASED ON CLUSTERING END######
###################################################


#################################
#####META-MODULE DETECTION END###
#################################




#########################################################
#######DIFFERENTIAL ANALYSIS OF EIGENGENE NETWORKS#######
#########################################################


#Plots several ways of eigengene network preservation statistics
#Always create the consMEsC variable 
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors)
MET = consensusOrderMEs(metaMEs)
sizeGrWindow(20,15)
par(cex = 0.8)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1),
                                    marHeatmap = c(2,3,2,1),
                      letterSubPlots = TRUE,
                      zlimPreservation = c(0.1, 1), 
                                    xLabelsAngle = 90)
dev.off()

#########################################################
######DIFFERENTIAL ANALYSIS OF EIGENGENE NETWORKS END####
#########################################################



#######################################################
##    MODULE TRAIT RELATIONSHIP AND DRIVER GENES     ##
#######################################################
#######################################################
### THREE PARTS EXIST:                              ### 
### 1. GET GENE SIGNIFICANCE AND CONNECTIVITY       ### 
### 2. CONSESUS OF GS AND GC                        ### 
### 3. PLOTS of GS vs GC ; ALL TISSUES AND CONSESUS ### 
#######################################################


###############################################
###  GET GENE SIGNIFICANCE AND CONNECTIVITY ###
###############################################

datKME<-vector(mode="list",length = nSets)
geneTraitSignificance<-vector(mode="list",length = nSets)
geneConnectivity<-vector(mode="list",length = nSets)
colorlevels=unique(changed.moduleColors)


#Calculate intramodular connectivity
for (set in 1:nSets){
ADJ=((1+cor(multiExpr[[set]]$data,use="p"))/2)^20
Alldegrees=intramodularConnectivity(ADJ, changed.moduleColors)
DegWithin<-Alldegrees$kWithin
geneConnectivity[[set]]<-DegWithin
}
names(geneConnectivity)<-setLabels
rm(ADJ)



for (set in 1:nSets) {
  
  #Gives the correlation between a gene and the respective module
  
  datKME[[set]]=signedKME(multiExpr[[set]]$data, metaMEs[[set]]$data, outputColumnName="ME")
  
  #Find gene significance accoridung to the Trait
  
  geneTraitSignificance[[set]] =list(data= as.data.frame
                                     (cor((Traits[[set]]$data),multiExpr[[set]]$data, use = "p")))
  geneTraitSignificance[[set]]$data =abs(geneTraitSignificance[[set]]$data)
}

names(datKME)<-setLabels
names(geneTraitSignificance)<-setLabels


###################################################
###  GET GENE SIGNIFICANCE AND CONNECTIVITY END ###
###################################################


#######################################################
###  CONSESUS OF GENE SIGNIFICANCE AND CONNECTIVITY ###
#######################################################


#Consesus of geneSignificance
size<-dim( geneTraitSignificance[[1]]$data)
ConsGS<-matrix(0,size[1],size[2])
rownames(ConsGS)<-rownames(geneTraitSignificance[[1]]$data)
colnames(ConsGS)<-colnames(geneTraitSignificance[[1]]$data)

for (i in 1:size[1]){
  for (j in 1:size[2]){
    buffer<-c(geneTraitSignificance[[1]]$data[i,j],geneTraitSignificance[[2]]$data[i,j],
              geneTraitSignificance[[3]]$data[i,j],geneTraitSignificance[[4]]$data[i,j])
    ConsGS[i,j]<-median(buffer)
  }
}
ConsGS<-t(ConsGS)

#Consesus of geneConnectivity
size<-length( geneConnectivity[[1]])
ConsGC<-matrix(0,size,1)
rownames(ConsGC)<-rownames(geneConnectivity[[1]])
colnames(ConsGC)<-colnames(geneConnectivity[[1]])

for (i in 1:size){
  buffer<-c(geneConnectivity[[1]][i],geneConnectivity[[2]][i],geneConnectivity[[3]][i],geneConnectivity[[4]][i])
  ConsGC[i]<-median(buffer)
}

###########################################################
###  CONSESUS OF GENE SIGNIFICANCE AND CONNECTIVITY END ###
###########################################################

#######################################################
####MODULE TRAIT RELATIONSHIP AND DRIVER GENES END#####
#######################################################



####################################################
### PLOTS of GS vs GC ; ALL TISSUES AND CONSESUS ### 
####################################################


# GS vs. GC for each module - tissue specific

#Set tissue index: 1. Meniscus, 2. Subchondral bone, 3.Synovium, 4. Cartilage 
tissueidx<-4

par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,10,1))
for (i in c(1:length(colorlevels))) {
  whichmodule=colorlevels[[i]]; 
  restrict1 = (changed.moduleColors==whichmodule);
  verboseScatterplot(geneConnectivity[[tissueidx]][restrict1], 
                     geneTraitSignificance[[tissueidx]]$data[restrict1], col=changed.moduleColors[restrict1],
                     main=whichmodule, displayAsZero = 1e-15,
                     xlab = "Gene degree", ylab = "Gene Significance", abline = TRUE)
  
}
mtext(names(geneConnectivity)[tissueidx], side = 3, line = -2, outer = TRUE,font=2,cex=1.1)


# GS vs. GC for each module - consensus
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,10,1))
for (i in c(1:length(colorlevels))) {
  whichmodule=colorlevels[[i]];
  restrict1 = (changed.moduleColors==whichmodule);
  verboseScatterplot(ConsGC[restrict1],
                     ConsGS[restrict1], col=changed.moduleColors[restrict1],
                     main=whichmodule, displayAsZero = 1e-15,
                     xlab = "Gene degree", ylab = "Gene Significance", abline = TRUE)
}
mtext("Consesus plot", side = 3, line = -2, outer = TRUE,font=2,cex=1.1)

########################################################
### PLOTS of GS vs GC ; ALL TISSUES AND CONSESUS END ### 
########################################################




########################################################
####    CHOICE OF TOP GENES AND DATA OUTPUT     ########
########################################################
####1. EXTRACT TOP GENES FROM IMPORTANT META-MODULES####
####2. META-MODULE EXPORT                           ####
########################################################


######################################################
#### EXTRACT TOP GENES FROM IMPORTANT META-MODULES####
######################################################
#Check first which MM you want to investigate by 
#looking at the Consesus GS vs. GC plot
#then look at the colorlevels => get index


gene.names<-names(multiExpr[[1]]$data)
genes.red<-vector(mode="list",length=nSets)
genes.turq<-vector(mode="list",length=nSets)

#threshold (top 20%)
thr.qn<-0.2

whichmodule=colorlevels[[7]]; 
restrict = (changed.moduleColors==whichmodule);
for (set in 1:nSets){
  thr.gs<-quantile(abs(geneTraitSignificance[[set]]$data)[restrict],thr.qn)
  thr.deg<-quantile(abs(geneConnectivity[[set]])[restrict],thr.qn)
  genes.in.module<-gene.names[restrict]
  FilterGenes<-abs(geneTraitSignificance[[set]]$data)[restrict]>as.numeric(thr.gs) & abs(geneConnectivity[[set]])[restrict] > as.numeric(thr.deg)
  genes.red[[set]]<-genes.in.module[FilterGenes]
}


whichmodule=colorlevels[[5]]; 
restrict = (changed.moduleColors==whichmodule);
for (set in 1:nSets){
  thr.gs<-quantile(abs(geneTraitSignificance[[set]]$data)[restrict],thr.qn)
  thr.deg<-quantile(abs(geneConnectivity[[set]])[restrict],thr.qn)
  genes.in.module<-gene.names[restrict]
  FilterGenes<-abs(geneTraitSignificance[[set]]$data)[restrict]>as.numeric(thr.gs) & abs(geneConnectivity[[set]])[restrict] > as.numeric(thr.deg)
  genes.turq[[set]]<-genes.in.module[FilterGenes]
}


x<-unlist(genes.turq)
y<-unlist(genes.red)

x<-unique(x)
y<-unique(y)


write.table(x,"DiseaseSig.txt",quote=FALSE,row.names = F,col.names = F)



############################
###META-MODULE EXPORT#######
############################


#This snippet gets the meta module colors and exports gene names 
#to file

MetaModuleColors<-rownames(table(changed.moduleColors))
MetaModuleColors<-MetaModuleColors[!MetaModuleColors %in% c("grey")]



for (i in 1:length(MetaModuleColors)){
genes.to.extr<-changed.moduleColors==MetaModuleColors[i]
genes<-multiExpr[[3]]$data
genes<-t(genes)
genes<-rownames(genes)
genes<-genes[genes.to.extr]
 filename<-paste('MM',i,'.txt',sep='',collapse=NULL)
 write.table(genes,filename,quote=FALSE,row.names = F,col.names = F)
}


################################
###META-MODULE EXPORT END#######
################################

save.image('WGCNAResults.RData')
