######################Authors : Michael Neidlin          ##################### 
#                               Smaragda Dimitrakopoulou                     #
#                               25-06-2019                                   #
##############################################################################

#INPUT NEEDED: Download from GEO
#Synovium=GSE46750, Subchondral=GSE51588, Meniscus=GSE98918, Cartilage=GSE57218

#OUTPUT: PreprocessingGSEA.RData - Expression data for GSEA analysis
#OUTPUT: multiExpr - List with normalized gene expression data

#THIS SCRIPT LOADS AND PRE-PROCESSES THE DATA FOR THE WGCNA ANALYSIS


####################################################################




#load packages

library(GEOquery)
library(limma)
library(tidyverse)
library(hgu133a.db)




##############################
###ADDITIONAL FUNCTIONS#######
##############################
Affy_function<-function (x) {
  
   normx<-rma(x)
  Exp<-exprs(normx)
  
  ## Annotation
  # try http:// if https:// URLs are not supported
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("hgu133a.db")
  #library(hgu133a.db)
  
  library(hgu133a.db)
  y<-rownames(Exp)
  Annot<-select(hgu133a.db, keys=y,keytype ="PROBEID",columns = c("ENTREZID","SYMBOL"))
  remove(y)
  
  ##Get rid of the probe sets that are not correspond to genes or correspond to many genes
  
  Annot<-Annot %>% filter(!is.na(Annot$ENTREZID))
  Annot<-Annot %>% group_by(PROBEID) %>% mutate(x=n_distinct(ENTREZID))%>% filter(x==1)
  Annot<-Annot %>% select(1:3)
  index<-match(Annot$PROBEID,rownames(Exp))
  any(is.na(index))
  Exp<-Exp[index,]
  rownames(Exp)<-Annot$SYMBOL
  
  ##Aggregate
  Exp<-aggregate(Exp, by=list(Annot$SYMBOL),FUN="median")
  rownames(Exp)<-Exp$Group.1
  Exp<-Exp %>% select(-1)
  
  return(Exp)
  
  
}  

Agilent_function<-function (x) {
  
  
  
  ##Background Correction and get rid of control probes
  str(x)
  x@.Data[[4]]$ControlType
  table(x@.Data[[4]]$ControlType)
  
  bacx<-backgroundCorrect(x,method ="normexp",normexp.method = "rma")
  bacx<-bacx[bacx@.Data[[3]]$ControlType== 0, ]
  
  ##Normalisation and average expression for replicate probes
  normx<-normalizeBetweenArrays(bacx,method = "quantile")
  final<- avereps(normx, ID=normx[[3]]$ProbeName)
  Exp<-final[[1]]
  Exp<-as.data.frame(Exp)
  
  ##Remove Probes that do not correspond to genes
  Annot<-final[[3]]
  any(is.na(Annot$GeneName))
  index<-!(Annot$ProbeName==Annot$GeneName)
  Annot<-Annot[index,]
  Exp<-Exp[index,]
  
  
  ##Agreggate
  Exp<-aggregate(Exp, by=list(Annot$GeneName),FUN="median")
  rownames(Exp)<-Exp$Group.1
  Exp<-Exp %>% select(-1)    
  
  return(Exp)
  
}


##############################
###ADDITIONAL FUNCTIONS END###
##############################


###################################
###DATA LOADING AND NORMALIZATION #
###################################

#------ Synovium ------#

##Change directory to Synovium Raw data directory (ADJUST IF NECESSARY!!!)
setwd("D:/R_projects/Paper/Preprocessing/GSE55235_RAW")
library(affy)
x<- ReadAffy()
x<-x[,1:20]

##Change the directory
#######################
pathwd<-getwd()

#!!Have to download source(""hgu133a.db"")
SynExp<-Affy_function(x)

##Sample information
SynSamples<-as.data.frame(c(rep("Normal",10),rep("OA",10)))
names(SynSamples)<-"Disease"


#------  Meniscus ------#

##Load the files - ADJUST IF NECESSARY
path.meniscus<-paste(pathwd,"/GSE98918_RAW",sep='')
file.names <- list.files(path =path.meniscus, 
                         pattern = ".txt")
x <- read.maimages(files =file.names, source="agilent",
           green.only=TRUE,path =path.meniscus ,
           columns=list(G="rMedianSignal",Gb="rBGMedianSignal"))
remove(file.names)
MenExp<-Agilent_function(x)

##Sample information
MeniscSamples<-as.data.frame(c(rep("Normal",12),rep("OA",12)))
names(MeniscSamples)<-"Disease"

#------Subchondral Bone------#

##Load the files - ADJUST IF NECESSARY

path.subc<-paste(pathwd,"/GSE51588_RAW",sep='')
file.names <- list.files(path =path.subc, 
                         pattern = ".txt.gz")
x <- read.maimages(files =file.names, source="agilent",
                   green.only=TRUE,path =path.subc)
                   
remove(file.names)
SubExp<-Agilent_function(x)

##Sample Information
Subc_Samples<-as.data.frame(c(rep("Normal",10),rep("OA",40)))
names(Subc_Samples)<-"Disease"

#------ Cartilage ------#

##Load the files - ADJUST IF NECESSARY

path.cart<-paste(pathwd,"/GSE117999_RAW",sep='')
file.names <- list.files(path =path.cart,
                         pattern = ".txt")
x <- read.maimages(files =file.names, source="agilent",
          green.only=TRUE,path =path.cart,
          columns=list(G="rMedianSignal",Gb="rBGMedianSignal"))
remove(file.names)
CarExp<-Agilent_function(x)

##Sample Information
CartSamples<-as.data.frame(c(rep("Normal",12),rep("OA",12)))
names(CartSamples)<-"Disease"

#######################################
###DATA LOADING AND NORMALIZATION END #
#######################################



####################################
##DATA SUMMARY AND EXPORT FOR GSEA##
####################################
##Save them for later use (GSEA)

save(CarExp,CartSamples,MenExp,MeniscSamples,SubExp,Subc_Samples,SynExp,
     SynSamples,file = "PreprocessingGSEA.RData")

########################################
##DATA SUMMARY AND EXPORT FOR GSEA END##
########################################



#########################################
####OUTLIER REMOVAL AND DATA SUMMARY#####
#########################################

#OUTLIER REMOVAL
MenExp<-MenExp[,-c(11,12)]
SubExp<-SubExp[,-c(1:5,11:30)]
CarExp<-CarExp[,-c(11,15,19)]

#DATA SUMMARY - HAVE 4 DATASETS WITH THE SAME GENES
keep1<-match(rownames(SynExp),rownames(MenExp))
any(is.na(keep1))
keep2<-match(rownames(SynExp),rownames(SubExp))
any(is.na(keep2))
keep3<-match(rownames(SynExp),rownames(CarExp))
any(is.na(keep3))

keep<-logical(length = length(SynExp[,1]))
for (i in 1:length(keep1)) {
  if (!is.na(keep1[i]) & !is.na(keep2[i]) & !is.na(keep3[i])) {
    keep[i]=TRUE
  }else{
    keep[i]=FALSE
  }
}

WGCNA_SynExp<-SynExp[keep,]

keep1<-match(rownames(WGCNA_SynExp),rownames(MenExp))
keep2<-match(rownames(WGCNA_SynExp),rownames(SubExp))
keep3<-match(rownames(WGCNA_SynExp),rownames(CarExp))

WGCNA_SubExp<-SubExp[keep2,]
WGCNA_MenisExp<-MenExp[keep1,]
WGCNA_CarExp<-CarExp[keep3,]
remove(keep,keep1,keep2,keep3)

# 4 DATASETS WGCNA_X for WGCNA ANALYSIS

#We have 4 datasets
nSets=4
setLabels = c("Meniscus", "Subchondral Bone","Synovium","Cartilage")

multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data =as.data.frame( t(WGCNA_MenisExp)))
multiExpr[[2]] = list(data = as.data.frame(t(WGCNA_SubExp)))
multiExpr[[3]] = list(data =as.data.frame( t(WGCNA_SynExp)))
multiExpr[[4]] = list(data =as.data.frame( t(WGCNA_CarExp)))
names(multiExpr)<-setLabels


save(multiExpr,file = "multiExpr.RData")
#############################################
####OUTLIER REMOVAL AND DATA SUMMARY END#####
#############################################












