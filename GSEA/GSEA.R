######################Authors : Michael Neidlin          ##################### 
#                               Smaragda Dimitrakopoulou                     #
#                               25-06-2019                                   #
##############################################################################

#1. THIS SCRIPT CAN BE USED TO PERFORM GENE SET ENRICHMENT ANALYSIS ON THE 
#GENE EXPRESSION DATA OF THE TISSUES (Individually or in combination)

#2. CODE TO CALCULATE T VALUES AND SUM THEM UP IN A LIST

#3. EXTRACTS THE T VALUES BASED ON THE DISEASE MODULE GENE NAMES AND
#EXPORTS A SORTED LIST (INPUT: MMi.txt - can be found in the WGCNA folder)

#Load files
# LOAD Preprocessing.Rdata from the Preprocessing folder



#OUTPUT: multi.t : List with t values


####################################################################



#######################################################

library(GEOquery)
library(limma)
library(piano)
library(dplyr)



###############################
####ADDITIONAL FUNCTIONS#######
###############################

#####DEGS######

DEGS<-function (x,design) {
  
  fit<-lmFit(x,design)
  fit<-eBayes(fit)
  genes<-topTable(fit,coef = "OAvsControl" , number=Inf, adjust.method = "BH",sort.by = c("logFC"))
  index<-which(genes$adj.P.Val<=0.05 & abs(genes$logFC)>=1.5)
  DEGS<-genes[index,]
  
  DEGS_Analysis = vector(mode = "list", length = 2)
  DEGS_Analysis[[1]] = list(genes=as.data.frame(genes))
  DEGS_Analysis[[2]] = list(DEGS=as.data.frame(DEGS))
  
  return(DEGS_Analysis)
} 
#####DEGS END######

#####DESIGN MATRIX############
designfunction<-function (x) {
  
  Group<-factor(x$Disease,levels=c("OA","Normal"))
  design<-model.matrix(~Group)
  colnames(design)<-c("OA","OAvsControl")
  return(design)
  
}  
#####DESIGN MATRIX END############


#####GSEA########################
Path_function<-function (tvalues,Annotation,GSC) {
  
  
  GeneSetStatist<-runGSA(geneLevelStats = tvalues
                         ,geneSetStat = "gsea",gsc =GSC,gsSizeLim = c(5,500)
                         ,signifMethod = "geneSampling", adjMethod = "fdr")
  flag0<-GSAsummaryTable(GeneSetStatist)
  flag1<-which(flag0$`p (dist.dir.up)`<=0.01 | flag0$`p (dist.dir.dn)`<=0.01)
  
  imp_paths<-as.data.frame(flag0[flag1,])
  imp_paths_names<-imp_paths$Name
  
  #Write the important pathways of each meta-module to .txt
    write.table(imp_paths_names,'Pathways of Meta-Module.txt',row.names=FALSE,quote = FALSE,sep='\t')
  
  return()
}  

#####GSEA END########################

###################################
####ADDITIONAL FUNCTIONS END#######
###################################



#########################################
## COMPUTE DEGS AND T-VALUES ALL TISSUES#
#########################################


#Meniscus

##Make the design matrix
design<-designfunction(MeniscSamples)

##Remove outliers according to the WGCNA hierarchical clustering plot
MenExp<-MenExp[,-c(11,12)]
design<-design[-c(11,12),]

##Differential Expression Analysis
Men_DEGS<-DEGS(MenExp,design)

##Make the gene-level statistics for the Pathway Analysis
Menisc_tvalues<-Men_DEGS[[1]]$genes$t
Menisc_tvalues<-as.data.frame(Menisc_tvalues)
rownames(Menisc_tvalues)<-rownames(Men_DEGS[[1]]$genes)

#Subchondral

##Make the design matrix
design<-designfunction(Subc_Samples)

##Remove outliers according to the WGCNA hierarchical clustering plot
SubExp<-SubExp[,-c(1:5,11:30)]
design<-design[-c(1:5,11:30),]

##Differential Expression Analysis 
Sub_DEGS<-DEGS(SubExp,design)

##Make the gene-level statistics for the Pathway Analysis
Subc_tvalues<-Sub_DEGS[[1]]$genes$t
Subc_tvalues<-as.data.frame(Subc_tvalues)
rownames(Subc_tvalues)<-rownames(Sub_DEGS[[1]]$genes)

#Cartilage

##Make the design matrix
design<-designfunction(CartSamples)

##Remove outliers according to the WGCNA hierarchical clustering plot
CarExp<-CarExp[,-c(11,15,19)]
design<-design[-c(11,15,19),]

##Differential Expression Analysis
Car_DEGS<-DEGS(CarExp,design)

##Make the gene-level statistics for the Pathway Analysis
Car_tvalues<-Car_DEGS[[1]]$genes$t
Car_tvalues<-as.data.frame(Car_tvalues)
rownames(Car_tvalues)<-rownames(Car_DEGS[[1]]$genes)

#Synovium

##Make the design matrix
design<-designfunction(SynSamples)


##Differential Expression Analysis 
Syn_DEGS<-DEGS(SynExp,design)

##Make the gene-level statistics for the Pathway Analysis
Syn_tvalues<-Syn_DEGS[[1]]$genes$t
Syn_tvalues<-as.data.frame(Syn_tvalues)
rownames(Syn_tvalues)<-rownames(Syn_DEGS[[1]]$genes)


#############################################
## COMPUTE DEGS AND T-VALUES ALL TISSUES END#
#############################################

####################################################
#### CREATE OVERLAP OF GENE NAMES FROM ALL DATASETS#
####################################################
keep1<-match(rownames(Syn_tvalues),rownames(Menisc_tvalues))
any(is.na(keep1))
keep2<-match(rownames(Syn_tvalues),rownames(Subc_tvalues))
any(is.na(keep2))
keep3<-match(rownames(Syn_tvalues),rownames(Car_tvalues))
any(is.na(keep3))

keep<-logical(length = length(Syn_tvalues[,1]))
for (i in 1:length(keep1)) {
  if (!is.na(keep1[i]) & !is.na(keep2[i]) & !is.na(keep3[i])) {
    keep[i]=TRUE
  }else{
    keep[i]=FALSE
  }
}

Syn_tvalues<-as.matrix(Syn_tvalues)
Syn_tvalues.red<-Syn_tvalues[keep,]
Syn_tvalues.red<-as.data.frame(Syn_tvalues.red)

keep1<-match(rownames(Syn_tvalues.red),rownames(Menisc_tvalues))
keep2<-match(rownames(Syn_tvalues.red),rownames(Subc_tvalues))
keep3<-match(rownames(Syn_tvalues.red),rownames(Car_tvalues))

Subc_tvalues.red<-Subc_tvalues[keep2,]
Menisc_tvalues.red<-Menisc_tvalues[keep1,]
Car_tvalues.red<-Car_tvalues[keep3,]

Subc_tvalues.red<-as.data.frame(Subc_tvalues.red)
rownames(Subc_tvalues.red)<-rownames(Syn_tvalues.red)
Menisc_tvalues.red<-as.data.frame(Menisc_tvalues.red)
rownames(Menisc_tvalues.red)<-rownames(Syn_tvalues.red)
Car_tvalues.red<-as.data.frame(Car_tvalues.red)
rownames(Car_tvalues.red)<-rownames(Syn_tvalues.red)

nSets=4
setLabels = c("Meniscus", "Subchondral Bone","Synovium","Cartilage")


multi.t = vector(mode = "list", length = nSets)
multi.t[[1]] = list(data =as.data.frame( Menisc_tvalues.red))
multi.t[[2]] = list(data = as.data.frame(Subc_tvalues.red))
multi.t[[3]] = list(data =as.data.frame( Syn_tvalues.red))
multi.t[[4]] = list(data =as.data.frame( Car_tvalues.red))
names(multi.t)<-setLabels
####################################################
#### CREATE OVERLAP OF GENE NAMES FROM ALL DATASETS#
####################################################


###############################
#### GET MEDIAN OF T VALUES####
###############################

##This code creates median of absolute t values of all four tissues

t.all<-cbind(multi.t[[1]]$data,multi.t[[2]]$data,multi.t[[3]]$data,multi.t[[4]]$data)
t.median<-apply(t.all,1,function(x) median(x))
t.median<-as.data.frame(t.median)
t.median<-multi.t[[3]]$data
#Read the genes of the meta-module
genes.to.extr<-read.table('MM5.txt')
genes.to.extr<-as.character(genes.to.extr$V1)

#Extract and sort the t-values of the meta-module genes
genes.t<-t.median[rownames(t.median) %in% genes.to.extr,,drop=F]
genes.t<-genes.t[order(-genes.t$t.median),,drop=F]

write.table(rownames(genes.t),'MM5_tval_sorted.txt',quote=FALSE,row.names = F,col.names = F)



###################################
#### GET MEDIAN OF T VALUES END####
###################################


# Pathway analysis in the meta-modules

##Load the GO and Pathway terms
###Pathways
GSC<-loadGSC(file = "c2.cp.v6.1.symbols.gmt")

###GO terms
library(org.Hs.eg.db)
library(AnnotationDbi)
columns(org.Hs.eg.db)
key<-rownames(SynExp)
Annot_File<-select(org.Hs.eg.db,keys = key,keytype = "SYMBOL",columns = c("GO"))
Annot_File<-Annot_File[,1:2]
keep<-which(!is.na(Annot_File[,2]))
Annot_File<-Annot_File[keep,]
Annot_File<-as.data.frame(Annot_File)


GSC_GO<-loadGSC(Annot_File)

library(GO.db)
Go_keys<-Annot_File$GO
GO_Annotation<-select(GO.db, keys=Go_keys, columns=c("TERM","ONTOLOGY"),
                      keytype="GOID")



#Pathway Analysis in each cluster
Path_function(genes.t,GO_Annotation,GSC_GO) 

GeneSetStatist<-runGSA(geneLevelStats = genes.t
                       ,geneSetStat = "gsea",gsc =GSC,gsSizeLim = c(5,500)
                       ,signifMethod = "geneSampling", adjMethod = "fdr")
flag0<-GSAsummaryTable(GeneSetStatist)
flag1<-which(flag0$`p (dist.dir.up)`<=0.05 | flag0$`p (dist.dir.dn)`<=0.05)

imp_paths<-as.data.frame(flag0[flag1,])

a<-GO_Annotation[GO_Annotation$GOID %in% imp_paths_names,2]
a<-unique(a)
