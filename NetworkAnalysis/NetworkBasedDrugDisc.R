######################Authors : Michael Neidlin          ##################### 
#                               Smaragda Dimitrakopoulou                     #
#                               11-06-2019                                   #
##############################################################################


#This script identifies agglomeration measures of a gene list in a global 
#interactome as presented by Menche et al. 2015 and Guney et al.2016

#FOLLOWING SCENARIOS ARE CONSIDERED:

#1. How agglomerated is a gene list in a background network?
# 1.1 - Largest connected component + z-score/p-value
# 1.2 - Mean shortest distance + z-score/p-value

#2. How separated are two groups of genes in a background network?
# 2.1 - Separation measure for gene list and list of disease genes (299)

#3. How close is a drug target to the gene list?
#3.1 - Mean shortest distance between two groups + z-score/p-value for drug-target list (238)

#INPUT needed:
#Background network, interactome.tsv - given in ENTREZ IDs, undirected network - interactome.tsv
#gene list of interest - DiseaseSig.txt
#drug-gene association list - drugbank2019.csv => Parsed from the drugbank xml file =>DrugbankParse.R

#OUTPUT: AGGLOMERATION MEASURES (LCC AND MSD) AND Z-SCORES
#OUTPUT: LIST OF DRUG DISEASE PROXIMITY

####################################################################

library(igraph)
#====================READ DATA===================
#READ INTERACTOME
df<-read.table('interactome.tsv',sep='\t',header=FALSE)
df<-df[,-3]
net <- graph_from_data_frame(d=df, directed=F) 
#Remove self-loops and multiple edges
net<-simplify(net)


#READ DRUG GENE LIST
df_dr<-read.table('drugbank2019.csv',sep='\t',header=T)
df_dr<-df_dr[,-1]

#READ GENE LIST
genes<-read.table('DiseaseSig.txt')
genes<-as.character(genes$V1)


# load
library('org.Hs.eg.db')

# use mapIds method to obtain Entrez IDs
gene_map<-mapIds(org.Hs.eg.db, genes,  'ENTREZID','SYMBOL',multiVals = list)
gene_map<-as.character(gene_map)


#get all node names of background network
nodes<-unique(c(as.character(df$V1),as.character(df$V2)))
nodes<-as.data.frame(nodes)

#reduce gene_map to just have IDs that are in the background network
gene_map<-gene_map[gene_map %in% nodes$nodes]

#====================READ DATA END===================


###############################################
###           ADDITIONAL FUNCTIONS          ###
###############################################


#===============LARGEST CONNECTED COMPONENT=============
#Takes background network and gene list, creates subgraph and computes
#largest connected components.

LargestConnectedComponent<-function(G,vertices){

sub_net<-induced_subgraph(G,vertices)
all_comp<-components(sub_net)
return<-max(all_comp$csize)
}


#Network sampling (vertices degree distribution) from background (G degreee distribution)
SampleNetworkBG<-function (G.dg,vertices.dg){
  rnd_genes<-vector()
  for (i in 1:length(vertices.dg)){
    dg.group<-names(vertices.dg[i])
    num.group<-as.numeric(vertices.dg[i])
    bg.with.deg<-names(G.dg[G.dg==dg.group])
    rnd_genes<-c(rnd_genes,sample(bg.with.deg, num.group, replace = FALSE, prob = NULL))
  }
  return(rnd_genes)
}

#Does the statistical comparison by running lcc computation choosing random
#nodes for node number same as gene list

LargestConnectedComponentRandom<-function (G,vertices,loops){
  lcc_rand_tot<-vector()
  dg.bg<-degree(G)
  dg.sub<-degree(G,vertices)
  dg.sub.dist<-table(dg.sub)
  
  for (i in(1:loops)){
    rnd_genes<-SampleNetworkBG(dg.bg,dg.sub.dist)
    lcc_rand<-LargestConnectedComponent(G,rnd_genes)
    lcc_rand_tot<-c(lcc_rand_tot,lcc_rand)
  }
  return(lcc_rand_tot)
}

#===============LARGEST CONNECTED COMPONENT END=========


#===============NETWORK BASED DISTANCE================
#Compute mean shortest distance of gene in a list 
#related to any of the proteins in the same gene list


MeanShortestDistance<-function(G,vertices){
  dist<-distances(G,v=vertices,to=vertices)
  dist[dist == 0] <- Inf
  ds<-apply(dist,1,function(x) min(x))
  msd<-mean(ds[is.finite(ds)])
  return(msd)
}

#Does the statistical comparison by running msd computation with random creation 
#of gene lists by sampling from a background network bg.nw

MeanShortestDistanceRandom<-function (G,vertices,loops){
  msd_rand_tot<-vector()
  dg.bg<-degree(G)
  dg.sub<-degree(G,vertices)
  dg.sub.dist<-table(dg.sub)
  for (i in(1:loops)){
    rnd_genes<-SampleNetworkBG(dg.bg,dg.sub.dist)
    msd_rand<-MeanShortestDistance(G,rnd_genes)
    msd_rand_tot<-c(msd_rand_tot,msd_rand)
  }
  return(msd_rand_tot)
}


#===============NETWORK BASED DISTANCE END=============


#============DRUG DISEASE PROXIMITY=====================

#Scripts to compute drug disease proximity as the average of shortest
#paths from drug protein set (verticesB) to disease protein set (verticesA)

DrugDiseaseProximity<-function(G,verticesA,verticesB){
  dist<-distances(G,v=verticesA,to=verticesB)
  dsb<-apply(dist,2,function(x) min(x))
  ds<-dsb
  msd<-mean(ds[is.finite(ds)])
  return(msd)
}

DrugDiseaseProximityRandom<-function(G,verticesA,verticesB,loops){
  d.C_rand_tot<-vector()
  dg.bg<-degree(G)
  dg.subA<-degree(G,verticesA)
  dg.sub.distA<-table(dg.subA)
  dg.subB<-degree(G,verticesB)
  dg.sub.distB<-table(dg.subB)
  for (i in(1:loops)){
    rnd_genesA<-SampleNetworkBG(dg.bg,dg.sub.distA)
    rnd_genesB<-SampleNetworkBG(dg.bg,dg.sub.distB)
    #rnd_genesA<-sample(as.character(nodes$nodes), length(verticesA), replace = FALSE, prob = NULL)
    #rnd_genesB<-sample(as.character(nodes$nodes), length(verticesB), replace = FALSE, prob = NULL)
    d.C_rand<-DrugDiseaseProximity(G,rnd_genesA,rnd_genesB)
    d.C_rand_tot<-c(d.C_rand_tot,d.C_rand)
  }
  return(d.C_rand_tot)
}

#============DRUG DISEASE PROXIMITY END=================



#============MISC FUNCTIONS=============================

ComputeZScore<-function(measurement,distribution){
  z.score<-(measurement-mean(distribution))/(sd(distribution))
}



#===========MISC FUNCTIONS END=========================


###############################################
###           ADDITIONAL FUNCTIONS  END      ###
###############################################


###############################################
#        RUN THE NETWORK ANALYSIS             #
###############################################


##################################
#LARGEST CONNECTED COMPONENT PART#
##################################

#LCC is computed and two distributions 1: sampling from bg interactome, 2: sampling
#from all wgcna nodes

lcc<-LargestConnectedComponent(net,gene_map)
lcc_rnd.bg<-LargestConnectedComponentRandom(net,gene_map,10000)

lcc.z_bg<-ComputeZScore(lcc,lcc_rnd.bg)
lcc.p_bg<-pnorm(-abs(lcc.z_bg))

cat('The z-score for the lcc against interactome is ', lcc.z_bg)
cat('The p-value for the lcc against interactome is ', lcc.p_bg)


##################################
#LARGEST CONNECTED COMPONENT END #
##################################


#############################
#MEAN SHORTEST DISTANCE PART#
#############################

msd<-MeanShortestDistance(net,gene_map)
msd_rnd.bg<-MeanShortestDistanceRandom(net,gene_map,10000)

msd.z_bg<-ComputeZScore(msd,msd_rnd.bg)
msd.p_bg<-pnorm(-abs(msd.z_bg))

cat('The z-score for the msd against background interactome is ', msd.z_bg)
cat('The p-value for the msd against background interactome is ', msd.p_bg)



#############################
#MEAN SHORTEST DISTANCE END #
#############################


############################
#DRUG DISEASE PROXIMITY    #
############################
#Drug disease distances
d.C_tot<-vector()

#Start looping through entire drug list
for (i in 1:dim(df_dr)[1]){
  cat(i)
  
  #reduce genes_dis to just have IDs that are in the background network
  
  genes_dis<-(df_dr$target[i])
  genes_dis<-strsplit(as.character(genes_dis)," ")
  genes_dis<-unlist(genes_dis)
  genes_dis<-genes_dis[genes_dis %in% nodes$nodes]
  genes_dis<-unique(genes_dis)
  
  d.C<-DrugDiseaseProximity(net,gene_map,genes_dis)
  
  d.C_tot<-c(d.C_tot,d.C)
}
df_dC<-as.data.frame(cbind(as.character(df_dr$name),d.C_tot))


#Take the lowest 5% of distance
x<-df_dC$d.C_tot
x<-as.matrix(x)
x<-as.numeric(x)
x<- x[!is.na(x)]
x<-quantile(x,0.05)

df_dr.red<-df_dC[as.character(df_dC$d.C_tot)<=x,]
colnames(df_dr.red)<-c('Drug Name','Distance')
df_dr.red<-df_dr.red[order(df_dr.red$Distance),] 

#Map df_dr to df_dr.red for the RANDOM distribution calculation (saves time)
df_dr2<-df_dr[as.factor(df_dr.red$`Drug Name`),]
d.C_z_tot<-vector()

#Start looping through entire drug list
for (i in 1:dim(df_dr2)[1]){
  cat(i)
  #reduce genes_dis to just have IDs that are in the background network
  
  genes_dis<-(df_dr2$target[i])
  genes_dis<-strsplit(as.character(genes_dis)," ")
  genes_dis<-unlist(genes_dis)
  genes_dis<-genes_dis[genes_dis %in% nodes$nodes]
  genes_dis<-unique(genes_dis)
  
  d.C_rnd<-DrugDiseaseProximityRandom(net,gene_map,genes_dis,1000)
  d.C_rnd<-d.C_rnd[is.finite(d.C_rnd)]
  d.C_z<-(d.C_tot[i]-mean(d.C_rnd))/sd(d.C_rnd)
  d.C_z_tot<-c(d.C_z_tot,d.C_z)
}


df_dC2<-as.data.frame(cbind(df_dr.red,d.C_z_tot))
df_dC2<-df_dC2[!is.na(df_dC2$d.C_z_tot),]

#=======================END=============================
df_dC3<-df_dC2
df_dC3<-df_dC3[df_dC3$d.C_z_tot< -2,]

save.image('NetworkAnalysis.RData')
write.table(df_dC3,'Drugs_from_db2019.csv',sep=',',row.names = F)



