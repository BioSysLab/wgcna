#This script parses the drugbank xml file and exports the information
#to a csv file in order to be postprocessed with network analysis
#Further on it changes the dataframe in order to have a data frame to 
#be used for network analysis - 

library(drugbankR)
library(stringr)

#########PARSE THE XML AND QUERY THE DATABASE###########
#OUTPUT: df
drugbank_dataframe <- dbxml2df(xmlfile="drugbank_db.xml", version="5.1.3")
df2SQLite(dbdf=drugbank_dataframe, version="5.1.3")


all <- queryDB(type = "getAll", db_path="drugbank_5.1.3.db") # get the entire drugbank dataframe

ids <- queryDB(type = "getIDs", db_path="drugbank_5.1.3.db") # get all the drugbank ids

# given drugbank ids, determine whether they are FDA approved,extract them
appr<-queryDB(ids = ids,type = "whichFDA", db_path="drugbank_5.1.3.db") 
appr<-appr[appr$whichFDA==TRUE,1]

# given drugbank ids, get their targets
trg<-queryDB(ids = ids,type = "getTargets", db_path="drugbank_5.1.3.db") 


#Do reordering of trg dataframe and mapping of Uniprot to ENTREZ ID
df<-drugbank_dataframe[drugbank_dataframe$`drugbank-id` %in% trg$q_db_id,1:2]
df<-cbind(df,trg$t_Uni_id)
colnames(df)<-c('drugbank_id','name','target')
df<-df[!is.na(df$target),]

#########PARSE THE XML AND QUERY THE DATABASE END###########


##########CHANGE STRUCTURE OF DF AND MAP FROM UNIPROT TO ENTREZ######
library(tidyr)
df.sep<-separate_rows(df, target)

# load
library('org.Hs.eg.db')
gene_map<-mapIds(org.Hs.eg.db, as.character(df.sep$target),  'ENTREZID','UNIPROT',multiVals = "first")

df.sep<-cbind(df.sep,gene_map)
df.sep.red<-df.sep[!is.na(df.sep$gene_map),-3]

library(dplyr)
df.group<-df.sep.red %>%
  group_by(drugbank_id,name) %>%
  summarise(gene_map = paste(gene_map, collapse = " "))

colnames(df.group)<-c('drugbank_id','name','target')

##########CHANGE STRUCTURE OF DF AND MAP FROM UNIPROT TO ENTREZ END######

write.table(df.group,'drugbank2019.csv',sep='\t',col.names = T,row.names = F,quote=F)


