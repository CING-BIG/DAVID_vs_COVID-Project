#######################################################################################
######  DEPARTMENT OF BIOINFORMATICS                                                
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By Margarita Zachariou
######  email: margaritaz@cing.ac.cy
######  
######  Takes as input the integration based on the node score table (df_int_filtered) and the 
######  Genemania network files
######  The edge score is callculated and combined to the node score to obtain the overall integration score.
######  The output of this script is used as an input to IntegrationNet.R to plot the integration network 
######  and as an input to CODRES and Pathwalks tools.
######  
######  Developed in R Studio R version 3.6.1 (2020-19-05)
#######################################################################################
# ----------------------------------------------------------------------------------
# --------------------------------EDGE SCORE----------------------------------------
# ----------------------------------------------------------------------------------


library(igraph)

df_int_filtered<-read.delim("integration_data/df_int_filtered.txt",header = TRUE,sep = "\t")
GM_genelist0<-read.table("integration_data/Inode_top1001_GMid.txt",header = FALSE, sep='\t',stringsAsFactors = FALSE)
GM_genelist<-GM_genelist0$V1

# LOAD GENEMANIA NETWORK   with entrez ids
GM<-read.delim("integration_data/GM_geneEdgelist_ids.tsv",header = FALSE,sep = "\t")

colnames(GM)<-c("from","to","weight")
dfGM<-cbind.data.frame(from=as.character(GM$from),to=as.character(GM$to),weight=as.numeric(GM$weight),stringsAsFactors=F)

GM_mapping<-read.delim("integration_data/genemania_mapping.tsv",header = T,sep = "\t")
colnames(GM_mapping)<-c("GMmapeid","GMmapSYS","GMmapGMid")

uniqueGM<-unique(c(dfGM$from,dfGM$to))

# Build the genemania network 
gGM<-graph_from_data_frame(dfGM,directed = F)

E(gGM)$weight<-dfGM$weight

length(uniqueGM)==length(V(gGM))

# calculate the weighted degree of the network 
IGGe<-strength(gGM)

ind_filt<-which(df_int_filtered$GMid %in% GM_genelist)
df_IGG<-df_int_filtered[ind_filt,]

IGGne_genelist_1<-cbind.data.frame(GMid=as.character(names(IGGe)),IGGe=as.numeric(IGGe),stringsAsFactors=F)
IGGne_genelist_2<-subset(df_IGG,select=c("eid","SYS","GMid","Inode"))
IGGne_genelist_3<-GM_mapping

IGGne_genelist_12<-merge(IGGne_genelist_1,IGGne_genelist_2, by.x="GMid",by.y="GMid",all.x=TRUE,all.y = TRUE)
IGGne_genelist_123<-merge(IGGne_genelist_12,IGGne_genelist_3, by.x="GMid",by.y="GMmapGMid",all.x=TRUE,all.y = TRUE)
IGGne_genelist<-IGGne_genelist_123[,c("eid","GMid","GMmapSYS","IGGe","Inode")]

IGGne_genelist$IGGe[which(is.na(IGGne_genelist$IGGe))]<-0.0
IGGne_genelist$Inode[which(is.na(IGGne_genelist$Inode))]<-0.0
maxn.att <- apply(IGGne_genelist[4:length(IGGne_genelist)], 2, max)
Features_norm<-sweep(data.matrix(IGGne_genelist[4:length(IGGne_genelist)]), 2, maxn.att,'/')
IGGne_genelist_norm<-IGGne_genelist
IGGne_genelist_norm[4:length(IGGne_genelist)]<-Features_norm

# ----------------------------------------------------------------------------------
# --------------------------------INTEGRATION SCORE---------------------------------
# ----------------------------------------------------------------------------------
w5=0.5
IGGne_05<-w5*IGGne_genelist_norm$Inode + (1.0-w5)*IGGne_genelist_norm$IGGe
IGGne_genelist_final<-IGGne_genelist_norm
IGGne_genelist_final$IGGne_05<-IGGne_05

write.table(IGGne_genelist_final,"integration_data/IGGne_genelist_final.txt",quote=FALSE, sep='\t', col.names = TRUE,row.names = FALSE)

# ----------------------------------------------------------------------------------
# -------------------------------- SAVE FILES FOR CODRES/PATHWALKS------------------
# ----------------------------------------------------------------------------------

IGGne_genelist_final_codres<-IGGne_genelist_norm
IGGne_genelist_final_codres$IGGne05<-IGGne_05


# Load genemania network with genesymbols
GMs<-read.delim("integration_data/GM_geneEdgelist_genesys.tsv",header = FALSE,sep = "\t")

colnames(GMs)<-c("from","to","weight")
dfGMs<-cbind.data.frame(from=as.character(GMs$from),to=as.character(GMs$to),weight=as.numeric(GMs$weight),stringsAsFactors=F)

gGMs<-graph_from_data_frame(dfGMs,directed = F)
E(gGMs)$weight<-dfGMs$weight

ordIGGne_genelist_final_codres<-IGGne_genelist_final_codres[order(IGGne_genelist_final_codres$IGGne05,decreasing = T),]
rownames(ordIGGne_genelist_final_codres) <- 1 : length(rownames(ordIGGne_genelist_final_codres))

ind300<-as.character(ordIGGne_genelist_final_codres[1:300,]$GMmapSYS)


gGMs300<-induced_subgraph(gGMs,v=ind300)

edgelist_gGMs300<-get.data.frame(gGMs300)

# Final Table used as input to Pathwalks
write.table(edgelist_gGMs300, file='integration_data/geneEdgelist_ids_genesys300.tsv', quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)

# Final Table used as input to CODRES
write.table(ordIGGne_genelist_final_codres,"integration_data/IGGne_genelist_final_codres_ordered.tsv",quote=FALSE, sep='\t', col.names = TRUE,row.names = FALSE)


