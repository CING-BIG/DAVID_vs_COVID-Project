#######################################################################################
######  DEPARTMENT OF BIOINFORMATICS                                                
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By Margarita Zachariou
######  email: margaritaz@cing.ac.cy
######  The script takes as an input lists of genes from multiple sources which include gene symbols, 
######  entrez ids, original scores and ranks (if available) and normalised scores (weights)
######  All gene lists have been manually curated to annotate which genes can not be found in the Genemania database
###### *  ranked genes lists in terms of absolute logFC from three serum transcriptomic (T) datasets (Series 15, BALF, PBMC)[6],
###### *  ranked genes lists in terms of absolute log fold from one proteomics (P) dataset[5]
###### *  ranked genes lists in terms of p-value from one metabolomics (M) dataset[5]
###### *  unranked list of host proteins (PPI) which interact with SARS-CoV-2 from Gordon et al.[8]
###### *  unranked unique gene list from HPA, excluding the genes identified in Gordon et al.
######
######  The output of this script is used as an input to Cytoscape to generate the genemania network 
######  and as an input to integration_part2 to complete the integration process.
######
######  Developed in R Studio R version 3.6.1 (2020-19-05)
#######################################################################################
# 
# setwd("~/Dropbox/margarita/CING_myshare/R_analysis/covid19_networks/github");

# ----------------------------------------------------------------------------------
# --------------------------------MERGE ALL LISTS-----------------------------------
# ----------------------------------------------------------------------------------
# load individual gene lists from multible sources 
dfBBf <-read.table("integration_data/dfBBfweighted.txt", header = TRUE)
dfBPf <-read.table("integration_data/dfBPfweighted.txt", header = TRUE)
dfB15f <-read.table("integration_data/dfB15fweighted.txt", header = TRUE)
uni_Bid_edited_f<-read.table("integration_data/uni_Bid_edited_f.txt", header = TRUE)
dfs_h_int <-read.table("integration_data/dfs_h_int.txt", header = TRUE)
dfsh_met_int <-read.table("integration_data/dfsh_met_int.txt", header = TRUE)
df_PPI_int <-read.table("integration_data/df_PPI_int.txt", header = TRUE)

dfBBf <-as.data.frame(dfBBf)
dfBPf <-as.data.frame(dfBPf)
dfB15f  <-as.data.frame(dfB15f)
uni_Bid_edited_f <-as.data.frame(uni_Bid_edited_f)
dfs_h_int <-as.data.frame(dfs_h_int)
dfsh_met_int <-as.data.frame(dfsh_met_int)
df_PPI_int <-as.data.frame(df_PPI_int)

# Merge the three transcriptomics lists 
df_merg_BBBP<-merge(dfBBf, dfBPf, by.x="eid",by.y="eid",all.x=TRUE,all.y = TRUE)
colnames(df_merg_BBBP)<-c("eid","gene.BB","abslogFC.BB", "rank.BB", "weightFC.BB","weightR.BB","gene.BP","abslogFC.BP","rank.BP","weightFC.BP", "weightR.BP")

df_merg_B<-merge(df_merg_BBBP, dfB15f, by.x="eid",by.y="eid",all.x=TRUE,all.y = TRUE)
colnames(df_merg_B)<-c("eid","gene.BB","abslogFC.BB", "rank.BB", "weightFC.BB","weightR.BB",
                       "gene.BP","abslogFC.BP","rank.BP","weightFC.BP", "weightR.BP",
                       "gene.B15","abslogFC.B15","rank.B15","weightFC.B15", "weightR.B15")
df_merg_B_int0<-df_merg_B[,c("eid","weightFC.BB","weightFC.BP","weightFC.B15")]
df_merg_B_int0[is.na(df_merg_B_int0)] <- 0

df_merg_B_int<-merge(df_merg_B_int0, uni_Bid_edited_f, by.x="eid",by.y="Beid",all.x=TRUE,all.y = TRUE)

wB=1.0/3.0
df_merg_B_int$BAvg<-wB*df_merg_B_int$weightFC.BB+wB*df_merg_B_int$weightFC.BP+wB*df_merg_B_int$weightFC.B15

# Merge Proteomics and Metabolomics
df_merg_P_M<-merge(dfs_h_int, dfsh_met_int, by.x="Peid",by.y="Meid",all.x=TRUE,all.y = TRUE)
colnames(df_merg_P_M)[1]<-"eid"

df_merg_P_M_int <-df_merg_P_M

# Merge Proteomics and Metabolomics and Transcriptomics
df_merg_B_P_M<-merge(df_merg_B_int, df_merg_P_M_int, by.x="eid",by.y="eid",all.x=TRUE,all.y = TRUE)

# Merge Proteomics and Metabolomics and Transcriptomics and PPI
df_merg_B_P_M_PPI<-merge(df_merg_B_P_M, df_PPI_int, by.x="eid",by.y="PPIeid",all.x=TRUE,all.y = TRUE)

df_merge_final<-subset(df_merg_B_P_M_PPI,select=c(eid,BidGM,Bgene,weightFC.BB,weightFC.BP,weightFC.B15,BAvg,Pprotein,PidGM,weightFC.P,MidGM,Upstream.regulator.M, weightR.M,PPIGenes,weight,PPIidGM ))

df_merge_final$SYS <-as.character(df_merge_final$Bgene)  # your new merged column start with x
df_merge_final$SYS[which(!is.na(df_merge_final$Pprotein))] <- as.character(df_merge_final$Pprotein[!is.na(df_merge_final$Pprotein)])
df_merge_final$SYS[which(!is.na(df_merge_final$Upstream.regulator.M))] <-as.character(df_merge_final$Upstream.regulator.M[!is.na(df_merge_final$Upstream.regulator.M)])  # merge with z
df_merge_final$SYS[which(!is.na(df_merge_final$PPIGenes))] <- as.character(df_merge_final$PPIGenes[!is.na(df_merge_final$PPIGenes)])

df_merge_final$GMid <-as.character(df_merge_final$BidGM)  # your new merged column start with x
df_merge_final$GMid[which(!is.na(df_merge_final$PidGM))] <- as.character(df_merge_final$PidGM[!is.na(df_merge_final$PidGM)])
df_merge_final$GMid[which(!is.na(df_merge_final$MidGM))] <-as.character(df_merge_final$MidGM[!is.na(df_merge_final$MidGM)])  # merge with z
df_merge_final$GMid[which(!is.na(df_merge_final$PPIidGM))] <- as.character(df_merge_final$PPIidGM[!is.na(df_merge_final$PPIidGM)])

df_merge_final_num<-subset(df_merge_final, select=c("weightFC.BB","weightFC.BP","weightFC.B15","BAvg","weightFC.P","weightR.M","weight"))
df_merge_final_num[is.na(df_merge_final_num)] <- 0.0

# Dataframe with merged lists of genes and associated scores
df_int<-cbind.data.frame(eid=as.character(trimws(df_merge_final$eid)),
                         GMid=as.character(trimws(df_merge_final$GMid)),
                         SYS=as.character(trimws(df_merge_final$SYS)),
                         weightFC.BB=as.numeric(df_merge_final_num$weightFC.BB),
                         weightFC.BP=as.numeric(df_merge_final_num$weightFC.BP),
                         weightFC.B15=as.numeric(df_merge_final_num$weightFC.B15),
                         BAvg=as.numeric(df_merge_final_num$BAvg),
                         weightFC.P=as.numeric(df_merge_final_num$weightFC.P),
                         weightR.M=as.numeric(df_merge_final_num$weightR.M),
                         weight=as.numeric(df_merge_final_num$weight),
                         stringsAsFactors=FALSE)

# Filter out the rows which are not recognized by Genemania
NONGM<-length(which(is.na(df_int$GMid)==TRUE))
row.has.na <- apply(df_int, 1, function(x){any(is.na(x))})

NONGM==sum(row.has.na)
df_int_filtered <- df_int[!row.has.na,]

# Add the Human Protein Atlas protein list as a separate columng
HPA4sys<-c("ACE2","TMPRSS2", "CTSB", "CTSL")
ind_HPA4s<-which(df_int_filtered$SYS %in% HPA4sys)
HPA4GMid<-df_int_filtered$GMid[ind_HPA4s]
HPA4eid<-df_int_filtered$eid[ind_HPA4s]
df_int_filtered$weight.HPA4s<-rep(0,length(df_int_filtered$eid))
df_int_filtered$weight.HPA4s[ind_HPA4s]<-1.0
df_int_filtered$weight[ind_HPA4s]<-0.0

# ----------------------------------------------------------------------------------
# --------------------------------NODE SCORE----------------------------------------
# ----------------------------------------------------------------------------------

# weights for the integration
wTf<-0.4
wPf<-0.35
wMf<-0.1
wPPIf<-0.05
wHPAf<-0.1

print(paste("weight sum =",wTf+wPf+wMf+wPPIf+wHPAf))

Inode<- wTf*df_int_filtered$BAvg +wPf*df_int_filtered$weightFC.P+wMf*df_int_filtered$weightR.M +wPPIf*df_int_filtered$weight+wHPAf*df_int_filtered$weight.HPA4s

df_int_filtered$Inode<-Inode

ordI_filt<-df_int_filtered[order(df_int_filtered$Inode,decreasing = T),]
rownames(ordI_filt) <- 1 : length(rownames(ordI_filt))

# Select the top 1001 genes to build the Genemania network in Cytoscape
GM_tableI<-ordI_filt[1:1001,]

GM_genelist<-ordI_filt$GMid[1:1001]

# to be used as input to IntegrationNet.R
write.table(GM_tableI,"integration_data/GM_table_Inode.txt",quote=FALSE, sep='\t', col.names = TRUE,row.names = FALSE)

# to be used as input to integration_part2
write.table(df_int_filtered,"integration_data/df_int_filtered.txt",quote=FALSE, sep='\t', col.names = TRUE,row.names = FALSE)

# to be used as input to Genemania
write.table(GM_genelist,"integration_data/Inode_top1001_GMid.txt",quote=FALSE, sep='\t', col.names = FALSE,row.names = FALSE)


