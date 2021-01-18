###################################################################################
######  DEPARTMENT OF BIOINFORMATICS                                                
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By Marios Tomazou
######  email: mariost@cing.ac.cy
######
######  This script takes as an input the drug list indicating the shortlisted   
######  drugs, a table of drugs indicating their originating list and their viral 
######  target proteins.                                                         
######  The output is a sankey plot visualising the mapping of lists to          
######  shortlisted drugs to genes to KEGG pathways.                             
######
######  Developed in R Studio R version 3.6.1 (2019-07-05)
###################################################################################

library(data.table)
library(dplyr)
library(igraph)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(ggalluvial)
library(stringr)
library(gridExtra)
library(ggrepel)


ldat<-read.csv("Top66Drugs-Organism.csv",header = T,stringsAsFactors = F)
sdr<-read.csv("SDRUG_Table.csv",header = T,stringsAsFactors = F)

ldat$Within.Top.20.List.ID[which(ldat2$Within.Top.20.List.ID=="")]<-ldat2$List.ID[which(ldat2$Within.Top.20.List.ID=="")]
ldat<-ldat[which(ldat$Target.Organism!="Human"),]
#Prepare a data frame wich maps Drugs to Gene targets to Pathawys and Pathways ID
DF<-data.frame(Drug=NA,List=NA,Target=NA)

for (i in 1:length(unique(ldat$Drug.Name))){
  d<-unique(ldat$Drug.Name)[i]
  print(d)
  
  tar<-unlist(strsplit(sdr$ALTtarget[which(sdr$Drug.Name==ldat$Drug.Name[i])],";"))

  if(length(tar)>0){
  for (j in 1:length(tar)){
    lis<-str_remove(unlist(strsplit(sdr$Within.Top.20.List.ID[which(sdr$Drug.Name==ldat$Drug.Name[i])],",")),"\\s+")
    for(z in 1:length(lis)){
    DF<-rbind.data.frame(DF,data.frame(Drug=ldat$Drug.Name[i],List=lis[z],Target=tar[j]))
    }
  }
  }
  
}

#Add the originating lists
DF<-DF[-1,]
DF$Freq<-1


DF3<-DF
DF3$Drug2<-DF3$Drug
DF3$do<-0
for (i in 1:nrow(DF3)){
  DF3$do[i]<-sdr$Max.Norm.Ranking[match(DF3$Drug2[i],sdr$Drug.Name)]
}
DF3<-DF3[order(DF3$do,decreasing = T),]

#extract and order lists,genes,drugs and pathways
list<-as.character(unique(DF3$List))
list<-list[order(list)]
tar<-as.character(unique(DF3$Target))
tar<-tar[order(tar)]
drug<-as.character(unique(DF3$Drug))
drug<-drug[order(sdr$Max.Norm.Ranking[match(drug,sdr$Drug.Name)],decreasing = T)]


bb <- to_lodes_form(DF3, key = "Demographic", axes = c(2,1,3))
# calculate within-group proportions
bb <- na.omit(bb)
bb <- group_by(bb, Demographic)
bb <- add_count(bb, Freq, name = "Total")
bb <- ungroup(bb)
bb$Freq<-as.numeric(bb$Freq)
bb <- transform(bb, Prop = Freq / Total)

# alluvial diagram of within-axis proportions
bb$stratum<-ordered(bb$stratum,levels=c(list,drug,tar))

for (j in unique(bb$Demographic)){
  us<-unique(bb$stratum[bb$Demographic==j])
  l<-length(unique(bb$stratum[bb$Demographic==j]))
  sl<-length((bb$stratum[bb$Demographic==j]))
  
  for(i in 1:l){
    
    bb$Prop[which(bb$stratum==us[i] & bb$Demographic==j)]<-1/(l*length(bb$Prop[which(bb$stratum==us[i] & bb$Demographic==j)]))
  }
}

bb$col<-as.character(bb$stratum)
c<-0

for (i in unique(bb$stratum[bb$Demographic=="List"])){
  ud<-unique(bb$stratum[which(!is.na(match(bb$alluvium,bb$alluvium[bb$stratum==i])) & bb$Demographic=="Drug")])
  if(length(ud)==1){
    c<-c+1
    bb$col[which(bb$Demographic=="List" & bb$stratum==i)]<-as.character(ud)
  }
  else{bb$col[which(bb$Demographic=="List" & bb$stratum==i)]<-"A"}
}

for (i in unique(bb$stratum[bb$Demographic=="Target"])){
  ud<-unique(bb$stratum[which(!is.na(match(bb$alluvium,bb$alluvium[bb$stratum==i])) & bb$Demographic=="Drug")])
  if(length(ud)==1){
    c<-c+1
    bb$col[which(bb$Demographic=="Target" & bb$stratum==i)]<-as.character(ud)
  }
  else{bb$col[which(bb$Demographic=="Target" & bb$stratum==i)]<-"A"}
}

# Plot alluvial

ggplot(bb, aes(x = Demographic, stratum = stratum, alluvium = alluvium, y = Prop, label = stratum, fill=col)) +
  geom_flow(aes(fill = Drug2)) +
  geom_stratum(aes(color="darkgrey")) +
  geom_text(stat = "stratum") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())+
  guides(fill = FALSE) +
  scale_fill_brewer(type = "qual", palette = "Paired")+
  scale_color_brewer(type = "qual", palette = "Paired")# +
