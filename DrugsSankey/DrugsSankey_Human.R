###################################################################################
######  DEPARTMENT OF BIOINFORMATICS                                                
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By Marios Tomazou
######  email: mariost@cing.ac.cy
######
######  This script takes as an input the drug list indicating the shortlisted   
######  drugs, the mapping of KEGG pathways with their IDs, the mapping of genes 
######  identified as drug targets with their KEGG pathways, a table of drugs and
######  their originating list and the output of the PWalks2net script.          
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

ldat2<-read.csv("Top66Drugs-Organism.csv",header = T,stringsAsFactors = F)
ldat<-read.csv("KEGG_drug2GeneNET.csv",header = T,stringsAsFactors = F)
kp<-read.csv("keggpathindex.csv",header=T,stringsAsFactors = F)
sdr<-read.csv("SDRUG_Table.csv",header = T,stringsAsFactors = F)
nods<-read.csv("nods_top300.csv",header = T)

ldat2$Within.Top.20.List.ID[which(ldat2$Within.Top.20.List.ID=="")]<-ldat2$List.ID[which(ldat2$Within.Top.20.List.ID=="")]

#Prepare a data frame wich maps Drugs to Gene targets to Pathawys and Pathways ID
DF<-data.frame(Drug=NA,Gene=NA,Pathway=NA,PathwayID=NA)

for (i in 1:length(unique(ldat$DrugName))){
  d<-unique(ldat$DrugName)[i]
  dft<-ldat[ldat$DrugName==d,]
  
  for (j in 1:nrow(dft)){
    g<-dft$GeneTarget[j]
    if(dft$Num[j]!=0){
      pstr<-as.vector(str_split(dft$PathwayName[j],";",simplify = T))
      psid<-as.vector(str_split(dft$PathwayID[j],";",simplify = T))
      #  pstr<-str_match(pstr,"(.+) - homo")[,2]
      for (z in 1:length(pstr)){
        DF<-rbind.data.frame(DF,c(d,g,pstr[z],psid[z]))
      }
    }
    else {
      DF<-rbind.data.frame(DF,c(d,g,NA,NA))
    }
    
  }
  
}

#Add the originating lists
DF<-DF[-1,]
DF$Freq<-1

DFF<-cbind.data.frame(DF,List=NA)
DFF<-DFF[1,]
for (i in 1:nrow(DF)){
  d<-DF$Drug[i]
  l<-ldat2$Within.Top.20.List.ID[which(ldat2$Drug.Name==d)]
  l<-as.vector(str_split(l,",",simplify = T))
  for (j in 1:length(l)){
    DFF<-rbind.data.frame(DFF,c(DF$Drug[i],DF$Gene[i],DF$Pathway[i],DF$PathwayID[i],DF$Freq[i],l[j]))
  }
}
DFF<-DFF[-1,]

DF<-DFF

#Keep the Shortlisted drugs
DF3<-DF[which(!is.na(match(DF$Drug,sdr$Drug.Name[sdr$Shortlisted>0]))),]


DF3$Drug2<-DF3$Drug
DF3$do<-0
for (i in 1:nrow(DF3)){
  DF3$do[i]<-sdr$Max.Norm.Ranking[match(DF3$Drug2[i],sdr$Drug.Name)]
}

DF3[which(DF3$Demographic=="Drug"),]<-DF3[which(DF3$Demographic=="Drug"),][order(1/DF3$do[which(DF3$Demographic=="Drug")],decreasing = T),]
DF3<-DF3[order(DF3$do,decreasing = T),]

#Map the KEGG pathways to IDs
nods$pathid<-NA
nods$pathid<-kp$ID[match(nods$Name,kp$Name)]

#extract and order lists,genes,drugs and pathways
list<-as.character(unique(DF3$List))
list<-list[order(list)]
gene<-as.character(unique(DF3$Gene))
gene<-gene[order(gene)]
drug<-as.character(unique(DF3$Drug))
drug<-drug[order(sdr$Max.Norm.Ranking[match(drug,sdr$Drug.Name)],decreasing = T)]
path<-as.character(unique(DF3$Pathway))
path<-path[order(nods$Size[match(path,tolower(as.character(nods$Name)))],decreasing = T)]

DF3<-DF3[which(!is.na(match(DF3$PathwayID,nods$pathid[nods$Size>1.5]))),]

bb <- to_lodes_form(DF3, key = "Demographic", axes = c(6,1,2,3))
# calculate within-group proportions
bb <- na.omit(bb)
bb <- group_by(bb, Demographic)
bb <- add_count(bb, Freq, name = "Total")
bb <- ungroup(bb)
bb$Freq<-as.numeric(bb$Freq)
bb <- transform(bb, Prop = Freq / Total)

# alluvial diagram of within-axis proportions
bb$stratum<-ordered(bb$stratum,levels=c(list,drug,gene,path))

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

for (i in unique(bb$stratum[bb$Demographic=="Gene"])){
  ud<-unique(bb$stratum[which(!is.na(match(bb$alluvium,bb$alluvium[bb$stratum==i])) & bb$Demographic=="Drug")])
  if(length(ud)==1){
    c<-c+1
    bb$col[which(bb$Demographic=="Gene" & bb$stratum==i)]<-as.character(ud)
  }
  else{bb$col[which(bb$Demographic=="Gene" & bb$stratum==i)]<-"A"}
}

for (i in unique(bb$stratum[bb$Demographic=="Pathway"])){
  ud<-unique(bb$stratum[which(!is.na(match(bb$alluvium,bb$alluvium[bb$stratum==i])) & bb$Demographic=="Drug")])
  if(length(ud)==1){
    c<-c+1
    bb$col[which(bb$Demographic=="Pathway" & bb$stratum==i)]<-as.character(ud)
  }
  else{bb$col[which(bb$Demographic=="Pathway" & bb$stratum==i)]<-"A"}
}


ggplot(bb, aes(x = Demographic, stratum = stratum, alluvium = alluvium, y = Prop, label = stratum, fill=col)) +
  geom_flow(aes(fill = Drug2)) +
  geom_stratum(aes(color="darkgrey")) +
  geom_text(stat = "stratum") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())+
  guides(fill = FALSE) +
  scale_fill_brewer(type = "qual", palette = "Paired")+
  scale_color_brewer(type = "qual", palette = "Paired")# +

