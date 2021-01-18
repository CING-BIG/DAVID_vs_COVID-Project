#######################################################################################
######  DEPARTMENT OF BIOINFORMATICS                                                
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By Marios Tomazou
######  email: mariost@cing.ac.cy
######
######  The script takes as an input a list of pathogens with their NCBI taxonomy IDs 
######  and the DrugBank drugs targeting their proteins. Loose and exact match refer  
######  to drugs that target the species in which the taxonomy ID belogns or exactly  
######  the taxonomy ID, respectively. In addition, a second with the pathogens and   
######  their taxonomy IDs and the number of the corresponding human-pathogen protein 
######  interactions is also required. The output is a circular dendrogram based on the     
######  taxonomy distance highlighting pathogens of interest and showing as a heatmaps
######  the number of pathogen-human protein interactions and number of drugs         
######  targeting these pathogens.
######
######  Developed in R Studio R version 3.6.1 (2019-07-05)
#######################################################################################

library(ggplot2)
library(stringr)
library(readxl)
library(xml2)
library(igraph)
library(taxize)
library(microbiome)
library(phyloseq)
library(treeio)
library(ape)
library(network)
library(ggtree)
library(phylobase)
library(microbiomeViz)
library(ggnewscale)
library(miscHelpers)
library(viridis)


#########################Files prep###############

gmids<-read.csv("GM_Classification_Drugs.csv",header = T,stringsAsFactors = F)
dfpp<-read.csv("VT.csv",stringsAsFactors = F,header = T)

#########################Taxa prep###############

#Get NCBI tax ids
uids<-gmids$taxID[gmids$Human_host>0]

#Get the taxonomy and convert to tree
tl5<-classification(uids,db="ncbi")
tree5<-class2tree(tl5,varstep = T,check = F)

#Convert to phylo4 object
g1 <- as(tree5$phylo, "phylo4")

#Parse into a dataframe
ndf<-data.frame(taxid=rep(NA,length(as.vector(g1@label))),Name=rep(NA,length(as.vector(g1@label))),Rank=rep(NA,length(as.vector(g1@label))))
for (i in 1:length(tl5)){
  ndf$taxid[i]<-tl5[[i]][nrow(tl5[[i]]),3]
  ndf$Name[i]<-tl5[[i]][nrow(tl5[[i]]),1]
  ndf$Rank[i]<-tl5[[i]][nrow(tl5[[i]]),2]
}

#get uids from NCBI for all nodes 
tdf<-get_uid(as.vector(g1@label[(nrow(dfpp)+1):length(g1@label)]))
ls<-tdf[1:length(tdf)]

#Replace manually the NAs from the get_uid() function - For other datasets the tax IDs must be manually entered in the vector below
ls[which(is.na(ls))]<-as.character(c(10533,12080,36375,370129,284213))

#Reclassify and fill the ndf data/frame
tls<-classification(ls,db="ncbi")
c<-0
for(i in (nrow(dfpp)+1):length(as.vector(g1@label))){
  c<-c+1
  ndf$taxid[i]<-tls[[c]][nrow(tls[[c]]),3]
  ndf$Name[i]<-tls[[c]][nrow(tls[[c]]),1]
  ndf$Rank[i]<-tls[[c]][nrow(tls[[c]]),2]
}

#Get number of drugs targeting proteins of each taxon and number of human-pathogen protein interactions
gmd<-gmids[gmids$Human_host>0,]

ndf$NDrugs<-NA
ndf$NDrugse<-NA
ndf$class<-NA
ndf$NDrugs[1:nrow(gmd)]<-str_count(gmd$drugs_loose,"DB")
ndf$NDrugse[1:nrow(gmd)]<-str_count(gmd$drugs_exact,"DB")

ndf$NPROT<-NA
ndf$col<-"darkgrey"

ind<-match(ndf$taxid,dfpp$taxID)
ndf$NPROT<-dfpp$proteinNum[ind]



#prepare the heatmap data.frame
nd<-ndf[1:625,c(4,5,7)]
for (i in 1:ncol(nd)){
  nd[which(is.na(nd[,i])),i]<-0
  
  nd[,i+3]<-as.character(nd[,i])
  nd[,i+3]<-values2colors(log(nd[,i]+0.1),col.pal="RdYlBu")
  nd[,i]<-log(nd[,i]+2)
  
}
nd$node<-c(1:nrow(nd))
tndf<-nd[,c(2,3)]
rownames(tndf)<-tree5$phylo$tip.label
hdf<-tndf
hdf<-as.data.frame(hdf[,-1])
rownames(hdf)<-rownames(tndf)

#viruses of interest for highlighting
voi<-c("SARS-CoV-2"=192,"Hepacivirus C"=745,"Ebolavirus"=815,"HIV 1"=854,"HIV 2"=853,"Influenza A"=741)


#create the plot
px<-ggtree(g1, layout='circular',aes(fill=log(c(ndf$NDrugse,mean(ndf$NDrugs)))))+labs(fill = "Anti-pathogen\nDrugbank Drugs\nlog(N)")+
  geom_hilight(node=1+which(ndf$Name=="Revtraviricetes"), fill="darkgrey", alpha=.2,extend = 10)+
  geom_hilight(node=1+which(ndf$Name=="cellular organisms"), fill="darkgreen", alpha=.2,extend = 10)+
  geom_hilight(node=1+which(ndf$Name=="Pisoniviricetes"), fill="steelblue", alpha=.2,extend = 15)+
  geom_hilight(node=1+which(ndf$Name=="Betacoronavirus"), fill="steelblue", alpha=.2,extend = 10)+
  geom_hilight(node=1+which(ndf$Name=="Monodnaviria"), fill="Red", alpha=.2,extend = 10)+
  geom_hilight(node=1+which(ndf$Name=="Bamfordvirae"), fill="Orange", alpha=.2,extend = 10)+
  geom_hilight(node=1+which(ndf$Name=="Orthornavirae"), fill="SteelBlue", alpha=.2,extend = 20)+
  geom_hilight(node=1+which(ndf$Name=="Heunggongvirae"), fill="Purple", alpha=.2,extend = 10)+
  geom_tippoint(size=(1.5*log(ndf$NPROT[1:625])/100+1),fill="darkgrey",shape=1, alpha=0.65)+theme(legend.position="right")+
  geom_point2(aes(subset=(node==as.vector(voi)[1])), shape=21, size=4, fill='green')+
  geom_text2(label=names(voi)[1],aes(subset=(node==as.vector(voi)[1])),nudge_x = 5)+
  geom_point2(aes(subset=(node==as.vector(voi)[2])), shape=21, size=4, fill='#40AEE4')+
  geom_text2(label=names(voi)[2],aes(subset=(node==as.vector(voi)[2])),nudge_x = 5)+
  geom_point2(aes(subset=(node==as.vector(voi)[3])), shape=21, size=4, fill='#E8384F')+
  geom_text2(label=names(voi)[3],aes(subset=(node==as.vector(voi)[3])),nudge_x = 5)+
  geom_point2(aes(subset=(node==as.vector(voi)[4])), shape=21, size=4, fill='#F59C49')+
  geom_text2(label=names(voi)[4],aes(subset=(node==as.vector(voi)[4])),nudge_x = 5)+
  geom_point2(aes(subset=(node==as.vector(voi)[5])), shape=21, size=4, fill='#FFE517')+
  geom_text2(label=names(voi)[5],aes(subset=(node==as.vector(voi)[5])),nudge_x = 5)+
  geom_point2(aes(subset=(node==as.vector(voi)[6])), shape=21, size=4, fill='#254597')+
  geom_text2(label=names(voi)[6],aes(subset=(node==as.vector(voi)[6])),nudge_x = 5)+
  geom_cladelabel(node=c(1,5,7,9), label="Riboviria", offset = 22,offset.text = 5) +
  geom_rootpoint(color="black", size=3)+
  scale_fill_gradient2(midpoint = mean(range(nd$NDrugse)), low = 'orange', high = 'red',mid = "orangered")
  
  ggplot(nd,aes(x=NPROT,y=NDrugse,fill=NDrugse))+geom_col()+




#add drug heatmap annotation
for (i in 1:nrow(nd)){

    if (length(descendants(g1,node=i))>=1){
    ind<-which(ndf$Name==ndf$Name[i])
    for(j in 1:length(ind)){
      if(isTRUE(ndf$Rank[ind[j]]!="superkingdom"|ndf$Rank[ind[j]]!="kingdom"|ndf$Rank[ind[j]]!="Phylum")){
        px<-px+geom_hilight(node=ind[j], fill=values2colors(col.start = "orange",col.end="red",nd$NDrugse)[i], alpha=.7,extend = log(ndf$NDrugse[i]))
      }
    }
  }
}

px<-px + new_scale_fill()


#add human-pathogen protein interactions heatmap
gheatmap(px, hdf, offset=15, width=.1,color = NA,
         colnames_offset_y = .25) +
  scale_fill_viridis_c(option="D", name="Human Protein\nInteractions\nlog2(N)")
