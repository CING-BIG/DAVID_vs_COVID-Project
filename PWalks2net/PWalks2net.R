################################################################################
######  DEPARTMENT OF BIOINFORMATICS                                                
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By Marios Tomazou
######  email: mariost@cing.ac.cy
######
######  This script takes as an input the result of the Pathwalks analysis,      
######  a fully random walk results, and optionally a layout matrix. Then it     
######  performs an Odds ratio analysis to identify the most enriched pathways   
######  which are highlighted in the resulting pathway community network and     
######  bar chart. Finally the script saves the results in table format.         
######
######  Developed in R Studio R version 3.6.1 (2019-07-05)
################################################################################

library(igraph)
library(readxl)
library(ggplot2)
library(RColorBrewer)

#Load input files
rand<-read.csv("FullyRandomWalks_GCFIC1500D_rankedPathways.tsv",sep = "\t",header = F,stringsAsFactors = F)
dat1<-read.csv("IF_PW_top300_rankedEdges.tsv",sep = "\t",header = F,stringsAsFactors = F)
dat2<-read.csv("IF_PW_top300_rankedPathways.tsv",sep = "\t",header = F,stringsAsFactors = F)
dat3<-read.csv("IF_PW_top300_clusters.txt",header = F,stringsAsFactors = F,sep = "\n")
lo<-read.csv("layout copy.csv")
lo<-as.matrix(lo[,c(2,3)])

DF<-merge.data.frame(rand,dat2,by="V1")
DF<-DF[,-3]
colnames(DF)<-c("Name","Topo","Map")

#Odds ratio analysis
DF$TopoF<-DF$Topo/sum(DF$Topo)
DF$MapF<-DF$Map/sum(DF$Map)
DF$TopoO<-(DF$Topo/sum(DF$Topo))/(1-DF$Topo/sum(DF$Topo))
DF$MapO<-(DF$Map/sum(DF$Map))/(1-DF$Map/sum(DF$Map))
DF$OR<-DF$MapF/DF$TopoF
DF$OR2<-DF$MapO/DF$TopoO
DF$x2<-NA
DF$pval<-NA
for (i in 1:nrow(DF)){
  mt<-matrix(ncol=2,nrow = 2)
  mt[1,1]<-DF$Topo[i]
  mt[2,1]<-sum(DF$Topo[-i])
  mt[1,2]<-DF$Map[i]
  mt[2,2]<-sum(DF$Map[-i])
  x2<-chisq.test(t(mt))
  DF$x2[i]<-as.vector(x2$statistic)
  DF$pval[i]<-as.vector(x2$p.value)
  }
DF$p_adj<-p.adjust(DF$pval,method="BH")
DF2<-DF[DF$OR>1,]


#Create the draft net
net<-graph_from_data_frame(dat1,directed = F)


#Extract and populate the network nodes data frame
nams<-V(net)$name
nods<-data.frame(Name=as.character(nams),Cluster=rep(NA,length(nams)),Color=rep(NA,length(nams)),Size=rep(NA,length(nams)))

c<-0
ncol<-length(grep("Cluster",dat3$V1))
my_pal<-brewer.pal(ncol,"Paired")
dat1$cl1<-NA
dat1$cl2<-NA
dat1$cross<-NA
for (i in 1:nrow(dat3)){
  if (length(grep("Cluster",dat3$V1[i]))>0){
    c<-c+1
  }
  if (length(grep("Cluster",dat3$V1[i]))==0){
  nods$Cluster[which(nods$Name==dat3$V1[i])]<-c
  nods$Color[which(nods$Name==dat3$V1[i])]<-my_pal[c]
  nods$Size[which(nods$Name==dat3$V1[i])]<-DF$OR2[which(DF$Name==dat3$V1[i])]
  }
}

#Refine colors
nods$color2<-nods$Color
nods$color2[nods$label==""]<-"#568774"
nods$color2[nods$color2=="#568774"& nods$Size<1]<-"#876856"

#Redine labels
nods$label<-rank(1/nods$Size)
nods$label[nods$label>30]<-""

 
#Calculate weighs according to Pathwalks clustering input
for (i in 1:nrow(dat1)){
dat1$cl1[i]<-nods$Cluster[which(nods$Name==dat1$V1[i])]
dat1$cl2[i]<-nods$Cluster[which(nods$Name==dat1$V2[i])]

if(dat1$cl1[i]==dat1$cl2[i]){
dat1$cross[i]<-F
}
else{dat1$cross[i]<-T}
}
weights <- ifelse(dat1$cross, 1, 200)

#Generate new layout if not provided externally
if(!exists("lo")){
lo<-layout.fruchterman.reingold(net, weights=weights)
}



#Bin node size
nods$Size2<-nods$Size
nods$Size2[nods$Size<1]<-0.5
nods$Size2[nods$Size>=1]<-1
nods$Size2[nods$Size>=2]<-1.5
nods$Size2[nods$Size>=3]<-2

#Cluster and refine membership according to Pathwalks defined clusters
c1<-cluster_louvain(net,weights = weights)
c1$membership<-nods$Cluster

#Plot Net
plot(c1,net,vertex.label=nods$label,vertex.size=5*(nods$Size2),vertex.color=nods$color2,layout=lo,edge.color = "darkgrey",mark.col=unique(nods$Color),vertex.label.color="black")

#Plot Bar Chart
DF<-DF[order(DF$OR2,decreasing = T),]
DF$Name2<-paste0(c(1:nrow(DF)),". ",DF$Name)
ggplot(DF[1:30,],aes(x=reorder(Name2,OR2),y=OR2,fill=log(DF$Map[1:30])/2))+geom_col()+coord_flip()+theme_bw()+ylab("Odds Ratio")+xlab("Pathway Name")+
  labs(fill = "log(Pathwalks Score)")

#Save OR results in csv format
write.csv(nods,"nods_top300.csv")
write.csv(DF2,"signpaths2_top300.csv")
