####################################################################################
######  DEPARTMENT OF BIOINFORMATICS                                                
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By Marios Tomazou
######  email: mariost@cing.ac.cy
######
######  This script takes as an input the resulting files from the integration step 
######  comprising lists of ranked genes, their MIGe, MIGn and MIG scores and       
######  creates a network with binned node sizes and coloured based on their        
######  originating dataset. In addition and calculates the degree, Betweeness and  
######  Closeness and plots their distributions.
######
######  Developed in R Studio R version 3.6.1 (2019-07-05)
####################################################################################

library(igraph)
library(readxl)
library(ggplot2)
library(RColorBrewer)

dat1<-read.csv("IGGne_genelist_final_codres_ordered.tsv",sep = "\t",header = T,stringsAsFactors = F)
dat2<-read.csv("geneEdgelist_ids_bonus3_genesys300.tsv",sep = "\t",header = F,stringsAsFactors = F)
dat3<-read.csv("geneEdgelist_ids_bonus3.tsv",sep = "\t",header = F,stringsAsFactors = F)
con<-read.csv("GM_tableI3f.txt",header = T,sep="\t")

#filter out very week edges (optional) to simplify the visualisation
dat3<-dat3[dat3$V3>0.0005,]

#Create graph
net<-graph_from_data_frame(dat3,directed = F)
E(net)$weight<-dat3$V3

#Normalise, rank and bin node sizes
s<-dat1$IGGne05[which(!is.na(match(dat1$GMid,unique(c(dat3$V1,dat3$V2)))))]
s1<-s/max(s)
ranks<-rank(s1)
q<-quantile(s)
qs<-s
qs[s<q[2]]<-0.3
qs[s>=q[2]]<-0.6
qs[s>=q[3]]<-0.9
qs[ranks>=(max(ranks)-300)]<-2

#Annotate nodes
nods<-data.frame(name=dat1$GMid[which(!is.na(match(dat1$GMid,unique(c(dat3$V1,dat3$V2)))))],Score=s,Ns=s1,rank=ranks,Q=qs)
nods$Color<-NA
nods$Color[nods$Q==2]<-"royalblue"
nods$stroke<-0
nods$stroke[nods$Q==2]<-4

con<-con[which(!is.na(match(con$GMid,unique(c(dat3$V1,dat3$V2))))),]


ls<-list()
sv<-c()
colls<-list()
colvect<-c()

my_pal<-brewer.pal(5,"Paired")

for (i in 1:nrow(con)){
  vect<-as.numeric(con[i,c(7:11)])
  vect[vect>0]<-1
  vect<-vect/length(vect[vect>0])
  ls[[i]]<-vect
  if(sum(vect==0)==4){sv<-c(sv,"circle")
  colls[[i]]<-my_pal[which(vect!=0)]
  }
  if(sum(vect==0)<4){sv<-c(sv,"pie")
  colls[[i]]<-my_pal
  }
colvect<-c(colvect,my_pal[which(vect!=0)][1])
}


V(net)$shape<-sv

##to generate layouts for different nets
#lo<-layout_with_lgl(net,repulserad = 100)
#lo2<-(0.9*lo+0.1*layout.drl(net))
#write.csv(lo2,"layout_mix_drl_lgl2.csv",row.names = F,col.names = F)

#for a fixed layout
lo2<-as.matrix(read.csv("layout_mix_drl_lgl.csv",header = T))[c(1:length(V(net)$name)),]

#plot network
plot.igraph(net,vertex.label="",
            #vertex.shape="pie", 
            vertex.pie=ls,
            vertex.color=colvect,
            vertex.pie.color=list(my_pal),
            vertex.pie.border=nods$stroke/2,
            vertex.border.color=nods$Color,
            vertex.size=3*nods$Q,layout=lo2,
            vertex.frame.color=nods$Color,
            vertex.frame.size=0.1,
            edge.width=10*dat3$V3/max(dat3$V3))



#create full network and derive analytics
dat3<-read.csv("geneEdgelist_ids_bonus3.tsv",sep = "\t",header = F,stringsAsFactors = F)
net2<-graph_from_data_frame(dat3)
av<-dat1$GMid[which(is.na(match(dat1$GMid,unique(c(dat3$V1,dat3$V2)))))]
net2<-add_vertices(net2,nv=length(av),name=av)
df<-data.frame(Name=rep(dat1$GMmapSYS,6))
df$Attr<-rep(c("Edge Score","Nodal Score","Combined Score","Degree","Betweeness","Closeness"),each=nrow(dat1))
df$V<-c(dat1$IGGe,dat1$Inode3f,dat1$IGGne05,degree(net2),betweenness(net2),closeness(net2))

#plot the distributions
library(gridExtra)
c<-0
pl<-list()
for (i in unique(df$Attr)){
c<-c+1
p<-ggplot(df[df$Attr==i,])+geom_density(aes(x=V,y=..scaled..),fill="lightgreen")+ggtitle(i)
pl[[c]]<-p
}
grid.arrange(grobs=pl)

ggplot(df[c(df$Attr=="Edge Score" | df$Attr=="Nodal Score"),],aes(x=reorder(Name,c(dat1$IGGne05,dat1$IGGne05)),y=V/2,fill=Attr))+geom_col()+
  coord_flip()

