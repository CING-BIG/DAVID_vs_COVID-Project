####################################################################################
######  DEPARTMENT OF BIOINFORMATICS                                                
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By Marios Tomazou
######  email: mariost@cing.ac.cy
######
######  This script takes as an input a list of pathogens and the targets,
######  polypeptides and drugs files from Drugbank. The output is a data frame of
######  anti-pathogen drugs scored with respect to the taxonomic distance of 
######  the target pathogens to the pathogen of interest (herein SARS-CoV-2).
######
######  Developed in R Studio R version 3.6.1 (2019-07-05)
####################################################################################


library(dbparser)
library(ggplot2)
library(stringr)
library(readxl)
library(xml2)
library(igraph)
library(taxize)
library(microbiome)
library(MASS)


#########################Files prep###############
tars<-read.csv("drug_targ.csv",header = T,stringsAsFactors = F)
tpol<-read.csv("drug_targ_polys.csv",header = T,stringsAsFactors = F)
inf<-read.csv("drugs.csv",header = T,stringsAsFactors = F)
gmids<-read.csv("pathogens_class.csv",header = T,stringsAsFactors = F)

#SARS-CoV-2 Tax ID
scv2id<-"2697049"

#Map drugs and their protein targets based on DrugBank's Targets and Polypeptides files
tars<-tars[which(!is.na(match(tars$organism,gmids$Name))),]
tpol<-tpol[which(!is.na(match(tpol$organism_ncbi_taxonomy_id,gmids$taxID))),]
tpol$drugs<-NA
for (i in 1:nrow(tpol)){tpol$drugs[i]<-paste0(tars$parent_key[tars$id==tpol$parent_id[i]],collapse = ";")}

#########################Taxa prep###############

#Get Classification from NCBI and convert into a tree
tl<-classification(c(scv2id,unique(tpol$organism_ncbi_taxonomy_id)),db="ncbi")
tree2<-class2tree(tl,varstep = T,check = F)
#Get taxonomy distance matrix
mat<-data.frame(as.matrix(tree2[["distmat"]]))
colnames(mat)<-names(tl)
rownames(mat)<-names(tl)

#Get the tax IDs at the species level
tpol$species<-NA
tpol$speciestaxid<-NA
for (i in 1:nrow(tpol)){
ind<-match(tpol$organism_ncbi_taxonomy_id[i],names(tl))
tpol$species[i]<-tl[[ind]]$name[tl[[ind]]$rank=="species"]
tpol$speciestaxid[i]<-tl[[ind]]$id[tl[[ind]]$rank=="species"]
}

ddr<-data.frame(V1=rep("SARS-CoV-2",nrow(mat)),V1tax=rep(rownames(mat)[1],nrow(mat)),V2=rownames(mat),V2tax=rownames(mat),Weight=mat[,1])
ddr$V2<-as.character(ddr$V2)

#Generate the dataframe for all pathogens
ds<-log(1/(0.01+ddr$Weight))+abs(min(log(1/(0.01+ddr$tvarNCTR))))+0.01
ds2<-1/(0.01+ddr$tvarNCTR)
ddnet<-data.frame(V1=as.character(ddr$V1),Virus=as.character(ddr$V2),Dist=1/(0.01+ddr$Weight),dista=ds)
ddnet$drugs<-NA

dfd<-data.frame(Drug=NA,Name=NA,Weight=NA,stdev=NA,maxdev=NA,Richness=NA,Shannon_entropy=NA,Closest_Pathogen=NA,Most_distant_pathogen=NA,All_Pathogens=NA,Description=NA,Indication=NA,Interactions=NA,Toxicity=NA,Bioavailability=NA)
dmap<-tpol
for (i in 1:nrow(ddr)){
  ind<-NA
  id<-ddr$V2tax[i]
  ind<-match(id,dmap$organism_ncbi_taxonomy_id,nomatch = NA)
  if(!is.na(ind)){
    dlist<-c()
    we<-c()
    ddnet$drugs[i]<-dmap$drugs[ind]
    dlist<-as.vector(str_split(dmap$drugs[ind],";",simplify = T))
    dlist<-gsub(" ","",dlist)
    if (dlist==""){dlist<-NA}
    we<-rep(ddr$Weight[i],length(dlist))
    cpath<-rep(ddr$V2[i],length(dlist))
    if (!is.na(dlist)){
      nam<-c()
      desc<-c()
      indic<-c()
      interactions<-c()
      toxicity<-c()
      bioa<-c()
      stdv<-c()
      mdv<-c()
      taxnams<-c()
      taxids<-c()
      entropy<-c()
      richness<-c()
      diversemat<-c()
      allpaths<-c()
      distpath<-c()
      for(j in 1:length(dlist)){
     
        diversemat<-c()
        did<-dlist[j]
        
        
        taxnams<-unique(dmap$organism[grep(did,dmap$drugs)])
        taxids<-unique(dmap$organism_ncbi_taxonomy_id[grep(did,dmap$drugs)])
        
        taxidsmat<-table(dmap$organism_ncbi_taxonomy_id[grep(did,dmap$drugs)])

        clstax<-taxids[which.min(ddr$Weight[match(taxids,ddr$V2tax)])]
        dists<-mat[match(taxids,rownames(mat)),match(clstax,colnames(mat))]
        
        #Calculate the richness and shannon entropy for all tax IDs
        diversemat<-as.vector(table(dists))
        richness<-c(richness,as.vector(microbiome::richness((diversemat),index="observed")[1,1]))
        if (is.null(richness)){richness<-1}
        entropy<-c(entropy,microbiome::diversity(as.vector(diversemat),index = "shannon")[1,1])
        
        stdv1<-sd(dists)
        if(is.na(stdv1)){stdv1<-0}
        stdv<-c(stdv,stdv1)
        mdv<-c(mdv,max(dists))
        
        ind<-match(dlist[j],inf$primary_key)
        nam<-c(nam,inf$name[ind])
        desc<-c(desc,inf$description[ind])
        indic<-c(indic,inf$indication[ind])
        interactions<-c(interactions,inf$drug_interactions_count[ind])
        toxicity<-c(toxicity,inf$toxicity[ind])
        bioa<-c(bioa,inf$absorption[ind])
        
        allpaths<-c(allpaths,as.character(paste0(taxnams,collapse = ";")))
        distpath<-c(distpath,taxnams[which.max(dists)])


      }
    }
    dfdt<-data.frame(Drug=dlist,Name=nam,Weight=we,stdev=stdv,maxdev=mdv,Richness=richness,Shannon_entropy=entropy,Closest_Pathogen=cpath,Most_distant_pathogen=distpath,All_Pathogens=allpaths,Description=desc,Indication=indic,Interactions=interactions,Toxicity=toxicity,Bioavailability=bioa)
    dfd<-rbind.data.frame(dfd,dfdt)
  }
}

dfd<-dfd[-1,]

#Sort the data frame
dfd<-dfd[order(dfd$Weight,decreasing = F),]
dfd<-dfd[-which(duplicated(dfd$Drug)),]
dfd$Weight<-dfd$Weight/max(dfd$Weight)



             dfd$invd<-(1)/((dfd$Weight))
             dfd$invd_max_dev<-dfd$invd+dfd$maxdev/max(dfd$maxdev)
             dfd$max_dev_entr<-(dfd$Shannon_entropy)/(max(dfd$Shannon_entropy))
             dfd$Score<-dfd$invd+dfd$max_dev_entr
             
             dfd<-dfd[order(dfd$Score,decreasing = T),]
             
             write.csv(dfd,"Taxonomy_based_ranking.csv")







