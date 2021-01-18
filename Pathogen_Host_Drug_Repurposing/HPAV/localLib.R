########################################################################################
######  DEPARTMENT OF BIOINFORMATICS 
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By George Minadakis 
######  email: georgem@cing.ac.cy
######
######  Developed in R Studio R version 3.6.1 (2019-07-05)
########################################################################################

getDrugScores<-function(drugIDs,targetDB,cProts,vProts,xProts,vProteome,xProteome, method="simple"){
  
  drugIDs<-trimws(drugIDs)
  drugNames<-matrix(nrow = 0,ncol = ncol(drugDB)+7)
  
  if(length(drugIDs)>0){
    idx<-which(toupper(targetDB[,"ID"]) %in% toupper(drugIDs))
    
    if(length(idx)>0){
      tmpTargetDB<-targetDB[idx,,drop=FALSE]
      
      # find which proteins each drug hits
      protCol<-c(); 
      nv_targets<-c();nx_targets<-c();
      nvp_targets<-c();nxp_targets<-c();
      nd_targets<-c(); ncom_targets<-c(); 
      #xp_scores<-c(); 
      
      nNewProts<-c()
      for(drug in drugIDs){
        idn<-which(tmpTargetDB[,"ID"] %in% drug)
        drugTargets<-toupper(tmpTargetDB[idn,"UniProt.ID"])
        
        ######## find which given protein hits ##########
        dtv <- which(drugTargets %in% vProts); v_targets<-drugTargets[dtv]
        dtx <- which(drugTargets %in% xProts); x_targets<-drugTargets[dtx]
        
        dtvp <- which(drugTargets %in% vProteome); vp_targets<-drugTargets[dtvp]
        dtxp <- which(drugTargets %in% xProteome); xp_targets<-drugTargets[dtxp]
        
        dtc <- which(drugTargets %in% cProts); c_targets<-drugTargets[dtc]
        
        nv_targets<-c(nv_targets,length(v_targets))
        nx_targets<-c(nx_targets,length(x_targets))
        nvp_targets<-c(nvp_targets,length(vp_targets))
        nxp_targets<-c(nxp_targets,length(xp_targets))
        ncom_targets<-c(ncom_targets,length(c_targets))
        
        nd_targets<-c(nd_targets,length(drugTargets))
        nNewProts<-c(nNewProts,length(xProts))
        
        all_targets<-c(v_targets,vp_targets)
        protCol<-c(protCol,paste(all_targets,collapse = ","))
        #print(paste(drug,"targets:",length(drugTargets),"src:",length(dtv),"dest:",length(dtx)))
      }
      
      # get drug score
      overlap<-paste(nv_targets,"/",length(vProts),sep = "")
    
      S0<-1/(10^(1/length(drugTargets)))
      
      S1<-nv_targets/length(vProts)
      S2<-nv_targets/nd_targets
      S3<-nx_targets/length(xProts)
      S4<-nx_targets/nd_targets
      
      S5<-nvp_targets/length(vProteome)
      S6<-nvp_targets/nd_targets
      S7<-nxp_targets/length(xProteome)
      S8<-nxp_targets/nd_targets
      
      #S6<- (xp_scores-min(xp_scores))/(max(xp_scores)-min(xp_scores))

      S9<-ncom_targets/length(cProts)
      
      if(method=="simple"){
        SDRUG<- S9
      }
      
      if(method=="4PH-fold"){
        SDRUG<- (S1 +  S3 + S5 + S7)/4
      }
     
      #print(paste(length(idx),length(overlap),length(nv_targets),length(SDRUG),length(protCol)))
      
      idx<-which(drugDB[,"ID"] %in% drugIDs)
      drugNames<-cbind(drugDB[idx,,drop=FALSE],overlap,nv_targets,nx_targets,nvp_targets,nxp_targets,SDRUG,protCol)
    }
  }
  colnames(drugNames)<-c(colnames(drugDB),"overlap","NVtar","NXtar","NVPtar","NXPtar","SDRUG","proteinIDs")
  drugNames<-drugNames[,-c(3,4,5),drop=FALSE]
  
  drugNames<-drugNames[order(as.double(drugNames[,c("SDRUG"),drop=FALSE]), decreasing=TRUE),,drop=FALSE]
  return(drugNames[,,drop=FALSE])
}

getSARSScanScores<-function(mynet,all_edge_drugs,nodeProteins,drugDB,sortby="WSDR"){
  
  Nedge<-nrow(mynet)
  NVirus<-length(nodeProteins)
  w1<-0.25; w2<-0.25; w3<-0.25;w4<-0.25
  uDrugs<-unique(all_edge_drugs[,"ID"])
  idx<-which(drugDB[,"ID"] %in% uDrugs)
  finalDrugs<-drugDB[idx,,drop=FALSE]
  
  ESDR<-c();PSDR<-c(); ESDR1<-c(); PSDR1<-c()
  appears<-c(); uTargets<-c(); uProts<-c()
  SS1<-c();SS2<-c();SS3<-c();SS4<-c();
  pathogens<-c()
  for(drug in finalDrugs[,"ID"]){
    idx <-which(all_edge_drugs[,"ID"] %in% drug)
    entries<-all_edge_drugs[idx,,drop=FALSE]
    
    Napp<-nrow(entries)
    prots<-unique(unlist(strsplit(entries[,"proteinIDs"],",")))
    Ntar<-length(prots)
    
    S1<-Napp/Nedge
    S2<- Ntar/NVirus
    S3<-sum(as.numeric(entries[,"SDRUG"])) / Napp
    S4<-sum(as.numeric(entries[,"WSEdge"])) / Napp
    
    #iSDR<-((S1 + S2 + S3)/3.0) * S4
    #iSDR<-S3*S4
    iSDR<-((S2 + S3)/2.0)* S4
    
    ESDR<-c(ESDR,iSDR)
    SS1<-c(SS1,S1);SS2<-c(SS2,S2);SS3<-c(SS3,S3);SS4<-c(SS4,S4)
    
    appears<-c(appears,Napp)
    uTargets<-c(uTargets,Ntar)
    uProts<-c(uProts,paste(prots,collapse = ","))
    pathogens<-c(pathogens,paste(unique(entries[,"node2"]),collapse = ",\r"))
    # print(paste(drug,unique(entries[,"name"])))
  }
  
  
  finalDrugs<-cbind(finalDrugs,appears,uTargets,ESDR,SS3,SS4,SS1,SS2,uProts)
  
  colnames(finalDrugs)<-c(colnames(drugDB),"Napp","UNtar","WSDR","SSDRUG","SSEDGE","Napp/Nedge","Utar/N1","proteinIDs")  
  finalDrugs<-finalDrugs[,-c(3,4,5),drop=FALSE]
  finalDrugs<-finalDrugs[order(as.numeric(finalDrugs[,c(sortby),drop=FALSE]), decreasing=TRUE),,drop=FALSE]
  rownames(finalDrugs)<-1:nrow(finalDrugs)
  
  return(as.matrix(finalDrugs[,c("ID","name","group",sortby),drop=FALSE]))
}

scanVirusNET<-function(myNet,drugDB,targetDB,virusTable,scanType="edges",method="simple"){
  
  print("Scanning Virus NET ...")
  all_edge_drugs<-c()
  for(i in 1:nrow(myNet)){
    print(paste(myNet[i,1],"<==>",myNet[i,2]))
    
    idx<-which(virusTable[,"Virus"] %in% myNet[i,1])
    vProts<-parseItem(virusTable[idx,"HProteome"],sep=",")
    vProteome<-parseItem(virusTable[idx,"Proteome"],sep=",")
    
    idx<-which(virusTable[,"Virus"] %in% myNet[i,2])
    xProts<-parseItem(virusTable[idx,"HProteome"],sep=",")
    xProteome<-parseItem(virusTable[idx,"Proteome"],sep=",")
   
    cProts<-parseItem(myNet[i,"V8"],sep=",")
    
    if(scanType=="COMMON"){
      all_prots<-c(cProts)
      idt<-which(toupper(targetDB[,"UniProt.ID"]) %in% cProts)
      drugIDs<-unique(targetDB[idt,"ID"])
    }
    
    ## create drug score matrix 
    drugMatrix<-getDrugScores(drugIDs,targetDB,cProts,vProts,xProts,vProteome,xProteome, method=method)
    #print(drugMatrix[,c("ID","Ntar","SDRUG")])
    
    ## attach net to drug matrix 
    if(nrow(drugMatrix)>0){
      inet<-myNet[i,,drop=FALSE]
      if(nrow(drugMatrix)>1){
        inet<-(rep.row(myNet[i,,drop=FALSE],nrow(drugMatrix)))
      }
      colnames(inet)<-c("node1","node2","edge","N1","N2","WSEdge","USEdge","commons","D")
      drugMatrix<-cbind(drugMatrix[,,drop=FALSE],inet[,c(1,2,3,4,5,6,7,9),drop=FALSE])
      
      all_edge_drugs<-rbind(all_edge_drugs,drugMatrix[,,drop=FALSE])
    }
  }
  
  ## attach net scores 
  return(all_edge_drugs[,,drop=FALSE])
}

scanSARSNET<-function(myNet,drugDB,targetDB,virusTable,sarsDrugs,method="simple"){
  
  print("Scanning NET ...")
  all_edge_drugs<-c()
  for(i in 1:nrow(myNet)){
    print(paste(myNet[i,1],"<==>",myNet[i,2]))
    
    idx<-which(virusTable[,"Virus"] %in% myNet[i,1])
    vProts<-parseItem(virusTable[idx,"HProteome"],sep=",")
    vProteome<-parseItem(virusTable[idx,"Proteome"],sep=",")
    
    idx<-which(virusTable[,"Virus"] %in% myNet[i,2])
    xProts<-parseItem(virusTable[idx,"HProteome"],sep=",")
    xProteome<-parseItem(virusTable[idx,"Proteome"],sep=",")
    
    cProts<-parseItem(myNet[i,"V8"],sep=",")
    
    ## get neighbour node drugs 
    idx<-which(sarsDrugs[,"Virus"] %in% myNet[i,2])
    drugIDs<-unique(parseItem(sarsDrugs[idx,"drugIDs"],sep=";"))
    
    ## create drug score matrix 
    drugMatrix<-getDrugScores(drugIDs,targetDB,cProts,vProts,xProts,vProteome,xProteome, method=method)
   
    ## attach net to drug matrix 
    if(nrow(drugMatrix)>0){
      inet<-myNet[i,,drop=FALSE]
      if(nrow(drugMatrix)>1){
        inet<-(rep.row(myNet[i,,drop=FALSE],nrow(drugMatrix)))
      }
      colnames(inet)<-c("node1","node2","edge","N1","N2","WSEdge","USEdge","commons","D")
      drugMatrix<-cbind(drugMatrix[,,drop=FALSE],inet[,c(1,2,3,4,5,6,7,9),drop=FALSE])
      
      all_edge_drugs<-rbind(all_edge_drugs,drugMatrix[,,drop=FALSE])
    }
  }
  
  ## attach net scores 
  return(all_edge_drugs[,,drop=FALSE])
}
getWeightedScore<-function(Ncom,N1,N2,D){
  
  S1<-(2*Ncom)/(N1+N2)
  S2 <- min(c(N1,N2))/max(c(N1,N2))
  S3<-1-D
  penalty<-1/10^(1/Ncom)
  #penalty<-Ncom/min(c(N1,N2))
  P<-S1*S2*S3*penalty
  return(P)
}

getUnweightedScore<-function(Ncom,N1,N2,D){
  S1<-2/(N1+N2)
  S2 <- min(c(N1,N2))/max(c(N1,N2))
  P<-S1*S2*(1-D) ;#* penalty
  return(P)
}

getEdgeScores<-function(srcNET,i,distance_table){
  
  D = 0;#runif(1,0.1,1)
  
  idx<-which(toupper(distance_table[,"Virus"]) %in% toupper(srcNET[i,2,drop=FALSE]))
  if(length(idx)>0){
    D<-as.numeric(distance_table[idx,"D"])
    if(is.na(D)){
      D=0.9999999999
    }
    if(is.infinite(D)){
      D=0.9999999999
    }
    if(D==1){
      D=0.9999999999
    }
  }
  print(paste(D,"---",srcNET[i,1],"---",srcNET[i,2]))
  
  weightedScore<-getWeightedScore(as.numeric(srcNET[i,3]),
                                  as.numeric(srcNET[i,4]),
                                  as.numeric(srcNET[i,5]),D)
  
  unweightedScore<-getUnweightedScore(as.numeric(srcNET[i,3]),
                                      as.numeric(srcNET[i,4]),
                                      as.numeric(srcNET[i,5]),D)
  
  return(c(weightedScore,1-D,unweightedScore))
}

attachScores<-function(srcNET,distance_table){
  N <-nrow(srcNET)
  scores1<-c()
  scores2<-c()
  scores3<-c()
  for(i in 1:N){
    res<-getEdgeScores(srcNET,i,distance_table)
    scores1<-c(scores1,res[1])
    scores2<-c(scores2,res[2])
    scores3<-c(scores3,res[3])
  }
  newNET<-cbind(srcNET[,c(1,2,3,4,5)],scores1,scores2,srcNET[,c(6)],scores3)
  colnames(newNET)<-paste("V",1:ncol(newNET),sep = "")
  
  ## remove edges with zero score 
  newNET<-newNET[as.numeric(newNET[,"V6"])>0,,drop=FALSE]
  
  ## sort the final network
  newNET<-newNET[order(as.numeric(newNET[,"V6"]),decreasing =TRUE),,drop=FALSE]
  return(newNET[,,drop=FALSE])
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}


protein2drug<-function(drugDB,targetDB,proteins){
  
  ## get unique drug targets 
  idx<-which(toupper(targetDB[,"UniProt.ID"]) %in% proteins)
  drugNames<-matrix(nrow = 0,ncol = ncol(drugDB)+4)
  
  if(length(idx)>0){ 
    tmpTargetDB<-targetDB[idx,,drop=FALSE]
    
    ## fing how many given proteins each drug hits 
    itable<-as.matrix(table(tmpTargetDB[,"ID"]))
    
    ## fing which proteins eadh drug hits 
    protCol<-c()
    dTarNum<-c()
    for(drug in rownames(itable)){
      idr<-which(tmpTargetDB[,"ID"] %in% drug)
      protCol<-c(protCol,paste(tmpTargetDB[idr,"UniProt.ID"],collapse = ","))
      
      idr<-which(targetDB[,"ID"] %in% drug)
      dTarNum<-c(dTarNum,length(idr))
    }
    
    # get drug score
    idx<-which(drugDB[,"ID"] %in% rownames(itable))
    overlap<-paste(itable[,1],"/",length(proteins),sep = "")
    hit_targets<-itable[,1]
    hit_score<-as.numeric(hit_targets)/length(proteins)
    
    drugNames<-cbind(drugDB[idx,,drop=FALSE],overlap,hit_targets,hit_score,protCol)
  }
  colnames(drugNames)<-c(colnames(drugDB),"overlap","Ntar","SDRUG","proteinIDs")
  drugNames<-drugNames[,-c(3,4,5),drop=FALSE]
  
  drugNames<-drugNames[order(as.double(drugNames[,c("SDRUG"),drop=FALSE]), decreasing=TRUE),,drop=FALSE]
  return(drugNames[,,drop=FALSE])
}

parseItem<-function(srcStr,sep=","){
  return(as.vector(trimws(unlist(strsplit(srcStr,sep,fixed = TRUE)))))
}




