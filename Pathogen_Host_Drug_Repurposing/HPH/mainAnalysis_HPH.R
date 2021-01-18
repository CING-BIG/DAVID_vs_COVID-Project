########################################################################################
######  DEPARTMENT OF BIOINFORMATICS 
######  THE CYPRUS INSTITUTE OF NEUROLOGY & GENETICS 
######  By George Minadakis 
######  email: georgem@cing.ac.cy
######
######  This script takes as an input a drug information file from DrugBank, a list of 
######  pathogens and their NCBI taxonomy IDs, the taxonomy distance matrix of these
######  pathogens, a network of host pathogen interactions file and a list of host
######  proteins interacting with SARS-CoV-2. The output is the HPH list which include
######  drugs targeting the host proteins commonly affected by SARS-CoV-2 and other
######  pathogens
######
######  Developed in R Studio R version 3.6.1 (2019-07-05)
########################################################################################
library(data.table)
## clear everything first 
rm(list=ls(all=TRUE));  cat("\f"); graphics.off()

## import local library functions 
suppressWarnings(suppressMessages(source("localLib.r")))

## import the latest drug information available at the DrugBank repository 
drugDB<- as.matrix(fread("drugDB.txt"))

## import the latest update of drug targets  available at the DrugBank repository
targetDB<- as.matrix(fread("targetDB.txt"))

## Obtain pathogen information that you have collected. 
## The pathogens.txt file includes pathogens related to SARS-COV-2 network
## according to the analysis described in the manuscript
pathogens<- as.matrix(fread("pathogens.txt"))

## import SARS-COV-2 taxonomic distance table
sars_tax_distance_table<- as.matrix(fread("sars_tax_distance_table.txt"))

## load the SARS-COV-2 pathogen-2-pathogen network described in the manuscript 
myNet<- as.matrix(fread("sars_net.txt"))

## Create a vector of SARS-COV-2 host proteins of interest. 
## Herein we use the 336 host proteins described in the manuscript
sars2HProteins<- as.matrix(fread("sars2HProteins.txt"))

## attach edge weights on the network based on given equation 
myNet<-attachScores(myNet,sars_tax_distance_table)

## scan the pathogen-2-pathogen network using 4PH-fold equation to bring related drugs
all_edge_drugs<-scanVirusNET(myNet,drugDB,targetDB,pathogens,scanType = "COMMON",method="4PH-fold")

## get the above results and score the final drugs
finalDrugs<-getSARSScanScores(myNet,all_edge_drugs,sars2HProteins,drugDB,sortby = "WSDR")

## show the top 10 scored drugs
print(finalDrugs[1:10,])
