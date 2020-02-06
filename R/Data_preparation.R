###Create a single PS.l with all the information from Foxes and Raccoon

###Libraries
library(MultiAmplicon)
library(ggplot2)
library(data.table)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(plyr)
library(vegan)
library(gridExtra)
library(grid)
library(lattice)
library("pheatmap")
library("viridisLite")

##Data
if(!exists("primerInput")){
  source("~/GitProjects/AA_Primer_Evaluation/R/General_Primer_information.R") ##   
}

###Load previously generated data from specific experiment

GetFox <- TRUE
GetRaccoon <- TRUE

if(GetFox){
  ####Fox data 
  ###Multiamplicon data 
  MAF <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF3.Rds")
  
  ##Phyloseq list
  #PS.l.f <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Fox/PhyloSeqList.Rds")
  PS.l.f <- toPhyloseq(MAF, samples=colnames(MAF), multi2Single=FALSE)
  
  rawcountsF <- rowSums(getRawCounts(MAF))
  rawcountsF <- data.frame(rawcountsF)
  rawcountsF[,2] <- rownames(rawcountsF)
  colnames(rawcountsF) <- c("Raw_counts", "Primer_name")
  rownames(rawcountsF) <- c(1:nrow(rawcountsF))
  rawcountsF <- data.frame(Primer_name = rawcountsF$Primer_name, Raw_counts = rawcountsF$Raw_counts) ###change the order of the columns
  rawcountsF$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcountsF$Primer_name)
  rawcountsF$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcountsF$Primer_name)
}

if(GetRaccoon){
  ####Raccoon data
  ##Multiamplicon data 
  MAR <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Raccoon/MAR2.Rds")
  
  ##Phyloseq list
  #PS.l.r <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/PhyloSeqList_Raccall.Rds")
  PS.l.r <- toPhyloseq(MAR, samples=colnames(MAF), multi2Single=FALSE)
  
  #names(PS.l.r) <- names(STNCR)[keepr] ###Raccoon has 39 primers that generate information :/ Check that again! 
  ##read count by data set 
  rawcountsR <- rowSums(getRawCounts(MAR))
  rawcountsR <- data.frame(rawcountsR)
  rawcountsR[,2] <- rownames(rawcountsR)
  colnames(rawcountsR) <- c("Raw_counts", "Primer_name")
  rownames(rawcountsR) <- c(1:nrow(rawcountsR))
  rawcountsR <- data.frame(Primer_name = rawcountsR$Primer_name, Raw_counts = rawcountsR$Raw_counts) ###change the order of the columns
  rawcountsR$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcountsR$Primer_name)
  rawcountsR$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcountsR$Primer_name)
}

###Eliminate the empty primers first and then merge the phyloseq lists... empty list make the next function bug
PS.l.r[sapply(PS.l.r, is.null)]<- NULL
PS.l.f[sapply(PS.l.f, is.null)]<- NULL
##It is necessary to eliminate those that are empty in both 
along<- names(PS.l.r)
PS.l <- lapply(along, function(i) merge_phyloseq(PS.l.r[[i]], PS.l.f[[i]])) ##Merge all the information from both experiments
names(PS.l) <- names(PS.l.r) ###Use the names from racconn list 

##Change names in the phyloseq list for short names
originalnames<- as.data.frame(names(PS.l))
colnames(originalnames)<- "Primer_name"

primerInput%>%
  select(1,15) -> shortnames

join(originalnames, shortnames, by= "Primer_name") -> shortnames

names(PS.l)<- as.character(shortnames$Primer_comb_ID) 

rm(originalnames, shortnames)

#saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/MergedPhyloSeqList.Rds") ###Just in case the merging doesn't work again! 
#PS.l <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/MergedPhyloSeqList.Rds")

###Merged seq reads counts
#seqcounts <- plyr::join(rawcountsF, rawcountsR, by= "Primer_name")
#colnames(seqcounts)<- c("Primer_name", "Raw_counts_Fox", "Raw_counts_Racc")

#seqcounts %>%
#  rowwise() %>%
#  mutate(Total_reads = sum(Raw_counts_Fox, Raw_counts_Racc)) -> seqcounts

rm(rawcountsF, rawcountsR, along)
