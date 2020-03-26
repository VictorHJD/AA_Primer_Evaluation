####In silico evaluation of primer used in Fox and Raccoon analysis
##Load Libraries
library("doMC")
library("primerTree")
library("ggplot2")
library("reshape2")
library("plyr")
library("ape")
library("TmCalculator")
library("dplyr")

##
registerDoMC(15)

##Data
if(!exists("primerInput")){
  source("~/GitProjects/AA_Primer_Evaluation/R/General_Primer_information.R") ##   
}

if(!exists("seqcounts")){
  seqcounts<- read.csv(file="~/AA_Primer_evaluation/Reads_per_Primer_Pair.csv") ##   
}

if(!exists("seqrawcounts")){
  seqrawcounts<- read.csv(file="~/AA_Primer_evaluation/ASVs_per_Primer_Pair.csv") ##   
}

seqcounts %>%
  join(seqrawcounts, by= "Primer_comb_ID")%>%
  join(primerInput, by= "Primer_comb_ID")%>%->seqcounts

seqcounts[,1]<-NULL
seqcounts[,1]<-NULL

##PrimerTree pipeline
##Try with one primer pair

Euk_18S_01<- search_primer_pair(name='Euk_18S_01', 'CAGCAGCCGCGGTAATTCC', 'TCCGTCAATTTCTTTAAGTTTCAG', .parallel = T)
plot(Euk_18S_01)

