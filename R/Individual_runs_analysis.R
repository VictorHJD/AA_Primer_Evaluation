##Make figures and analysis for Fox and Raccoon datasets separated
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

sumSeqByTax <- function(PS.l, rank) {
  lapply(PS.l, function(x){ 
    counts <- data.frame(cbind(asvCount=colSums(otu_table(x)), tax_table(x)))
    counts$asvCount <- as.numeric(as.character(counts$asvCount))
    tapply(counts$asvCount, counts[, rank], sum)})
}

unique.taxa.cumsum <- function(ps.numord.l, rank){
  taxonNa <- lapply(ps.numord.l, function (x) 
    get_taxa_unique(x, taxonomic.rank = rank))
  names(taxonNa) <- names(ps.numord.l)
  unique.taxonNa.cumsum <- sapply(0:length(taxonNa),
                                  function (i) length(unique(unlist(taxonNa[0:i]))))
  return(unique.taxonNa.cumsum)
}

if(!exists("primerInput")){
  source("~/GitProjects/AA_Primer_Evaluation/R/General_Primer_information.R") ##   
}

if(!exists("PS.l.f")){
  source("~/GitProjects/AA_Primer_Evaluation/R/Data_preparation.R") ##   
}

if(!exists("PS.l.r")){
  source("~/GitProjects/AA_Primer_Evaluation/R/Data_preparation.R") ##   
}

####Data frames by different taxonomic rank count of ASVs

readNumByPhylumF <- sumSeqByTax(PS.l = PS.l.f, rank = "phylum")
readNumByPhylumR <- sumSeqByTax(PS.l = PS.l.r, rank = "phylum")

AbPhyF <- data.frame() ###Create the data frame 
  
for (i in 1: length(readNumByPhylumF)) ### Start a loop: fro every element in the list ...
{ 
  phyla <- data.frame() #### make an individual data frame ...
    
  if(nrow(as.data.frame(readNumByPhylumF[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    phyla[1,1] <- 0    ### Add a zero in the first column
    phyla[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    phyla <- as.data.frame((readNumByPhylumF[[i]]))  ###Make a data frame with the data included in each element of the list 
    phyla[,2] <- rownames(phyla) ### And use the rownames as information of the second column 
  }
    
  phyla[,3] <- names(readNumByPhylumF)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(phyla) <- c("ASV", "Phyla", "Primer_name") ### change the names for the columns 
  AbPhyF <- rbind(AbPhyF, phyla) ### Join all the "individual" data frames into the final data frame 
    
}  ### close loop
  
rownames(AbPhyF) <- c(1:nrow(AbPhyF)) ### change the rownames to consecutive numbers 
AbPhyF <- data.frame(Primer_name = AbPhyF$Primer_name, Phyla = AbPhyF$Phyla, ASV = AbPhyF$ASV) ###change the order of the columns
  
AbPhyF %>%
  group_by(Primer_name) %>% 
  mutate(Total_ASV = sum(ASV)) -> AbPhyF ### Add a new variable that will contain the sum of all the sequencing reads by primer pair
  
Relative_abundance = AbPhyF$ASV/AbPhyF$Total_ASV ### create a vector with the result of the operation 
  
AbPhyF[,5] <- Relative_abundance ### And put it in the same order in the colum 5
  
colnames(AbPhyF)[5] <- "Relative_abundance" ### Change the name of the column 
  
AbPhyF$Primer_name <- gsub(pattern = " ", replacement = "", x = AbPhyF$Primer_name) ### check if the primer names have extra spaces
  
AbPhyF$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbPhyF$Primer_name)

#intersect(AbPhyF$Primer_name, primers$Primer_name)
#setdiff(AbPhyF$Primer_name, primers$Primer_name)
#setdiff(primers$Primer_name, AbPhyF$Primer_name)
  
AbPhyF <- merge(AbPhyF, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 
  
AbPhyF <- plyr::join(AbPhyF, rawcountsF, by= "Primer_name")
  
AbPhyF$numberCycles <- rep(35,nrow(AbPhyF))


###
AbPhyR <- data.frame() ###Create the data frame 
  
  for (i in 1: length(readNumByPhylumR)) ### Start a loop: fro every element in the list ...
  { 
    phyla <- data.frame() #### make an individual data frame ...
    
    if(nrow(as.data.frame(readNumByPhylumR[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
    {
      phyla[1,1] <- 0    ### Add a zero in the first column
      phyla[1,2] <- "NA" ### And a NA in the second column.
    }else               ### For the rest of the elements: 
    {
      phyla <- as.data.frame((readNumByPhylumR[[i]]))  ###Make a data frame with the data included in each element of the list 
      phyla[,2] <- rownames(phyla) ### And use the rownames as information of the second column 
    }
    
    phyla[,3] <- names(readNumByPhylumR)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
    colnames(phyla) <- c("ASV", "Phyla", "Primer_name") ### change the names for the columns 
    AbPhyR <- rbind(AbPhyR, phyla) ### Join all the "individual" data frames into the final data frame 
    
}  ### close loop
  
rownames(AbPhyR) <- c(1:nrow(AbPhyR)) ### change the rownames to consecutive numbers 
AbPhyR <- data.frame(Primer_name = AbPhyR$Primer_name, Phyla = AbPhyR$Phyla, ASV = AbPhyR$ASV) ###change the order of the columns
  
AbPhyR %>%
    group_by(Primer_name) %>% 
    mutate(Total_ASV = sum(ASV)) -> AbPhyR ### Add a new variable that will contain the sum of all the sequencing reads by primer pair
  
  Relative_abundance = AbPhyR$ASV/AbPhyR$Total_ASV ### create a vector with the result of the operation 
  
  AbPhyR[,5] <- Relative_abundance ### And put it in the same order in the colum 5
  
  colnames(AbPhyR)[5] <- "Relative_abundance" ### Change the name of the column 
  
  AbPhyR$Primer_name <- gsub(pattern = " ", replacement = "", x = AbPhyR$Primer_name) ### check if the primer names have extra spaces
  
  AbPhyR$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbPhyR$Primer_name)
  
  #intersect(AbPhyF$Primer_name, primers$Primer_name)
  #setdiff(AbPhyF$Primer_name, primers$Primer_name)
  #setdiff(primers$Primer_name, AbPhyF$Primer_name)
  
  AbPhyR <- merge(AbPhyR, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 
  
  AbPhyR <- plyr::join(AbPhyR, rawcountsR, by= "Primer_name")
  
  AbPhyR$numberCycles <- rep(35,nrow(AbPhyR))


###Cumulative graph
###Fox
num.taxa.f <-sapply(c("species","genus", "family", "order", "class", "phylum", "superkingdom"), function (rank){
    lapply(PS.l.f, function (x) 
      length(get_taxa_unique(x, taxonomic.rank = rank)))
})
  
num.reads.f <- unlist(lapply(PS.l.f, function (x) 
  sum(otu_table(x))))
  
PrimTaxFox <- as.data.frame(cbind(num.taxa.f, num.reads.f))
  
rm(num.reads.f, num.taxa.f)
  
PrimTaxFox <- as.data.frame(apply(PrimTaxFox, 2, unlist))
PrimTaxFox[,9] <- row.names(PrimTaxFox)  
colnames(PrimTaxFox)[9] <- "Primer_name"
  
ps.numord.l <- PS.l.f[order(PrimTaxFox$num.reads, decreasing=TRUE)]
  
###To get the total unique taxa "identified" in the mix 
cum.tax.fox <- sapply(c("species","genus", "family", "class", "order", "phylum", "superkingdom"), function (rank){
    unique.taxa.cumsum(ps.numord.l, rank)
})
  
total.fox <- cum.tax.fox[42,1:7] ## save the totals 
PrimTaxFox <- rbind(PrimTaxFox, total.fox)
PrimTaxFox[42,8]<- NA
PrimTaxFox[42,9]<- "Total"
  
##Estimate taxonomic coverage by primer pair as total species detected by primer divided with the number of unique species identified by all primers
Taxonomic_coverage = (PrimTaxFox$species/PrimTaxFox[42,1])*100 ### create a vector with the result of the operation 
  
PrimTaxFox[,10] <- Taxonomic_coverage ### And put it in the same order in the colum 10
  
colnames(PrimTaxFox)[10] <- "Taxonomic_coverage" ### Change the name of the column
  
##Estimate "taxonomic resolution" by primer pair as total species detected divided by the number of total genus (this definition is not correct)
Taxonomic_resolution = (PrimTaxFox$species/PrimTaxFox$genus)*100
PrimTaxFox[,11] <- Taxonomic_resolution ### And put it in the same order in the colum 11
colnames(PrimTaxFox)[11] <- "Taxonomic_resolution" ### Change the name of the column
  
###Linear model
  
summary(lm(formula= Taxonomic_coverage~num.reads.f, data = PrimTaxFox))
  
ggplot(PrimTaxFox, aes(x=num.reads.f, y=Taxonomic_coverage, color=Taxonomic_resolution)) + 
    geom_point()+
    theme_classic()+
    #geom_hline(yintercept=50, linetype="dashed", color = "red")+
    scale_color_gradient(low="red", high="blue", name= "Taxonomic\nResolution")+
    xlab("Number of reads per primer pair")+
    ylab("Taxonomic coverage")+
    scale_x_log10()+
    ylim(0,20)+
    geom_smooth()
  
summary(lm(formula= Taxonomic_resolution~num.reads.f, data = PrimTaxFox))
  
ggplot(PrimTaxFox, aes(x=num.reads.f, y=Taxonomic_resolution, color=Taxonomic_coverage)) + 
    geom_point()+
    theme_classic()+
    #geom_hline(yintercept=50, linetype="dashed", color = "red")+
    scale_color_gradient(low="red", high="blue", name= "Taxonomic\nCoverage")+
    xlab("Number of reads per primer pair")+
    ylab("Taxonomic resolution")+
    scale_x_log10()+
    ylim(0,250)+
    geom_smooth()
  
summary(lm(formula= PrimTaxFox$Taxonomic_coverage~Taxonomic_resolution, data = PrimTaxFox))
  
ggplot(PrimTaxFox, aes(x=Taxonomic_coverage, y=Taxonomic_resolution, color=num.reads.f)) + 
    geom_point()+
    theme_classic()+
    #geom_hline(yintercept=50, linetype="dashed", color = "red")+
    scale_color_gradient(low="red", high="blue", name= "Number of\nreads")+
    xlab("Taxonomic coverage")+
    ylab("Taxonomic resolution")+
    ylim(0,250)+
    xlim(0,17)+
    geom_smooth(method = lm)
  
ggplot(PrimTaxFox, aes(x=reorder(Primer_name, -num.reads.f), y=num.reads.f))+
    geom_point()+ 
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 1))+
    xlab("Primer name")+
    ylab("Number of reads")
  #scale_y_log10()
  
colnames(cum.tax.fox) <- paste0("cum.", colnames(cum.tax.fox))
  
cum.tax.fox <- as.data.frame(cbind(1:nrow(cum.tax.fox), cum.tax.fox))
prnames <- names(ps.numord.l) ## get names from the primers
#cum.tax[2:42,1] <- prnames ## add primer names 
colnames(cum.tax.fox)[1] <- "Primer_name"
  
cum.plot.fox <- ggplot(cum.tax.fox) +
    geom_step(aes(Primer_name, cum.species), color="red") +
    geom_text(aes(20, 2500, label="Species"), color="red") +
    geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
    geom_text(aes(20, 2100, label="Genera"), color="#0D0887FF") +
    geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
    geom_text(aes(20, 850, label="Families"), color="#7E03A8FF") +
    geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
    geom_text(aes(20, 350, label="Orders"), color="#CC4678FF") +
    geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
    geom_text(aes(20, 150, label="Classes"), color="#F89441FF") +
    geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
    geom_text(aes(20, 53, label="Phyla"), color="#F0F921FF") +
    scale_y_log10("Cummulative count of taxa") +
    scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
    annotation_logticks(sides="l") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

##Raccoon
num.taxa.r <-sapply(c("species","genus", "family", "order", "class", "phylum", "superkingdom"), function (rank){
    lapply(PS.l.r, function (x) 
      length(get_taxa_unique(x, taxonomic.rank = rank)))
})
  
num.reads.r <- unlist(lapply(PS.l.r, function (x) 
    sum(otu_table(x))))
  
PrimTaxRac <- as.data.frame(cbind(num.taxa.r, num.reads.r))
  
rm(num.reads.r, num.taxa.r)

rm(MAF, MAR, PS.l.f, PS.l.r)
