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

##Modification to make it work! 
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

###Create a single PS.l with all the information from Foxes and Raccoon

###Load previously generated data from specific experiment

GetFox <- TRUE
GetRaccoon <- TRUE
#GetTest <- FALSE

Indrun <- FALSE ##Make figures and analysis for Fox and Raccoon datasets separated
SSU_18S <- TRUE ## Analysis just with 18S primers
SSU_LSU <- TRUE ## Analysis with 18S and 28S primers
SSU_LSU_COI <- TRUE ## Analysis with 18S, 28S and COI primers

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
seqcounts <- plyr::join(rawcountsF, rawcountsR, by= "Primer_name")
colnames(seqcounts)<- c("Primer_name", "Raw_counts_Fox", "Raw_counts_Racc")

seqcounts %>%
  rowwise() %>%
  mutate(Total_reads = sum(Raw_counts_Fox, Raw_counts_Racc)) -> seqcounts

rm(rawcountsF, rawcountsR, along)
 
####Data frames by different taxonomic rank count of ASVs
readNumByPhylum <- sumSeqByTax(PS.l = PS.l, rank = "phylum")
readNumByGenus <- sumSeqByTax(PS.l = PS.l, rank = "genus")
readNumByFamily <- sumSeqByTax(PS.l = PS.l, rank = "family")
readNumByOrder <- sumSeqByTax(PS.l = PS.l, rank = "order")
readNumBySpecies <- sumSeqByTax(PS.l = PS.l, rank = "species")

####All primers both datasets 
AbPhy <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByPhylum)) ### Start a loop: fro every element in the list ...
{ 
  phyla <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByPhylum[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    phyla[1,1] <- 0    ### Add a zero in the first column
    phyla[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    phyla <- as.data.frame((readNumByPhylum[[i]]))  ###Make a data frame with the data included in each element of the list 
    phyla[,2] <- rownames(phyla) ### And use the rownames as information of the second column 
  }
  
  phyla[,3] <- names(readNumByPhylum)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(phyla) <- c("ASV", "Phyla", "Primer_name") ### change the names for the columns 
  AbPhy <- rbind(AbPhy, phyla) ### Join all the "individual" data frames into the final data frame 
  
}   ### close loop

rm(phyla)

rownames(AbPhy) <- c(1:nrow(AbPhy)) ### change the rownames to consecutive numbers 
AbPhy <- data.frame(Primer_name = AbPhy$Primer_name, Phyla = AbPhy$Phyla, ASV = AbPhy$ASV) ###change the order of the columns
AbPhy$ASV <- as.numeric(AbPhy$ASV)

Total_ASV<- aggregate(AbPhy$ASV, by=list(AbPhy$Primer_name), FUN= sum)  ### Add a new variable that will contain the sum of all the sequencing reads by primer pair
colnames(Total_ASV) <- c("Primer_name", "Total_ASV")
AbPhy <- join(AbPhy, Total_ASV, by= "Primer_name")

Relative_abundance = AbPhy$ASV/AbPhy$Total_ASV ### create a vector with the result of the operation 

AbPhy[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbPhy)[5] <- "Relative_abundance" ### Change the name of the column 

AbPhy$Primer_name <- gsub(pattern = " ", replacement = "", x = AbPhy$Primer_name) ### check if the primer names have extra spaces

AbPhy$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbPhy$Primer_name)

#intersect(AbPhyF$Primer_name, primers$Primer_name)
#setdiff(AbPhyF$Primer_name, primers$Primer_name)
#setdiff(primers$Primer_name, AbPhyF$Primer_name)

#AbPhy <- merge(AbPhy, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

#AbPhy <- plyr::join(AbPhy, rawcounts, by= "Primer_name")

###What about genus level?
AbGen <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByGenus)) ### Start a loop: for every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByGenus[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((readNumByGenus[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(readNumByGenus)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("ASV", "Genus", "Primer_name") ### change the names for the columns 
  AbGen <- rbind(AbGen, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

rm(genus)

rownames(AbGen) <- c(1:nrow(AbGen)) ### change the rownames to consecutive numbers 
AbGen <- data.frame(Primer_name = AbGen$Primer_name, Genus = AbGen$Genus, ASV = AbGen$ASV) ###change the order of the columns

Total_ASV<- aggregate(AbGen$ASV, by=list(AbGen$Primer_name), FUN= sum)  ### Add a new variable that will contain the sum of all the sequencing reads by primer pair
colnames(Total_ASV) <- c("Primer_name", "Total_ASV")
AbGen <- join(AbGen, Total_ASV, by= "Primer_name")

Relative_abundance = AbGen$ASV/AbGen$Total_ASV ### create a vector with the result of the operation 

AbGen[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbGen)[5] <- "Relative_abundance" ### Change the name of the column 

AbGen$Primer_name <- gsub(pattern = " ", replacement = "", x = AbGen$Primer_name) ### check if the primer names have extra spaces

AbGen$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbGen$Primer_name)

#AbGen <- merge(AbGen, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

#AbGen <- plyr::join(AbGen, rawcounts, by= "Primer_name")

rm(Relative_abundance, Total_ASV)

###What are the most amplified taxa?
TopPhylum <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByPhylum)) ### Start a loop: for every element in the list ...
{ 
  top <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByPhylum[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    top[1,1] <- 0    ### Add a zero in the first column
    top[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    top <- as.data.frame(readNumByPhylum[[i]])  ###Make a data frame with the data included in each element of the list 
    colnames(top)<- "j"
    top <-subset(top, top$j == max(top))
    top[,2] <- rownames(top) ### And use the rownames as information of the second column 
  }
  
  top[,3] <- names(readNumByPhylum)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(top) <- c("ASV_Phylum", "Phylum", "Primer_name") ### change the names for the columns 
  TopPhylum <- rbind(TopPhylum, top) ### Join all the "individual" data frames into the final data frame 
}

TopOrder <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByOrder)) ### Start a loop: for every element in the list ...
{ 
  top <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByPhylum[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    top[1,1] <- 0    ### Add a zero in the first column
    top[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    top <- as.data.frame(readNumByOrder[[i]])  ###Make a data frame with the data included in each element of the list 
    colnames(top)<- "j"
    top <-subset(top, top$j == max(top))
    top[,2] <- rownames(top) ### And use the rownames as information of the second column 
  }
  
  top[,3] <- names(readNumByOrder)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(top) <- c("ASV_Order", "Order", "Primer_name") ### change the names for the columns 
  TopOrder <- rbind(TopOrder, top) ### Join all the "individual" data frames into the final data frame 
}

TopFam <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByFamily)) ### Start a loop: for every element in the list ...
{ 
  top <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByFamily[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    top[1,1] <- 0    ### Add a zero in the first column
    top[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    top <- as.data.frame(readNumByFamily[[i]])  ###Make a data frame with the data included in each element of the list 
    colnames(top)<- "j"
    top <-subset(top, top$j == max(top))
    top[,2] <- rownames(top) ### And use the rownames as information of the second column 
  }
  
  top[,3] <- names(readNumByFamily)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(top) <- c("ASV_Family", "Family", "Primer_name") ### change the names for the columns 
  TopFam <- rbind(TopFam, top) ### Join all the "individual" data frames into the final data frame 
}

TopGen <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByGenus)) ### Start a loop: for every element in the list ...
{ 
  top <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByGenus[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    top[1,1] <- 0    ### Add a zero in the first column
    top[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    top <- as.data.frame(readNumByGenus[[i]])  ###Make a data frame with the data included in each element of the list 
    colnames(top)<- "j"
    top <-subset(top, top$j == max(top))
    top[,2] <- rownames(top) ### And use the rownames as information of the second column 
  }
  
  top[,3] <- names(readNumByGenus)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(top) <- c("ASV_Genus", "Genus", "Primer_name") ### change the names for the columns 
  TopGen <- rbind(TopGen, top) ### Join all the "individual" data frames into the final data frame 
}

TopByRank <- join(TopPhylum, TopOrder, by= "Primer_name")
TopByRank <- join(TopFam, TopByRank, by= "Primer_name")
TopByRank <- join(TopByRank, TopGen, by= "Primer_name")

TopByRank<- TopByRank[c("Primer_name", "Phylum", "ASV_Phylum", "Order", "ASV_Order", "Family", "ASV_Family", "Genus", "ASV_Genus")]

rm(top, TopPhylum, TopOrder, TopFam, TopGen)
#write.csv(TopByRank, "~/AA_Primer_evaluation/Top_by_Rank_Primers.csv")

#join(TopByRank, primerInput, by= "Primer_name")

######Group specific analysis####
####By marker 

###Amplicon size vs raw counts
require(ggpubr)

ampsiz1<- ggplot(PrimTax, aes(x = Expected, y = num.reads), geom=c("point", "smooth")) +
  scale_x_continuous(name = "Expected amplicon size") +
  scale_y_log10(name = "log10 Total sequencing reads")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), aes(size= 25))+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20),legend.position = "none")+
  labs(tag = "A)")+
  stat_cor(label.x = 500, label.y = log10(3000000), aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  stat_regline_equation(label.x = 500, label.y = log10(2200000))

ampsiz2 <-ggplot(PrimTax, aes(x = Expected, y =num.reads, color= Gen), geom=c("point", "smooth")) +
  scale_x_continuous(name = "Expected amplicon size") +
  scale_y_log10(name = "log10 Total sequencing reads")+
  scale_color_manual(values = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"))+ 
  geom_jitter(shape=16, position=position_jitter(0.2), aes(size= 25))+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(2,"line")) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20))+
  labs(tag = "B)")

pdf(file = "~/AA_Primer_evaluation/Figures/Supplementary_1.pdf", width = 8, height = 10)
grid.arrange(ampsiz1, ampsiz2, nrow= 2, ncol= 1)
dev.off()

##Correlation with obtained amplicon size

ampsiz3<- ggplot(PrimTax, aes(x = Obtained, y = num.reads), geom=c("point", "smooth")) +
  scale_x_continuous(name = "Obtained amplicon size") +
  scale_y_log10(name = "log10 Total sequencing reads")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), aes(size= 25))+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20),legend.position = "none")+
  labs(tag = "A)")+
  stat_cor(label.x = 400, label.y = log10(3000000), aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  stat_regline_equation(label.x = 400, label.y = log10(2200000))

ampsiz4 <-ggplot(PrimTax, aes(x = Obtained, y =num.reads, color= Gen), geom=c("point", "smooth")) +
  scale_x_continuous(name = "Obtained amplicon size") +
  scale_y_log10(name = "log10 Total sequencing reads")+
  scale_color_manual(values = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"))+ 
  geom_jitter(shape=16, position=position_jitter(0.2), aes(size= 25))+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(2,"line")) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20))+
  labs(tag = "B)")

pdf(file = "~/AA_Primer_evaluation/Figures/Supplementary_1.2.pdf", width = 8, height = 10)
grid.arrange(ampsiz3, ampsiz4, nrow= 2, ncol= 1)
dev.off()

##Correlation between obtained and expected amplicon size

ampsiz5 <-ggplot(PrimTax, aes(x = Obtained, y = Expected), geom=c("point", "smooth")) +
  scale_x_continuous(name = "Obtained amplicon size") +
  scale_y_continuous(name = "Expected amplicon size")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), aes(size= 25))+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20),legend.position = "none")+
  labs(tag = "A)")+
  stat_cor(label.x = 400, label.y = 650, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  stat_regline_equation(label.x = 400, label.y = 600)

ampsiz6 <-ggplot(PrimTax, aes(x = Obtained, y =Expected, color= Gen), geom=c("point", "smooth")) +
  scale_x_continuous(name = "Obtained amplicon size") +
  scale_y_continuous(name = "Expected amplicon size")+
  scale_color_manual(values = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"))+ 
  geom_jitter(shape=16, position=position_jitter(0.2), aes(size= 25))+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(2,"line")) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme(text = element_text(size=20))+
  labs(tag = "B)")

pdf(file = "~/AA_Primer_evaluation/Figures/Supplementary_1.3.pdf", width = 8, height = 10)
grid.arrange(ampsiz5, ampsiz6, nrow= 2, ncol= 1)
dev.off()

###Redundancy for parasites, fungi and diet

Primtax.comb2 <- subset(PrimTax, Gen=="28S" | Gen== "18S" | Gen== "COI")

ps.numord.l.comb2 <- PS.l[c(Primtax.comb2$Primer_name)][order(Primtax.comb2$num.reads, decreasing=TRUE)]

####Parasites 18S+28S+COI
ps.numord.l.comb2[["ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R"]]<- NULL ### No parasites or fungi 

PM.para.comb2 <- lapply(ps.numord.l.comb2, function (x) {
  subset_taxa(x, phylum%in%c("Nematoda",
                             "Apicomplexa",
                             "Platyhelminthes"))
})

PG.para.comb2 <- lapply(PM.para.comb2, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})

#lapply(PG.para.comb2, function (x) {
#  table(tax_table(x)[,5])
#}) 


mat.para.comb2 <- lapply(PG.para.comb2, function (x) {
  as.matrix(tax_table(x))
}) 


P.comb2 <- lapply(PG.para.comb2, function (y) {
  make.names(tax_table(y)[, 5])
}) 

Parasites.comb2 <- data.frame() ###Create the data frame 

for (i in 1: length(P.comb2)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(P.comb2[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((P.comb2[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(P.comb2)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Genus", "Count", "Primer_name") ### change the names for the columns 
  Parasites.comb2 <- rbind(Parasites.comb2, genus) ### Join all the "individual" data frames into the final data frame 
  
} 

Parasites.comb2 <- unique(Parasites.comb2[c("Genus", "Primer_name")]) ##Take unique combinations 
Parasites.comb2 <- join(Parasites.comb2, Primtax.comb2) 

paraprimer <- ggplot(Parasites.comb2, aes(Primer_name, color=Gen, fill=Gen))+
  geom_bar()+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(name = "Number of Parasites\n genera amplified")+
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))+
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF")) +
  scale_x_discrete(name = "Primer combination")+
  theme(legend.position = "none", text = element_text(size=15))+
  labs(tag = "A)")

Parasites <- Parasites.comb2$Genus
Parasites <- as.data.frame(table(Parasites))
Parasites <- Parasites[order(-Parasites$Freq),]

###Add a bar plot (Not necessary now... maybe supplementary)

ggplot(Parasites, aes(reorder(Parasites, Freq), Freq))+
  geom_point()+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(name = "Number of primer pairs")+
  scale_x_discrete(name = "Parasite genera")

###Genus of parasites in 18S + 28S + COI

Count.para.comb2 <- data.frame() ###Create the data frame 

for (i in 1: length(P.comb2)) ### Start a loop: fro every element in the list ...
{ 
  cnt <- data.frame() #### make an individual data frame ...
  
  cnt <- as.data.frame(length(unique(P.comb2[[i]])))  ###Make a data frame with the data included in each element of the list 
  cnt[,2] <- names(P.comb2[i]) ### Take the names of every list and use them to fill column 2 as many times the logitude of the column 1
  colnames(cnt) <- c("Count_para", "Primer_name") ### change the names for the columns 
  Count.para.comb2 <- rbind(Count.para.comb2, cnt) ### Join all the "individual" data frames into the final data frame 
} 

####compair all the primers for multiple intersections

P.intersection <- sapply(seq_len(length(P.comb2)), function(x)
  sapply(seq_len(length(P.comb2)), function(y) length(intersect(unlist(P.comb2[x]), unlist(P.comb2[y]))))
)

rownames(P.intersection)<- names(P.comb2)
colnames(P.intersection)<- names(P.comb2)

###In order to add the annotations in good order, it is necessary to have the same order in the intersection matrix and in the annotation table

P.clust <- hclust(dist(P.intersection), method = "complete") ##Dendogram
#library(dendextend)
as.dendrogram(P.clust) %>%
       plot(horiz = TRUE)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))
P.col$Primer_name <- rownames(P.col)

P.col <- merge(P.col, Primtax.comb2, by="Primer_name", sort= F)

col_groups <- P.col %>%
  select("Primer_name", "Gen", "Region") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Primer_name

col_groups$Primer_name<- NULL

colour_groups <- list( Gen= c("18S"= "#440154FF", "28S"= "#21908CFF", "COI"= "#FDE725FF"))

paraheatmap <- pheatmap(P.intersection, 
         color = plasma(100),
         border_color = NA,
         annotation_col = col_groups, 
         #annotation_row = col_groups,
         annotation_colors = colour_groups,
         #cutree_rows = 2,
         #cutree_cols = 2,
         show_rownames = F,
         show_colnames = F,
         main= "Redundant genus of parasites amplified")

pdf(file = "~/AA_Primer_evaluation/Figures/Figure_2.pdf", width = 10, height = 8)
paraheatmap
dev.off()

####Cummulative curve for parasites
num.taxa.para <-sapply(c( "species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.para.comb2, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.para <- unlist(lapply(PM.para.comb2, function (x) 
  sum(otu_table(x))))

PrimTax.para <-as.data.frame(cbind(num.taxa.para,as.data.frame(num.reads.para)))
PrimTax.para$Primer_name<- rownames(PrimTax.para)
markers <- PrimTax[,c(8,16)]
PrimTax.para <- join(PrimTax.para, markers)
PrimTax.para <- PrimTax.para[order(-PrimTax.para$num.reads.para),] 

ps.numord.l.para <- PM.para.comb2[c(PrimTax.para$Primer_name)][order(as.vector(unlist(PrimTax.para$num.reads.para)), decreasing=TRUE)]


cum.tax.para <- sapply(c("species","genus", "family","order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.para, rank)
})

colnames(cum.tax.para) <- paste0("cum.", colnames(cum.tax.para))

cum.tax.para <- as.data.frame(cbind(1:nrow(cum.tax.para), cum.tax.para))
cum.tax.para[2:28,7] <- PrimTax.para$Gen ## add marker amplified 
colnames(cum.tax.para)[1] <- "Primer_name"
colnames(cum.tax.para)[7] <- "Gen"

###Fungi 18s + 28S + COI

PM.fungi.comb2 <- lapply(ps.numord.l.comb2, function (x) {
  subset_taxa(x, phylum%in%c("Ascomycota",
                             "Basidiomycota",
                             "Zygomycota", "Mucoromycota", "Zoopagomycota", "Blastocladiomycota", "Cryptomycota", "Microsporidia"))
})

PG.fungi.comb2 <- lapply(PM.fungi.comb2, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})

#lapply(PG.fungi.comb2, function (x) {
#  table(tax_table(x)[,5])
#}) 

mat.fungi.comb2 <- lapply(PG.fungi.comb2, function (x) {
  as.matrix(tax_table(x))
}) 

F.comb2 <- lapply(PG.fungi.comb2, function (y) {
  make.names(tax_table(y)[,5])
}) 

Fungi.comb2 <- data.frame() ###Create the data frame 

for (i in 1: length(F.comb2)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(F.comb2[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((F.comb2[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(F.comb2)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Genus", "Count", "Primer_name") ### change the names for the columns 
  Fungi.comb2 <- rbind(Fungi.comb2, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

Fungi.comb2 <- unique(Fungi.comb2[c("Genus", "Primer_name")]) ##Take unique combinations 
Fungi.comb2 <- join(Fungi.comb2, Primtax.comb2) 

fungiprimer <-ggplot(Fungi.comb2, aes(Primer_name, color=Gen, fill=Gen))+
  geom_bar()+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(name = "Number of Fungi\n genera amplified")+
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))+
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF")) +
  #scale_x_discrete(name = "Primer combination")
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(), text = element_text(size=15))+
  labs(tag = "B)")

#pdf(file = "~/AA_Primer_evaluation/Figures/Supplementary_2.pdf", width = 10, height = 10)
#grid.arrange(paraprimer, fungiprimer, nrow= 1, ncol= 2, widths= c(2,1))
#dev.off()

Fungi <- Fungi.comb2$Genus
Fungi <- as.data.frame(table(Fungi))
Fungi <- Fungi[order(-Fungi$Freq),]

### plot (To improve) 

ggplot(Fungi, aes(reorder(Fungi, Freq), Freq))+
  geom_point()+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(name = "Number of primer pairs")+
  scale_x_discrete(name = "Fungi genera")

Count.fungi.comb2 <- data.frame()

for (i in 1: length(F.comb2)) ### Start a loop: fro every element in the list ...
{ 
  cnt <- data.frame() #### make an individual data frame ...
  
  cnt <- as.data.frame(length(unique(F.comb2[[i]])))  ###Make a data frame with the data included in each element of the list 
  cnt[,2] <- names(F.comb2[i]) ### Take the names of every list and use them to fill column 2 as many times the logitude of the column 1
  colnames(cnt) <- c("Count_fungi", "Primer_name") ### change the names for the columns 
  Count.fungi.comb2 <- rbind(Count.fungi.comb2, cnt) ### Join all the "individual" data frames into the final data frame 
} 

####compair all the primers for multiple intersections

F.intersection <- sapply(seq_len(length(F.comb2)), function(x)
  sapply(seq_len(length(F.comb2)), function(y) length(intersect(unlist(F.comb2[x]), unlist(F.comb2[y]))))
)
rownames(F.intersection)<- names(F.comb2)
colnames(F.intersection)<- names(F.comb2)

F.clust <- hclust(dist(F.intersection), method = "complete") ##Dendogram
#library(dendextend)
as.dendrogram(F.clust) %>%
  plot(horiz = TRUE)

F.col <- cutree(tree = F.clust, k = 2)
F.col  <- data.frame(cluster = ifelse(test = F.col  == 1, yes = "cluster 1", no = "cluster 2"))
F.col$Primer_name <- rownames(F.col)

F.col <- merge(F.col, Primtax.comb2, by="Primer_name", sort= F)

col_groups <- F.col %>%
  select("Primer_name", "Gen", "Region") ##"Expected"

row.names(col_groups)<- col_groups$Primer_name

col_groups$Primer_name<- NULL

colour_groups <- list( Gen= c("18S"= "#440154FF", "28S"= "#21908CFF", "COI"= "#FDE725FF"))

fungiheatmap <- pheatmap(F.intersection, 
         color = plasma(100),
         border_color = NA,
         annotation_col = col_groups,
         #annotation_row = col_groups,
         annotation_colors = colour_groups,
         #cutree_rows = 2,
         #cutree_cols = 2,
         show_rownames = F,
         show_colnames = F,
         main= "Redundant genus of Fungi amplified")

pdf(file = "~/AA_Primer_evaluation/Figures/Figure_3.pdf", width = 10, height = 8)
fungiheatmap
dev.off()

####Cummulative curve for fungi
num.taxa.fungi <-sapply(c("species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.fungi.comb2, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.fungi <- unlist(lapply(PM.fungi.comb2, function (x) 
  sum(otu_table(x))))

PrimTax.fungi <-as.data.frame(cbind(num.taxa.fungi,as.data.frame(num.reads.fungi)))
PrimTax.fungi$Primer_name<- rownames(PrimTax.fungi)
markers <- PrimTax[,c(8,16)]
PrimTax.fungi <- join(PrimTax.fungi, markers)
PrimTax.fungi <- PrimTax.fungi[order(-PrimTax.fungi$num.reads.fungi),] 

ps.numord.l.fungi <- PM.fungi.comb2[c(PrimTax.fungi$Primer_name)][order(as.vector(unlist(PrimTax.fungi$num.reads.fungi)), decreasing=TRUE)]


cum.tax.fungi <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.fungi, rank)
})

colnames(cum.tax.fungi) <- paste0("cum.", colnames(cum.tax.fungi))


cum.tax.fungi <- as.data.frame(cbind(1:nrow(cum.tax.fungi), cum.tax.fungi))
cum.tax.fungi[2:28,7] <- PrimTax.fungi$Gen ## add marker amplified 
colnames(cum.tax.fungi)[1] <- "Primer_name"
colnames(cum.tax.fungi)[7] <- "Gen"

###Other eukaryotic elements ("Diet" and passing material)
##Working approach 
PM.diet.comb3 <- lapply(ps.numord.l, function (x) {
  subset_taxa(x, phylum%in%c("Arthropoda", "Chlorophyta", "Streptophyta", "Chordata", 
                             "Annelida", "Mollusca", "Bryozoa", "Picozoa", "Porifera",
                             "Rotifera", "Tardigrada"))
})

##Also working approach to select primers but it doesn't eliminate bacterial ASVs
###Take just eukariotes 
#PM.diet.comb3 <- lapply(PM.diet.comb3, function (x) {
#  subset_taxa(x, superkingdom%in%"Eukaryota")
#})

##Then select those phyla that are NOT fungi or parasite
#PM.diet.comb3 <- lapply(ps.numord.l, function (x) {
#  subset_taxa(x, !(phylum%in%c("Ascomycota", "Basidiomycota", "Zygomycota", 
#                             "Mucoromycota", "Zoopagomycota", "Blastocladiomycota", 
#                             "Cryptomycota", "Microsporidia", "Nematoda",
#                             "Apicomplexa", "Platyhelminthes")))
#})

###Not working to genus leve, removing NAs is making the tax_glom bug... Solution tax_glom to family :S 
##This will make the comparison a bit complicated :S... Keep finding a solution 
PG.diet.comb3 <- lapply(PM.diet.comb3, function (x) {
  tax_glom(x, "family", NArm = T)
})

lapply(PG.diet.comb3, function (x) {
  table(tax_table(x)[,4])
}) 

mat.diet.comb3 <- lapply(PG.diet.comb3, function (x) {
  as.matrix(tax_table(x))
}) 

D.comb3 <- lapply(PG.diet.comb3, function (y) {
  make.names(tax_table(y)[,4])
}) 

Diet.comb3 <- data.frame() ###Create the data frame 

for (i in 1: length(D.comb3)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(D.comb3[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((D.comb3[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(D.comb3)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Genus", "Count", "Primer_name") ### change the names for the columns 
  Diet.comb3 <- rbind(Diet.comb3, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

Diet.comb3 <- unique(Diet.comb3[c("Genus", "Primer_name")]) ##Take unique combinations ##Column Genus has family names 
Diet.comb3 <- join(Diet.comb3, PrimTax) 

dietprimer <-ggplot(Diet.comb3, aes(Primer_name, color=Gen, fill=Gen))+
  geom_bar()+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(name = "Number of non fungi or parasites\n families amplified")+
  scale_color_manual(values = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"))+
  scale_fill_manual(values = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14")) +
  scale_x_discrete(name = "Primer combination")+
  theme(text = element_text(size=15))+
  labs(tag = "C)")

pdf(file = "~/AA_Primer_evaluation/Figures/Supplementary_2.pdf", width = 10, height = 20)
grid.arrange(paraprimer, fungiprimer, dietprimer, widths= c(2,1), layout_matrix= rbind(c(1,2), c(3,3)))
dev.off()

Diet <- Diet.comb3$Genus
Diet <- as.data.frame(table(Diet))
Diet <- Diet[order(-Diet$Freq),]

### Not necessary plot (To improve in case it's necessary) 

ggplot(Diet, aes(reorder(Diet, Freq), Freq))+
  geom_point()+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(name = "Number of primer pairs")+
  scale_x_discrete(name = "Not fungi or parasite families")

####compair all the primers for multiple intersections

D.intersection <- sapply(seq_len(length(D.comb3)), function(x)
  sapply(seq_len(length(D.comb3)), function(y) length(intersect(unlist(D.comb3[x]), unlist(D.comb3[y])))))

rownames(D.intersection)<- names(D.comb3)
colnames(D.intersection)<- names(D.comb3)

D.clust <- hclust(dist(D.intersection), method = "complete") ##Dendogram
#library(dendextend)
as.dendrogram(D.clust) %>%
  plot(horiz = TRUE)

D.col <- cutree(tree = D.clust, k = 2)
D.col  <- data.frame(cluster = ifelse(test = D.col  == 1, yes = "cluster 1", no = "cluster 2"))
D.col$Primer_name <- rownames(D.col)

D.col <- merge(D.col, PrimTax, by="Primer_name", sort= F)

col_groups <- D.col %>%
  select("Primer_name", "Gen", "Region") ##"Expected"

row.names(col_groups)<- col_groups$Primer_name

col_groups$Primer_name<- NULL

colour_groups <- list( Gen= c("18S"= "#440154FF", 
                              "28S"= "#21908CFF", 
                              "COI"= "#FDE725FF", 
                              "16S"= "pink",
                              "12S"= "#E3DAC9",
                              "ITS"= "#C46210",
                              "rbcL" = "#D0FF14"))

dietheatmap <- pheatmap(D.intersection, 
                         color = plasma(100),
                         border_color = NA,
                         annotation_col = col_groups,
                         #annotation_row = col_groups,
                         annotation_colors = colour_groups,
                         #cutree_rows = 2,
                         #cutree_cols = 2,
                         show_rownames = F,
                         show_colnames = F,
                         main= "Redundant families of other eukaryotes amplified")

pdf(file = "~/AA_Primer_evaluation/Figures/Figure_4.pdf", width = 10, height = 8)
dietheatmap
dev.off()

###Cummulative plot Diet 

num.taxa.diet <-sapply(c("species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.diet.comb3, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.diet <- unlist(lapply(PM.diet.comb3, function (x) 
  sum(otu_table(x))))

PrimTax.diet <-as.data.frame(cbind(num.taxa.diet,as.data.frame(num.reads.diet)))
PrimTax.diet$Primer_name<- rownames(PrimTax.diet)
markers <- PrimTax[,c(8,16)]
PrimTax.diet <- join(PrimTax.diet, markers)
PrimTax.diet <- PrimTax.diet[order(-PrimTax.diet$num.reads.diet),] 

ps.numord.l.diet <- PM.diet.comb3[c(PrimTax.diet$Primer_name)][order(as.vector(unlist(PrimTax.diet$num.reads.diet)), decreasing=TRUE)]


cum.tax.diet <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.diet, rank)
})

colnames(cum.tax.diet) <- paste0("cum.", colnames(cum.tax.diet))


cum.tax.diet <- as.data.frame(cbind(1:nrow(cum.tax.diet), cum.tax.diet))
cum.tax.diet[2:40,7] <- PrimTax.diet$Gen ## add marker amplified 
colnames(cum.tax.diet)[1] <- "Primer_name"
colnames(cum.tax.diet)[7] <- "Gen"

