library(MultiAmplicon)
library(ggplot2)
library(data.table)
library(phyloseq)
library(dplyr)
library(vegan)
library(gridExtra)
library(grid)
library(lattice)
#library("ranacapa")
library("pheatmap")
library("viridisLite")

####Emanuel's function (don't work)
#sumSeqByTax2 <- function (Phy, tax) {
#  counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
#  counts$asvCount <- as.numeric(as.character(counts$asvCount))
#  tapply(counts$asvCount, counts[, tax], sum)
#}

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

#function for calculating relative amp efficiency for each taxon (Wrong! it can't be calculated for our results, at least not in this way)
##Mod from Kelly et al. 2019
EFFICIENCY <- function(Relative_abundance, Total_ASV, numberCycles, ...){
  temp <- unlist(
    (Total_ASV/Relative_abundance)^(1/numberCycles) - 1)
  temp[which(is.infinite(temp))] <- NA
  temp <- temp/(max(temp, na.rm=T))  #normalize to maximum, to create relative index of amp efficiency
  temp[temp<=0] <- 0
  as.numeric(temp)
}
###Create a single PS.l with all the information from Foxes and Raccoon

###Load previously generated data from specific experiment

GetFox <- TRUE

GetRaccoon <- FALSE

GetTest <- FALSE

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

##old pipeline
#MAF <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF_complete.RDS") ##Without taxa data
###No chimeras data 
#STNCF <- getSequenceTableNoChime(MAF)
###Annotation list 
#annot.list.f <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Fox/Fox_blast_tax_complete.Rds") 
##Discard primers without annotation 
#keepf <- unlist(lapply(annot.list.f, nrow))>0
#names(STNCF)[keepf]<-gsub(pattern = "-", replacement = "_", x= names(STNCF[keepf]))
#names(STNCF)[keepf]<-gsub(pattern = " ", replacement = "", x= names(STNCF[keepf]))
#names(PS.l.f) <- names(STNCF)[keepf]

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

##Old pipeline
#MAR <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Raccoon/MAR_complete.RDS") 
###No chimeras data 
#STNCR <- getSequenceTableNoChime(MAR)
###Annotation list 
#annot.list.r <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/Racc_blast_tax_all.Rds")
##Discard primers without annotation 
#keepr <- unlist(lapply(annot.list.r, nrow))>0
#names(STNCR)[keepr]<-gsub(pattern = "-", replacement = "_", x= names(STNCR[keepr]))
#names(STNCR)[keepr]<-gsub(pattern = " ", replacement = "", x= names(STNCR[keepr]))

###Merge both phyloseq lists ... Use as reference Raccoon list 
#(in theory the same primer pairs worked but "515F_99.R1073_103" and "16SA_L_110.16SB_H_110" didn't work for raccoons)

##Let's eliminate them to merge them correclty
##if both lists have the same primer should work properly 
#PS.l.f[["16SA_L_110_F.16SB_H_110_R"]]<- NULL
#PS.l.f[["515F_99_F.R1073_103_R"]]<- NULL

###Eliminate the empty primers first and then merge the phyloseq lists... empty list make the next function bug
PS.l.r[sapply(PS.l.r, is.null)]<- NULL
PS.l.f[sapply(PS.l.f, is.null)]<- NULL
##It is necessary to eliminate those that are empty in both 
along<- names(PS.l.r)
PS.l <- lapply(along, function(i) merge_phyloseq(PS.l.r[[i]], PS.l.f[[i]])) ##Merge all the information from both experiments
names(PS.l) <- names(PS.l.r) ###Use the names from racconn list 
#saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_Fox/MergedPhyloSeqList.Rds") ###Just in case the merging doesn't work again! 

#PS.l <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Fox/MergedPhyloSeqList.Rds")

PS.l <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/")

#PS.l[sapply(PS.l, is.null)]<- NULL

###Merged rawcounts
rawcounts <- plyr::join(rawcountsF, rawcountsR, by= "Primer_name")
colnames(rawcounts)<- c("Primer_name", "Raw_counts_Fox", "Raw_counts_Racc")

rawcounts %>%
  rowwise() %>%
  mutate(Total_reads = sum(Raw_counts_Fox, Raw_counts_Racc)) -> rawcounts

rm(rawcountsF, rawcountsR, along)
 
#n16S <- subset(primerInput, primerInput$Gen== "16S")
#n18S <- subset(primerInput, primerInput$Gen== "18S")
#n28S <- subset(primerInput, primerInput$Gen== "28S")
#COI <-  subset(primerInput, primerInput$Gen== "COI")
#plants <- subset(primerInput, primerInput$Target== "Viridiplantae")

####Data frames by different taxonomic rank count of ASVs
##readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum") Not working with new PS.l
readNumByPhylumF <- sumSeqByTax(PS.l = PS.l.f, rank = "phylum")
readNumByPhylumR <- sumSeqByTax(PS.l = PS.l.r, rank = "phylum")

readNumByPhylum <- sumSeqByTax(PS.l = PS.l, rank = "phylum")
readNumByGenus <- sumSeqByTax(PS.l = PS.l, rank = "genus")
readNumByFamily <- sumSeqByTax(PS.l = PS.l, rank = "family")
readNumByOrder <- sumSeqByTax(PS.l = PS.l, rank = "order")
readNumBySpecies <- sumSeqByTax(PS.l = PS.l, rank = "species")

#readNumByPhylumF<- lapply(getTaxonTable(MAF3), function (x) table(as.vector(x[, "phylum"])))
#readNumByPhylumR<- lapply(getTaxonTable(MAR2), function (x) table(as.vector(x[, "phylum"])))
# Then it is necessary to eliminate all empty primer pairs to make the rest of the code work (annoying)

#Function is not working anymore :( find out why
#readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") 
#readNumByFam <- lapply(PS.l, sumSeqByTax, "family")
#readNumbySpecies <-lapply(PS.l, sumSeqByTax, "species")

####
if(GetFox){
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
}

###
if(GetRaccoon){
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
}
####
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

rownames(AbPhy) <- c(1:nrow(AbPhy)) ### change the rownames to consecutive numbers 
AbPhy <- data.frame(Primer_name = AbPhy$Primer_name, Phyla = AbPhy$Phyla, ASV = AbPhy$ASV) ###change the order of the columns
AbPhy$ASV <- as.numeric(AbPhy$ASV)

AbPhy %>%
  group_by(Primer_name) %>% 
  mutate(Total_ASV = sum(ASV)) -> AbPhy ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

Relative_abundance = AbPhy$ASV/AbPhy$Total_ASV ### create a vector with the result of the operation 

AbPhy[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbPhy)[5] <- "Relative_abundance" ### Change the name of the column 

AbPhy$Primer_name <- gsub(pattern = " ", replacement = "", x = AbPhy$Primer_name) ### check if the primer names have extra spaces

AbPhy$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbPhy$Primer_name)


#intersect(AbPhyF$Primer_name, primers$Primer_name)
#setdiff(AbPhyF$Primer_name, primers$Primer_name)
#setdiff(primers$Primer_name, AbPhyF$Primer_name)

AbPhy <- merge(AbPhy, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

AbPhy <- plyr::join(AbPhy, rawcounts, by= "Primer_name")

#AbPhy$numberCycles <- rep(35,nrow(AbPhy))


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

rownames(AbGen) <- c(1:nrow(AbGen)) ### change the rownames to consecutive numbers 
AbGen <- data.frame(Primer_name = AbGen$Primer_name, Genus = AbGen$Genus, ASV = AbGen$ASV) ###change the order of the columns

AbGen %>%
  group_by(Primer_name) %>% 
  mutate(Total_ASV = sum(ASV)) -> AbGen ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

Relative_abundance = AbGen$ASV/AbGen$Total_ASV ### create a vector with the result of the operation 

AbGen[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbGen)[5] <- "Relative_abundance" ### Change the name of the column 

AbGen$Primer_name <- gsub(pattern = " ", replacement = "", x = AbGen$Primer_name) ### check if the primer names have extra spaces

AbGen$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbGen$Primer_name)


AbGen <- merge(AbGen, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

AbGen <- plyr::join(AbGen, rawcounts, by= "Primer_name")

#AbGen$numberCycles <- rep(35,nrow(AbGen))


####Estimate efficiency (function from Kelly et al. 2019) ###This is not possible with our data :(

#AbPhy$Efficieny <- EFFICIENCY(Relative_abundance = AbPhy$Relative_abundance, Total_ASV = AbPhy$Total_ASV, numberCycles = AbPhy$numberCycles)
#AbPhy$New_efficieny <- EFFICIENCY(Relative_abundance = AbPhy$Relative_abundance, Total_ASV = AbPhy$Total_reads, numberCycles = AbPhy$numberCycles)
#AbPhyF$Efficieny <- EFFICIENCY(Relative_abundance = AbPhyF$Relative_abundance, Total_ASV = AbPhyF$Total_ASV, numberCycles = AbPhyF$numberCycles)
#AbPhyR$Efficieny <- EFFICIENCY(Relative_abundance = AbPhyR$Relative_abundance, Total_ASV = AbPhyR$Total_ASV, numberCycles = AbPhyR$numberCycles)
#AbGen$Efficiency <- EFFICIENCY(Relative_abundance = AbGen$Relative_abundance, Total_ASV = AbGen$Total_ASV, numberCycles = AbGen$numberCycles)

###Plot distribution of efficiency by primer (By Phylum)

#require("ggridges")
#pdf(file = "~/AA_Primer_evaluation/Figures/Amplification_efficiency_FoxvsRaccoon.pdf", width = 10, height = 20)
#a <- ggplot(AbPhy, aes(x = Efficieny, y = reorder(Primer_name, Total_reads))) +
#      geom_density_ridges(aes(fill = Gen, alpha=0.5))+
#      geom_jitter(aes(alpha= 0.5))+
#      scale_x_continuous(limits = c(0,1))+
#      geom_vline(aes(xintercept=0.5),
#             linetype="dashed", size=1, colour="red")+
#      xlab("Amplification efficiency")+
#      ylab("Primer name")+
#      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#      labs(tag= "A)")

#b <- ggplot(AbPhyF, aes(x = Efficieny, y = reorder(Primer_name, Raw_counts))) +
#      geom_density_ridges(aes(fill = Gen, alpha=0.5))+
#      geom_jitter(aes(alpha= 0.5))+
#      scale_x_continuous(limits = c(0,1))+
#      geom_vline(aes(xintercept=0.5),
#             linetype="dashed", size=1, colour="red")+
#      xlab("Amplification efficiency")+
#      ylab("Primer name")+
#      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#      labs(tag = "B)")

#c <- ggplot(AbPhyR, aes(x = Efficieny, y = reorder(Primer_name, Raw_counts))) +
#      geom_density_ridges(aes(fill = Gen, alpha=0.5))+
#      geom_jitter(aes(alpha= 0.5))+
#      scale_x_continuous(limits = c(0,1))+
#      geom_vline(aes(xintercept=0.5),
#             linetype="dashed", size=1, colour="red")+
#      xlab("Amplification efficiency")+
#      ylab("Primer name")+
#      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#      labs(tag = "C)")

#grid.arrange(a, b, c, ncol= 1, nrow= 3)
#dev.off()

#ggplot(AbPhy, aes(x = New_efficieny, y = reorder(Primer_name, Total_reads))) +
#  geom_density_ridges(aes(fill = Gen, alpha=0.5))+
  #geom_jitter(aes(alpha= 0.5))+
#  scale_x_continuous(limits = c(0,1))+
#  geom_vline(aes(xintercept=.5),
#             linetype="dashed", size=1, colour="red")+
#  xlab("Amplification efficiency")+
#  ylab("Primer name")+
#  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#  labs(tag= "A)")

#d <- ggplot(AbGen, aes(x = Efficiency, y = reorder(Primer_name, Total_reads))) +
#      geom_density_ridges(aes(fill = Gen, alpha=0.5))+
#      geom_jitter(aes(alpha= 0.005))+
#      scale_x_continuous(limits = c(0,1))+
#      geom_vline(aes(xintercept=0.5),
#             linetype="dashed", size=1, colour="red")+
#      xlab("Amplification efficiency")+
#      ylab("Primer name")+
#      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#      labs(tag= "B)") #+ annotate("density", x = 1, y = "L2513_64_F.H2714_64_R", label = "OK")


#lapply(AbPhy, unZeroOne)
#AbPhy <- AbPhy[-c(373),]

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

#write.csv(TopByRank, "~/AA_Primer_evaluation/Top_by_Rank_Primers.csv")

#join(TopByRank, primerInput, by= "Primer_name")

###The lower efficiency, the higher amplification probability 

#to fit to a Beta, make strictly non-zero, non-one
#unZeroOne <- function(x){
#  x[x<1e-5] <- 1e-5  #set zero values to near zero
#  x[x> (1 - 1e-5)] <- (1 - 1e-5) #set one values to near one
#  x
#}


###Cumulative graph
###Fox

if(GetFox){
num.taxa.f <-sapply(c("species","genus", "family", "order", "class", "phylum", "superkingdom"), function (rank){
  lapply(PS.l.f, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.f <- unlist(lapply(PS.l.f, function (x) 
  sum(otu_table(x))))

PrimTaxFox <- as.data.frame(cbind(num.taxa.f, num.reads.f))

rm(num.reads.f)
rm(num.taxa.f)

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
}
##Raccoon
if(GetRaccoon){
num.taxa.r <-sapply(c("species","genus", "family", "order", "class", "phylum", "superkingdom"), function (rank){
  lapply(PS.l.r, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.r <- unlist(lapply(PS.l.r, function (x) 
  sum(otu_table(x))))

PrimTaxRac <- as.data.frame(cbind(num.taxa.r, num.reads.r))

rm(num.reads.r)
rm(num.taxa.r)
}

## All

num.taxa <-sapply(c("species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PS.l, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads <- unlist(lapply(PS.l, function (x) 
  sum(otu_table(x))))

PrimTax <- as.data.frame(cbind(num.taxa, num.reads))

PrimTax <- as.data.frame(apply(PrimTax, 2, unlist))
PrimTax[,8] <- row.names(PrimTax)  
colnames(PrimTax)[8] <- "Primer_name"

ps.numord.l <- PS.l[order(PrimTax$num.reads, decreasing=TRUE)]

cum.tax <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l, rank)
})

colnames(cum.tax) <- paste0("cum.", colnames(cum.tax))

cum.tax <- as.data.frame(cbind(1:nrow(cum.tax), cum.tax))
prnames <- names(ps.numord.l) ## get names from the primers
#cum.tax[2:42,1] <- prnames ## add primer names 
colnames(cum.tax)[1] <- "Primer_name"

cum.plot.all <- ggplot(cum.tax) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(20, 4800, label="Species"), color="red")+ 
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 2600, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 1200, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 500, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  #geom_text(aes(20, 200, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(20, 60, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa") +
  #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank())+
  labs(tag = "A)")


###richeness analysis 
###Calculate alpha diversity indexes for every primer pair 
#richness <- data.frame()

#for (i in 1:length(PS.l)){
  
#  a <- data.frame()
  
#  b <- data.frame(estimate_richness(PS.l[[i]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson")))
  
#  a<-colMeans(b)
  
#  richness<- rbind(richness, a)
#}

#richness[,7] <- names(PS.l)

#colnames(richness) <- c("Observed","Chao1","se_Chao1", "Shannon", "Simpson", "InvSimpson", "Primer_name")


### Obtain the number of samples 
#nsamp <- data.frame()

#nsamp[1:41,1] <-names(PS.l)

#for (i in 1:length(PS.l))
#  nsamp[i,2]<- nrow(otu_table(PS.l.f[[i]])) + nrow(otu_table(PS.l.r[[i]]))

#colnames(nsamp) <- c("Primer_name", "Number_samples")


###Check both elements and merge! 
#richness$Primer_name <- gsub(pattern = " ", replacement = "", x = richness$Primer_name) ### check if the primer names have extra spaces
#nsamp$Primer_name <- gsub(pattern = " ", replacement = "", x= nsamp$Primer_name)
#rawcounts <- as.data.frame(num.reads)
#rawcounts[,2] <- rownames(rawcounts)
#colnames(rawcounts)[2] <- "Primer_name"


#primers$Primer_name <- gsub(pattern = " ", replacement = "", x= primers$Primer_name)

#setdiff(richness$Primer_name, nsamp$Primer_name)

#richness <- merge(richness, nsamp, by= "Primer_name")

#richness <- richness %>% distinct(Primer_name, .keep_all = TRUE)

#richness$Primer_name <- gsub(pattern = "-", replacement = "_", x = richness$Primer_name)

#richness <- merge(richness, rawcounts, "Primer_name", all= T) %>% distinct(Primer_name, .keep_all = TRUE)

#richness <- merge(richness, primerInput, "Primer_name")

#richness <- richness[which(richness$num.reads>8),] 

#Div <- richness$Shannon
#Rich <- log(richness$Observed + 1)
#evenness = Div/Rich
#richness$Evenness = evenness

###Plots

#ob <- ggplot(data = richness, aes(x = Gen, y = na.omit(Observed), color= Target)) +
#  geom_point() +
  #facet_grid(.~Gen) +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Observed richness")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#ggplot(data = richness, aes(x = Gen, y = na.omit(Observed), fill= Target)) +
#  geom_boxplot() +
#  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  #facet_grid(.~Gen) +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Observed richness")+
  #scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#ch1<- ggplot(data = richness, aes(x = Gen, y = na.omit(Chao1), color= Target)) +
#  geom_point() +
  #facet_grid(.~Gen) +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Chao-1 richness")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#sh <- ggplot(data = richness, aes(x = Gen, y = na.omit(Shannon), color= Target)) +
#  geom_point() +
  #facet_grid(.~Gen) +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Diversity (Shannon Index)")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#sim<- ggplot(data = richness, aes(x = Gen, y = na.omit(Simpson), color= Target)) +
#  geom_point() +
  #facet_grid(.~Gen) +
 # theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Diversity (Simpson Index)")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#isim<-ggplot(data = richness, aes(x = Gen, y = na.omit(InvSimpson), color= Target)) +
#  geom_point() +
  #facet_grid(.~Gen) +
 # theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Diversity (Inverse Simpson Index)")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#eve<- ggplot(data = richness, aes(x = Gen, y = na.omit(Evenness), color= Target)) +
#  geom_point() +
  #facet_grid(.~Gen) +
 # theme_classic() +
 # theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
 # scale_y_continuous(name = "Evenness (Shannon index/ log(Observed richness +1))")+
 # scale_colour_brewer(palette = "Set1") +
 # theme(legend.position="right") 


#grid.arrange(ob, ch1, sh, sim, isim, eve, nrow= 2, ncol= 3)

#linob <- ggplot(data = richness, aes(x = log10(num.reads), y = Observed, color= Target)) +
 #       geom_point() +
        #facet_grid(.~Gen) +
 #       theme_classic() +
 #       geom_smooth(method = "lm", col = "black") +
        #theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#        scale_x_continuous(name = "log10 (Sequencing reads)")+
 #       scale_colour_brewer(palette = "Set1") +
 #       theme(legend.position="right") 

#linshan <- ggplot(data = richness, aes(x = log10(num.reads), y = Shannon, color= Target)) +
 #         geom_point() +
          #facet_grid(.~Gen) +
 #         theme_classic() +
 #         geom_smooth(method = "lm", col = "black") +
          #theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
 #         scale_x_continuous(name = "log10 (Sequencing reads)")+
 #         scale_colour_brewer(palette = "Set1") +
 #         theme(legend.position="right") 


#lineve <-  ggplot(data = richness, aes(x = log10(num.reads), y = Evenness, color= Target)) +
#          geom_point() +
          #facet_grid(.~Gen) +
#          theme_classic() +
#          geom_smooth(method = "lm", col = "black") +
          #theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
 #         scale_x_continuous(name = "log10 (Sequencing reads)")+
#          scale_colour_brewer(palette = "Set1") +
#          theme(legend.position="right") 

#summary(lm(log10(richness$Raw_counts)~ richness$Observed))$adj.r.squared

###Plot raw counts by primer pair 
#ggplot(richness, aes(x= reorder(Primer_name, -num.reads), y=num.reads, color= Gen, fill= Gen)) + 
#  geom_bar(stat = "identity")+
#  scale_y_log10(name = "log10 Total sequencing reads")+
#  theme_minimal() +
#  labs(x= "Primer name")+
#  geom_hline(yintercept=100, linetype="dashed", color = "black") +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1))


###Amplicon size vs raw counts
require(ggpubr)
ggplot(AbPhy, aes(x = Expected, y = Total_reads), geom=c("point", "smooth")) +
  scale_x_continuous(name = "Expected amplicon size") +
  #scale_y_continuous(name = "Sequencing reads") + #, limits=c(0.95, 1.60)) + ggtitle("Oocysts L/W ratio by group")+
  scale_y_log10(name = "log10 Total sequencing reads")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), aes(size= 25))+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  #geom_hline(yintercept=10000, linetype="dashed", color = "black") + 
  theme(text = element_text(size=20))+
  labs(tag = "A)")#+
  #stat_cor(label.x = 500, label.y = log10(3000000), aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  #stat_regline_equation(label.x = 500, label.y = log10(2200000))

ggplot(AbPhy, aes(x = Expected, y = Total_reads, color= Gen), geom=c("point", "smooth")) +
     scale_x_continuous(name = "Expected amplicon size") +
    #scale_y_continuous(name = "Sequencing reads") + #, limits=c(0.95, 1.60)) + ggtitle("Oocysts L/W ratio by group")+
     scale_y_log10(name = "log10 Total sequencing reads")+ 
     geom_jitter(shape=16, position=position_jitter(0.2), aes(size= 25))+
     theme_bw() +
     theme(legend.text=element_text(size=20)) +
     theme(legend.key.size = unit(3,"line")) +
     geom_smooth(method = "lm", se = FALSE, col = "red") +
     guides(colour = guide_legend(override.aes = list(size=10))) +
     #geom_hline(yintercept=10000, linetype="dashed", color = "black") + 
     theme(text = element_text(size=20))+
     labs(tag = "B)")

#######cool function to analyse diversity by primer pair########

#names(PS.l)<- names(STNC[keep])


#setdiff(names(PS.l), as.vector(primerInput$Primer_name))
#setdiff(names(PS.l.f), as.vector(primerInput$Primer_name))
#setdiff(names(PS.l.r), as.vector(primerInput$Primer_name))

#Index<- "Shannon" ###Change according to the index is required 
#vector_of_primers<- as.vector(primerInput$Primer_name)

#Pdiversity <- function(vector_of_primers, PS.l, Index)
#{
#  finalPrime <- list(list())
#  for(i in vector_of_primers)
#  {
#    if(i %in% names(PS.l))
#    {
#      id <- which(i == names(PS.l))
#      finalPrime[[i]][["richness_index"]] <- data.frame(estimate_richness(PS.l[[id]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson")))
#      Div <- finalPrime[[i]][["richness_index"]]$Shannon
#      Rich <- log(finalPrime[[i]][["richness_index"]]$Observed + 1)
#      evenness = Div/Rich
#      finalPrime[[i]][["richness_index"]]$Evenness = evenness
#      finalPrime[[i]][["number_of_samples"]] <- nrow(otu_table(PS.l[[id]])) ##Not correct 
#      mean_ind <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = mean)
#      std_dev <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = sd)
#      ci <- qnorm(0.95)*std_dev/sqrt(finalPrime[[i]][["number_of_samples"]])
#      upper <- mean_ind+ci
#      lower <- mean_ind-ci
#      finalPrime[[i]][["central_tendency"]] <- data.frame(Index_name = c("Observed","Chao1", "se.chao1", "Shannon", "Simpson", "InvSimpson", "Evenness"), 
#                                                          Mean = mean_ind, 
#                                                          SD = std_dev,
#                                                          Upper = upper, 
#                                                          Lower = lower)
 #     finalPrime[[i]][["primer_info"]] <- data.frame(primerInput[ primerInput$Primer_name %in% i, c(6:ncol(primerInput))])
#    }
#  }
#  finalPrime[[1]] <- NULL
  # plots
#  plot_df <- data.frame()
#  x <- data.frame()
#  for(j in names(finalPrime))
#  {
#    x <- finalPrime[[j]][["richness_index"]]
#    x[1:nrow(x),8] <- j
#    x[1:nrow(x),9] <- finalPrime[[j]][["primer_info"]]$Gen
#    x[1:nrow(x),10] <- finalPrime[[j]][["primer_info"]]$Target
#    plot_df <- rbind(plot_df, x)
#  }
#  colnames(plot_df)[8] <- c("primer_name")
#  colnames(plot_df)[9] <- "gene"
#  colnames(plot_df)[10] <- "target"
  
#  ggplot(plot_df, aes(x= primer_name, y = plot_df[,Index], fill = gene))+
#    geom_boxplot(outlier.colour="black", alpha = 0.5)+
#    geom_jitter(shape=1)+
#    xlab("Primer name")+
#    ylab(paste0(Index, " Diversity Index", collapse = ''))+ ###Change according to the index that is ploted 
#    theme_classic()+
#    theme(axis.text.x = element_text(angle = 90, vjust = 1))
  
  
  # upsetr
#  require(UpSetR)
  # We need the list of sample names for a given vector of primers
#  upset_list <- list()
#  for(k in vector_of_primers)
#  {
#    upset_list[[k]] <- rownames(finalPrime[[k]][["richness_index"]])
#  }
#  upset(fromList(upset_list), order.by = "freq", mainbar.y.label = "Number of samples", sets.x.label = "Sample size", nintersects = NA, nsets = 41)
  
#}


######Group specific analysis####
####By marker 
###First merge PrimTax with primerInput
PrimTax$Primer_name <- gsub(pattern = " ", replacement = "", x = PrimTax$Primer_name) ### check if the primer names have extra spaces
intersect(PrimTax$Primer_name, primerInput$Primer_name) 
#setdiff(PrimTax$Primer_name, primerInput$Primer_name)

library(plyr)
PrimTax <-join(PrimTax, primerInput) ###merge the selected information with the origial data frame created 
PrimTax$X<-  NULL

####Just 18S
if(SSU_18S){
Primtax.comb18 <- subset(PrimTax, Gen== "18S")

ps.numord.l.comb18 <- PS.l[c(Primtax.comb18$Primer_name)][order(Primtax.comb18$num.reads, decreasing=TRUE)]


cum.tax.comb18 <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.comb18, rank)
})

colnames(cum.tax.comb18) <- paste0("cum.", colnames(cum.tax.comb18))


cum.tax.comb18 <- as.data.frame(cbind(1:nrow(cum.tax.comb18), cum.tax.comb18))
#prnames28s <- names(ps.numord.l.28S) ## get names from the primers
#cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
colnames(cum.tax.comb18)[1] <- "Primer_name"

cum.plot.comb18 <- ggplot(cum.tax.comb18) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(20, 2200, label="Species"), color="red")+ 
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 1600, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 840, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 380, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 135, label="Classes"), color="#F0F921FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(20, 40, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa") +
  #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank())+
  labs(tag = "A)")+
  coord_cartesian(ylim = c(20, 3800))
}
####18S and 28S
if(SSU_LSU){
Primtax.comb28 <- subset(PrimTax, Gen=="28S" | Gen== "18S")

ps.numord.l.comb28 <- PS.l[c(Primtax.comb28$Primer_name)][order(Primtax.comb28$num.reads, decreasing=TRUE)]


cum.tax.comb28 <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.comb28, rank)
})

colnames(cum.tax.comb28) <- paste0("cum.", colnames(cum.tax.comb28))


cum.tax.comb28 <- as.data.frame(cbind(1:nrow(cum.tax.comb28), cum.tax.comb28))
#prnames28s <- names(ps.numord.l.28S) ## get names from the primers
#cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
colnames(cum.tax.comb28)[1] <- "Primer_name"

cum.plot.comb28 <- ggplot(cum.tax.comb28) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(22, 3300, label="Species"), color="red")+ 
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(22, 2000, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(22, 1000, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(22, 420, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 135, label="Classes"), color="#F0F921FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(22, 42, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa") +
  #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank())+
  labs(tag = "B)")+
  coord_cartesian(ylim = c(20, 3800))

}
####28S, 18S, COI####
if(SSU_LSU_COI){
Primtax.comb2 <- subset(PrimTax, Gen=="28S" | Gen== "18S" | Gen== "COI")

ps.numord.l.comb2 <- PS.l[c(Primtax.comb2$Primer_name)][order(Primtax.comb2$num.reads, decreasing=TRUE)]


cum.tax.comb2 <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.comb2, rank)
})

colnames(cum.tax.comb2) <- paste0("cum.", colnames(cum.tax.comb2))


cum.tax.comb2 <- as.data.frame(cbind(1:nrow(cum.tax.comb2), cum.tax.comb2))
#prnames28s <- names(ps.numord.l.28S) ## get names from the primers
#cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
colnames(cum.tax.comb2)[1] <- "Primer_name"

cum.plot.comb2 <- ggplot(cum.tax.comb2) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(26, 3800, label="Species"), color="red")+ 
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(26, 2200, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(26, 1100, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(26, 450, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 135, label="Classes"), color="#F0F921FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(26, 43, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa") +
  #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank())+
  labs(tag = "C)")+
  coord_cartesian(ylim = c(20, 3800))

grid.arrange(cum.plot.comb18, cum.plot.comb28, cum.plot.comb2, nrow= 1, ncol= 3,
             bottom= textGrob("Number of primers considerd (starting with the one with highest read count)"))
}
#####Parasites All primers (Dodgy solution, but working!!!!)
Primtax.comb <- subset(PrimTax, Gen=="28S" | Gen== "18S" | Gen== "COI" | Gen=="16S")
ps.numord.l.comb <- PS.l[c(Primtax.comb$Primer_name)][order(Primtax.comb$num.reads, decreasing=TRUE)]

ps.numord.l.comb[["L2513_64_F.H2714_64_R"]] <- NULL
ps.numord.l.comb[["Ins16S_1shortF_70_F.Ins16S_1shortR_70_R"]] <- NULL
ps.numord.l.comb[["Coleop_16Sc_68_F.Coleop_16Sd_68_R"]] <- NULL
ps.numord.l.comb[["F1073_101_F.1492R_100_R"]] <- NULL
ps.numord.l.comb[["ADM330_18_F.Klin0785_CR_19_R"]] <- NULL
ps.numord.l.comb[["16S.1100.F16_100_F.1492R_100_R"]] <- NULL
ps.numord.l.comb[["ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R"]]<- NULL

PM.para.comb <- lapply(ps.numord.l.comb, function (x) {
  subset_taxa(x, phylum%in%c("Nematoda",
                             "Apicomplexa",
                             "Platyhelminthes"))
})


PG.para.comb <- lapply(PM.para.comb, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


mat.para.comb <- lapply(PG.para.comb, function (x) {
  as.matrix(tax_table(x))
}) 

P.comb <- lapply(PG.para.comb, function (y) {
  make.names(tax_table(y)[, 5])
}) 

P.comb$L2513_64_F.H2714_64_R<- list()
P.comb$Ins16S_1shortF_70_F.Ins16S_1shortR_70_R<- list()
P.comb$Coleop_16Sc_68_F.Coleop_16Sd_68_R<- list()
P.comb$Ave12F_91_F.Ave12R_92_R<- list()
P.comb$F1073_101_F.1492R_100_R<- list()
P.comb$ADM330_18_F.Klin0785_CR_19_R<- list()
P.comb$`16S.1100.F16_100_F.1492R_100_R`<- list()
P.comb$S2F_80_F.S3R_81_R<- list()
P.comb$F52_69_F.R193_69_R<- list()
P.comb$ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R<- list()

####compair all the primers for multiple intersections

P.intersection <- sapply(seq_len(length(P.comb)), function(x)
  sapply(seq_len(length(P.comb)), function(y) length(intersect(unlist(P.comb[x]), unlist(P.comb[y]))))
)

rownames(P.intersection)<- names(P.comb)
colnames(P.intersection)<- names(P.comb)

P.clust <- hclust(dist(P.intersection), method = "complete") ##Dendogram
#library(dendextend)
as.dendrogram(P.clust) %>%
  plot(horiz = TRUE)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))
P.col$Primer_name <- rownames(P.col)

P.col <- merge(P.col, PrimTax, by="Primer_name", sort= F)

col_groups <- P.col %>%
  select("Primer_name", "Gen", "Region") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Primer_name

col_groups$Primer_name<- NULL

colour_groups <- list( Gen= c("18S"= "#440154FF", 
                              "28S"= "#21908CFF", 
                              "COI"= "#FDE725FF", 
                              "16S"= "pink",
                              "12S"= "#E3DAC9",
                              "ITS"= "#C46210",
                              "rbcL" = "#D0FF14"))

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

ggsave("~/AA_Primer_evaluation/Figures/Parasites_redundancy_allprimers.pdf", plot = paraheatmap, dpi = 450, width = 14, height = 10)

####Cummulative curve for parasites (All primers)
num.taxa.para <-sapply(c( "species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.para.comb, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.para <- unlist(lapply(PM.para.comb, function (x) 
  sum(otu_table(x))))

PrimTax.para <-as.data.frame(cbind(num.taxa.para,as.data.frame(num.reads.para)))
PrimTax.para$Primer_name<- rownames(PrimTax.para)
markers <- PrimTax[,c(8,16)]
PrimTax.para <- join(PrimTax.para, markers)
PrimTax.para <- PrimTax.para[order(-PrimTax.para$num.reads.para),] 

ps.numord.l.para <- PM.para.comb[c(PrimTax.para$Primer_name)][order(as.vector(unlist(PrimTax.para$num.reads.para)), decreasing=TRUE)]

cum.tax.para <- sapply(c("species","genus", "family","order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.para, rank)
})

colnames(cum.tax.para) <- paste0("cum.", colnames(cum.tax.para))

cum.tax.para <- as.data.frame(cbind(1:nrow(cum.tax.para), cum.tax.para))
cum.tax.para[2:30,7] <- PrimTax.para$Gen ## add marker amplified 
colnames(cum.tax.para)[1] <- "Primer_name"
colnames(cum.tax.para)[7] <- "Gen"
cum.tax.para<- add_row(cum.tax.para, Primer_name= 31:40)
cum.tax.para[31:40,2:6]<-cum.tax.para[30,2:6] 
cum.tax.para[31:40,7]<-col_groups[30:39,1] 

cum.plot.para <- ggplot(cum.tax.para) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(26, 380, label="Species"), color="red") +
  geom_point(aes(Primer_name, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("#E3DAC9","pink","#440154FF","#21908CFF","#FDE725FF","#C46210","#D0FF14"))+
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(26, 180, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(26, 110, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(26, 30, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 10, label="Classes"), color="#F0F921FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(26, 4, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa (Parasites)") +
  #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = 'none', axis.title.x = element_blank())+
  labs(tag = "B)")#+
#coord_cartesian(ylim = c(20, 3800))
#ggsave("~/AA_Primer_evaluation/Figures/Parasites_cummulative.pdf", plot = cum.plot.para, dpi = 450, width = 14, height = 10)

#####Fungi with all primers 

ps.numord.l.comb[["Klin0341_19_F.Klin0785_CR_19_R"]]<- NULL
ps.numord.l.comb[["27M_F_98_F.CCR_98_R"]]<- NULL

PM.fungi.comb <- lapply(ps.numord.l.comb, function (x) {
  subset_taxa(x, phylum%in%c("Ascomycota",
                             "Basidiomycota",
                             "Zygomycota", "Mucoromycota", "Zoopagomycota", "Blastocladiomycota", "Cryptomycota", "Microsporidia"))
})

PG.fungi.comb <- lapply(PM.fungi.comb, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})

mat.fungi.comb <- lapply(PG.fungi.comb, function (x) {
  as.matrix(tax_table(x))
}) 


F.comb <- lapply(PG.fungi.comb, function (y) {
  make.names(tax_table(y)[,5])
})

###Add empty primers
F.comb$L2513_64_F.H2714_64_R<- list()
F.comb$Ins16S_1shortF_70_F.Ins16S_1shortR_70_R<- list()
F.comb$Coleop_16Sc_68_F.Coleop_16Sd_68_R<- list()
F.comb$Ave12F_91_F.Ave12R_92_R<- list()
F.comb$F1073_101_F.1492R_100_R<- list()
F.comb$ADM330_18_F.Klin0785_CR_19_R<- list()
F.comb$`16S.1100.F16_100_F.1492R_100_R`<- list()
F.comb$S2F_80_F.S3R_81_R<- list()
F.comb$F52_69_F.R193_69_R<- list()
F.comb$ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R<- list()
F.comb$Klin0341_19_F.Klin0785_CR_19_R<- list()
F.comb$`27M_F_98_F.CCR_98_R`<- list()

####compair all the primers for multiple intersections

F.intersection <- sapply(seq_len(length(F.comb)), function(x)
  sapply(seq_len(length(F.comb)), function(y) length(intersect(unlist(F.comb[x]), unlist(F.comb[y]))))
)

rownames(F.intersection)<- names(F.comb)
colnames(F.intersection)<- names(F.comb)

F.clust <- hclust(dist(F.intersection), method = "complete") ##Dendogram
#library(dendextend)
as.dendrogram(F.clust) %>%
  plot(horiz = TRUE)

F.col <- cutree(tree = F.clust, k = 2)
F.col  <- data.frame(cluster = ifelse(test = F.col  == 1, yes = "cluster 1", no = "cluster 2"))
F.col$Primer_name <- rownames(F.col)

F.col <- merge(F.col, PrimTax, by="Primer_name", sort= F)

col_groups <- F.col %>%
  select("Primer_name", "Gen", "Region") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Primer_name

col_groups$Primer_name<- NULL

colour_groups <- list( Gen= c("18S"= "#440154FF", 
                              "28S"= "#21908CFF", 
                              "COI"= "#FDE725FF", 
                              "16S"= "pink",
                              "12S"= "#E3DAC9",
                              "ITS"= "#C46210",
                              "rbcL" = "#D0FF14"))

fungheatmap <- pheatmap(F.intersection, 
                        color = plasma(100),
                        border_color = NA,
                        annotation_col = col_groups, 
                        #annotation_row = col_groups,
                        annotation_colors = colour_groups,
                        #cutree_rows = 2,
                        #cutree_cols = 2,
                        show_rownames = F,
                        show_colnames = F,
                        main= "Redundant genus of fungi amplified")

ggsave("~/AA_Primer_evaluation/Figures/Fungi_redundancy_allprimers.pdf", plot = fungheatmap, dpi = 450, width = 14, height = 10)

####Cummulative curve for fungi (All primers)
num.taxa.fungi <-sapply(c( "species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.fungi.comb, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.fungi <- unlist(lapply(PM.fungi.comb, function (x) 
  sum(otu_table(x))))

PrimTax.fungi <-as.data.frame(cbind(num.taxa.fungi,as.data.frame(num.reads.fungi)))
PrimTax.fungi$Primer_name<- rownames(PrimTax.fungi)
markers <- PrimTax[,c(8,16)]
PrimTax.fungi <- join(PrimTax.fungi, markers)
PrimTax.fungi <- PrimTax.fungi[order(-PrimTax.fungi$num.reads.fungi),] 

ps.numord.l.fungi <- PM.fungi.comb[c(PrimTax.fungi$Primer_name)][order(as.vector(unlist(PrimTax.fungi$num.reads.fungi)), decreasing=TRUE)]

cum.tax.fungi <- sapply(c("species","genus", "family","order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.fungi, rank)
})

colnames(cum.tax.fungi) <- paste0("cum.", colnames(cum.tax.fungi))

cum.tax.fungi <- as.data.frame(cbind(1:nrow(cum.tax.fungi), cum.tax.fungi))
cum.tax.fungi[2:28,7] <- PrimTax.fungi$Gen ## add marker amplified 
colnames(cum.tax.fungi)[1] <- "Primer_name"
colnames(cum.tax.fungi)[7] <- "Gen"
cum.tax.fungi<- add_row(cum.tax.fungi, Primer_name= 29:40)
cum.tax.fungi[29:40,2:6]<-cum.tax.fungi[28,2:6] 
cum.tax.fungi[29:40,7]<-col_groups[28:39,1] 

cum.plot.fungi <- ggplot(cum.tax.fungi) +
    geom_step(aes(Primer_name, cum.species), color="red") +
    geom_text(aes(26, 2000, label="Species"), color="red") +
    geom_point(aes(Primer_name, cum.species, colour = Gen, size= 2))+
    scale_color_manual(values = c("#E3DAC9","pink","#440154FF","#21908CFF","#FDE725FF","#C46210","#D0FF14"))+
    geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
    geom_text(aes(26, 1000, label="Genera"), color="#0D0887FF") +
    geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
    geom_text(aes(26, 440, label="Families"), color="#7E03A8FF") +
    geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
    geom_text(aes(26, 160, label="Orders"), color="#CC4678FF") +
    #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
    #geom_text(aes(20, 50, label="Classes"), color="") +
    geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
    geom_text(aes(26, 8, label="Phyla"), color="#F89441FF") +
    scale_y_log10("Cummulative count of taxa (Fungi)") +
    #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
    annotation_logticks(sides="l") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = 'none', axis.title.x = element_blank())+
    labs(tag = "C)")#+
  
  CummParaFung <- grid.arrange(cum.plot.all, cum.plot.para, cum.plot.fungi, nrow= 1, ncol= 3,  
                               bottom= textGrob("Number of primers considerd (starting with the one with highest read count)"))
  ###Put parasites and fungi together 
  ggsave("~/AA_Primer_evaluation/Figures/ParaFung_cummulative_allprimers.pdf", plot = CummParaFung, dpi = 450, width = 14, height = 10)
  
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

lapply(PG.para.comb2, function (x) {
  table(tax_table(x)[,5])
}) 


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
  theme(legend.position = "none")+
  labs(tag = "A)")


Parasites <- Parasites.comb2$Genus
Parasites <- as.data.frame(table(Parasites))
Parasites <- Parasites[order(-Parasites$Freq),]

###Add a bar plot (To improve)

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
  select("Primer_name", "Gen", "Region", "num.reads") ##Here It is possible to add the expected size 

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

ggsave("~/AA_Primer_evaluation/Figures/Parasites_redundancy.pdf", plot = paraheatmap, dpi = 450, width = 14, height = 10)


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

cum.plot.para <- ggplot(cum.tax.para) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(26, 380, label="Species"), color="red") +
  geom_point(aes(Primer_name, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))+
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(26, 180, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(26, 110, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(26, 30, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 10, label="Classes"), color="#F0F921FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(26, 4, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa (Parasites)") +
  #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = 'none', axis.title.x = element_blank())+
  labs(tag = "A)")#+
  #coord_cartesian(ylim = c(20, 3800))
#ggsave("~/AA_Primer_evaluation/Figures/Parasites_cummulative.pdf", plot = cum.plot.para, dpi = 450, width = 14, height = 10)


###Fungi 18s + 28S + COI

PM.fungi.comb2 <- lapply(ps.numord.l.comb2, function (x) {
  subset_taxa(x, phylum%in%c("Ascomycota",
                             "Basidiomycota",
                             "Zygomycota", "Mucoromycota", "Zoopagomycota", "Blastocladiomycota", "Cryptomycota", "Microsporidia"))
})

PG.fungi.comb2 <- lapply(PM.fungi.comb2, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


lapply(PG.fungi.comb2, function (x) {
  table(tax_table(x)[,5])
}) 


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
        axis.text.y=element_blank())+
  labs(tag = "B)")

grid.arrange(paraprimer, fungiprimer, nrow= 1, ncol= 2)

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

PFPlot <- grid.arrange(paraprimer, fungiprimer, nrow= 1, ncol= 2) ##Parasites vs Fungi
ggsave("~/AA_Primer_evaluation/Figures/Parasites_Fungi_Primers.pdf", plot = PFPlot, dpi = 450, width = 14, height = 10)


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
  select("Primer_name", "Gen", "Region", "num.reads") ##, "Expected"

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

ggsave("~/AA_Primer_evaluation/Figures/Fungi_redundancy.pdf", plot = fungiheatmap, dpi = 450, width = 14, height = 10)



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

cum.plot.fungi <- ggplot(cum.tax.fungi) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(26, 2000, label="Species"), color="red") +
  geom_point(aes(Primer_name, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))+
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(26, 1000, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(26, 440, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(26, 160, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 50, label="Classes"), color="") +
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(26, 8, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa (Fungi)") +
  #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = 'none', axis.title.x = element_blank())+
  labs(tag = "B)")#+

CummParaFung <- grid.arrange(cum.plot.para, cum.plot.fungi, nrow= 1, ncol= 2,  
             bottom= textGrob("Number of primers considerd (starting with the one with highest read count)"))
 ###Put parasites and fungi together 
ggsave("~/AA_Primer_evaluation/Figures/ParaFung_cummulative.pdf", plot = CummParaFung, dpi = 450, width = 14, height = 10)


####Turtles project 

###Good primers for turtles
ps.chlorord.l <- ps.numord.l

ps.chlorord.l[["Ave12F_91_F.Ave12R_92_R"]]<- NULL
ps.chlorord.l[["Coleop_16Sc_68_F.Coleop_16Sd_68_R"]]<- NULL ### No chlorophyta
ps.chlorord.l[["fwhF2_106_F.fwhR2n_106_R"]]<- NULL
ps.chlorord.l[["S2F_80_F.S3R_81_R"]]<- NULL
ps.chlorord.l[["F52_69_F.R193_69_R"]]<- NULL
ps.chlorord.l[["Prot1702_32_F.EukB_MH_24_R"]]<- NULL
ps.chlorord.l[["27M_F_98_F.CCR_98_R"]]<- NULL
ps.chlorord.l[["ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R"]]<- NULL
ps.chlorord.l[["515F_99_F.R1073_103_R"]]<- NULL
ps.chlorord.l[["F1073_101_F.1492R_100_R"]]<- NULL
ps.chlorord.l[["ADM330_18_F.Klin0785_CR_19_R"]]<- NULL
ps.chlorord.l[["L2513_64_F.H2714_64_R"]]<- NULL
ps.chlorord.l[["16SA_L_110_F.16SB_H_110_R"]]<- NULL

ps.chlorord.l <- lapply(ps.chlorord.l, function (x) {
  subset_taxa(x, phylum%in%"Chlorophyta")
})

ps.chlorord.l <- lapply(ps.chlorord.l, function (x) {
  tax_glom(x, "species", NArm = TRUE)
})

sumSeqByTax(PS.l = ps.chlorord.l, rank = "species")

num.taxa.chlorophyta <-sapply(c("species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(ps.chlorord.l, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.chlorophyta <- unlist(lapply(ps.chlorord.l, function (x) 
  sum(otu_table(x))))


PrimTax.chloro <-as.data.frame(cbind(num.taxa.chlorophyta,as.data.frame(num.reads.chlorophyta)))
PrimTax.chloro$Primer_name<- rownames(PrimTax.chloro)
markers <- PrimTax[,c(8,16)]
PrimTax.chloro <- join(PrimTax.chloro, markers)
PrimTax.chloro <- PrimTax.chloro[order(-PrimTax.chloro$num.reads.chlorophyta),] 

ps.numord.l.chloro <- ps.chlorord.l[c(PrimTax.chloro$Primer_name)][order(as.vector(unlist(PrimTax.chloro$num.reads.chlorophyta)), decreasing=TRUE)]


cum.tax.chloro <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.chloro, rank)
})

colnames(cum.tax.chloro) <- paste0("cum.", colnames(cum.tax.chloro))


cum.tax.chloro <- as.data.frame(cbind(1:nrow(cum.tax.chloro), cum.tax.chloro))
cum.tax.chloro[2:29,7] <- PrimTax.chloro$Gen ## add marker amplified 
colnames(cum.tax.chloro)[1] <- "Primer_name"
colnames(cum.tax.chloro)[7] <- "Gen"

cum.plot.chloro <- ggplot(cum.tax.chloro) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(20, 250, label="Species"), color="red") +
  geom_point(aes(Primer_name, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("pink","#440154FF", "#21908CFF", "#FDE725FF"))+
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 120, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 50, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 20, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  #geom_text(aes(20, 9, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
  geom_text(aes(20, 2, label="Phyla"), color="#F0F921FF") +
  scale_y_log10("Cummulative count of taxa (Chlorophyta)") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

####compair all the primers for multiple intersections

ps.numord.l.chloro2 <- lapply(ps.numord.l.chloro, function (x) {
  tax_glom(x, "species", NArm = TRUE)
})

###These list is used for heatmap later!!!!
chloro.ord <- lapply(ps.numord.l.chloro, function (y) {
  make.names(tax_table(y)[,5]) #5 Genus
}) 

Chlorophyta <- data.frame() ###Create the data frame 

for (i in 1: length(chloro.ord)) ### Start a loop: fro every element in the list ...
{ 
  orderch <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(chloro.ord[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    orderch[1,1] <- 0    ### Add a zero in the first column
    orderch[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    orderch <- as.data.frame((chloro.ord[[i]]))  ###Make a data frame with the data included in each element of the list 
    orderch[,2] <- rownames(orderch) ### And use the rownames as information of the second column 
  }
  
  orderch[,3] <- names(chloro.ord)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(orderch) <- c("Order", "Count", "Primer_name") ### change the names for the columns 
  Chlorophyta <- rbind(Chlorophyta, orderch) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

Chlorophyta <- unique(Chlorophyta[c("Order", "Primer_name")]) ##Take unique combinations 
Chlorophyta <- join(Chlorophyta, PrimTax.chloro) 

Chloroprimer <-ggplot(Chlorophyta, aes(Primer_name, color=Gen, fill=Gen))+
  geom_bar()+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(name = "Number of Chlorophyta orders amplified")+
  scale_color_manual(values = c("pink","#440154FF", "#21908CFF", "#FDE725FF"))+
  scale_fill_manual(values = c("pink","#440154FF", "#21908CFF", "#FDE725FF")) +
  scale_x_discrete(name = "Primer combination")+
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank())

####compair all the primers for multiple intersections

C.intersection <- sapply(seq_len(length(chloro.ord)), function(x)
  sapply(seq_len(length(chloro.ord)), function(y) length(intersect(unlist(chloro.ord[x]), unlist(chloro.ord[y]))))
)
rownames(C.intersection)<- names(chloro.ord)
colnames(C.intersection)<- names(chloro.ord)


C.clust <- hclust(dist(C.intersection), method = "complete") ##Dendogram
#library(dendextend)
as.dendrogram(C.clust) %>%
  plot(horiz = TRUE)

C.col <- cutree(tree = C.clust, k = 2)
C.col  <- data.frame(cluster = ifelse(test = C.col  == 1, yes = "cluster 1", no = "cluster 2"))

col_groups <- PrimTax.chloro %>%
  select("Gen")

row.names(col_groups)<- colnames(C.intersection)

colour_groups <- list( Gen= c("18S"= "#440154FF", "28S"= "#21908CFF", "COI"= "#FDE725FF", "16S"= "pink"))

c<- pheatmap(C.intersection, 
         color = plasma(100),
         border_color = "black",
         annotation_col = col_groups,
         #annotation_row = col_groups,
         annotation_colors = colour_groups,
         #cutree_rows = 2,
         #cutree_cols = 2,
         show_rownames = T,
         show_colnames = F,
         main= "Redundant genus of Chlorophyta amplified")

grid.arrange(cum.plot.chloro, c, nrow= 2, ncol= 1) 



