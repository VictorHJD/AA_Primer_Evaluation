####Analysis of primers from Fox and Raccoon access array experiment

###Loading libraries 
library(MultiAmplicon)
library(ggplot2)
library(data.table)
library(phyloseq)
library(dplyr)
library(vegan)
library(gridExtra)
library(grid)
library(lattice)
library("TmCalculator")
library("ranacapa")

##########Functions############
tabulate.taxa <- function(taxtab, taxon, phylumsubset){
  if(nrow(taxtab)>0){
    t <- taxtab[taxtab[, "phylum"]%in%phylumsubset, ]
    if(!is.null(ncol(t))){
      table(t[, taxon])
    } else {NULL} 
  }else {NULL} 
}

sumSeqByTax <- function (Phy, tax) {
  counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
  counts$asvCount <- as.numeric(as.character(counts$asvCount))
  tapply(counts$asvCount, counts[, tax], sum)
}

###########Sequencing and annotation data#########
###Load previously generated data from specific experiment

GetFox <- TRUE

GetRaccoon <- FALSE

GetTest <- FALSE

####Foxes

if(GetFox){
  ##Multiamplicon data 
  MA <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF_all.RDS")
  ###No chimeras data 
  STNC <- getSequenceTableNoChime(MA)
  ###Annotation list 
  annot.list <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Fox/Fox_blast_tax.Rds") 
  ##Discard primers without annotation 
  keep <- unlist(lapply(annot.list, nrow))>0
  names(STNC)[keep]<-gsub(pattern = "-", replacement = "_", x= names(STNC[keep]))
  names(STNC)[keep]<-gsub(pattern = " ", replacement = "", x= names(STNC[keep]))
  ### Data from samples 
  sample.data <- read.csv("/SAN/Victors_playground/Metabarcoding/AA_Fox/Metabarcoding_Fox_Final.csv")
  sample.data$samples <- gsub("_S\\d+_L001", "\\1", sample.data$samples)
  sample.data$samples <- gsub("-", "_", sample.data$samples)
  
  rownames(sample.data) <- sample.data$samples 
  ##Phyloseq list
  PS.l <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Fox/PhyloSeqList.Rds")
  names(PS.l) <- names(STNC)[keep]
}

####Raccoon

if(GetRaccoon){
  ##Multiamplicon data 
  MA <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Raccoon/MAR_all.RDS") 
  ###No chimeras data 
  STNC <- getSequenceTableNoChime(MA)
  ###Annotation list 
  annot.list <- annot.list <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/Racc_blast_tax_all.Rds")
  ##Discard primers without annotation 
  keep <- unlist(lapply(annot.list, nrow))>0
  names(STNC)[keep]<-gsub(pattern = "-", replacement = "_", x= names(STNC[keep]))
  names(STNC)[keep]<-gsub(pattern = " ", replacement = "", x= names(STNC[keep]))
  ### Data from samples 
  #sample.data <- read.csv("/SAN/Victors_playground/Metabarcoding/AA_Raccoon/")
  #rownames(sample.data) <- sample.data$samples 
  ##Phyloseq list
  PS.l <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/PhyloSeqList_Raccall.Rds")
  names(PS.l) <- names(STNC)[keep]
}

#### Test 
if(GetTest){
  ##Multiamplicon data 
  MA <- readRDS("/SAN/Victors_playground/Metabarcoding/MA_final.RDS")
  ###No chimeras data 
  STNC <- getSequenceTableNoChime(MA)
  ###Annotation list 
  annot.list <- readRDS(file="/SAN/Victors_playground/Metabarcoding/Fox_blast_tax.Rds")
  ##Discard primers without annotation 
  keep <- unlist(lapply(annot.list, nrow))>0
  names(STNC)[keep]<-gsub(pattern = "-", replacement = "_", x= names(STNC[keep]))
  names(STNC)[keep]<-gsub(pattern = " ", replacement = "", x= names(STNC[keep]))
  ### Data from samples 
  sample.data <- read.csv("/SAN/Victors_playground/Metabarcoding/Metabarcoding_Chip1_Chip2_CS_20180824.csv")
  ##Phyloseq list
  PS.l <- readRDS("/SAN/Victors_playground/Metabarcoding/PhyloSeqList.Rds")
  names(PS.l) <- names(STNC)[keep]
}
###Number of sequencing reads per primer per taxa 

readNumByPhylum <- lapply(PS.l, sumSeqByTax, "phylum")
names(readNumByPhylum) <- names(STNC)[keep]


readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") 
names(readNumByGenus) <- names(STNC)[keep]


readNumbyFam <- lapply(PS.l, sumSeqByTax, "family")
names(readNumbyFam) <- names(STNC)[keep]

###Summarize raw counts
### Total amount of read per primer pair 
rawcounts <- rowSums(getRawCounts(MA))
#rawcounts.F[order(rawcounts.F)]
rawcounts <- data.frame(rawcounts)
rawcounts[,2] <- rownames(rawcounts)

colnames(rawcounts) <- c("Raw_counts", "Primer_name")
rownames(rawcounts) <- c(1:nrow(rawcounts))
rawcounts <- data.frame(Primer_name = rawcounts$Primer_name, Raw_counts = rawcounts$Raw_counts) ###change the order of the columns
rawcounts$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcounts$Primer_name)
rawcounts$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcounts$Primer_name)
#sum(rawcounts$Raw_counts)

#write.csv(rawcounts, file = "/SAN/Victors_playground/Metabarcoding/AA_Fox/Raw_counts_primers_Fox.csv")

### Total amount of reads per sample 
rawsamp <- colSums(getRawCounts(MA)) 
#rawsamp[order(rawsamp)]
rawsamp <- data.frame(rawsamp)
rawsamp[,2] <- rownames(rawsamp)
colnames(rawsamp) <- c("Raw_counts", "Sample_ID")
rownames(rawsamp) <- c(1:nrow(rawsamp))
rawsamp <- data.frame(Sample_ID = rawsamp$Sample_ID, Raw_counts = rawsamp$Raw_counts) ###change the order of the columns
rawsamp$Sample_ID <- gsub(pattern = " ", replacement = "", x = rawsamp$Sample_ID)

#rawsamp$Sample_ID<- gsub("S\\d+_", "\\1", rawsamp$Sample_ID)
#rawsamp$Sample_ID<- gsub("d+_", "\\1", rawsamp$Sample_ID)

#sum(rawsamp$Raw_counts)

#write.csv(rawsamp, file = "/SAN/Victors_playground/Metabarcoding/AA_Fox/Raw_counts_per_Fox.csv")

################Primer information ##################
#primerInput <- read.csv("~/AA_Primer_evaluation/primerInput.csv") ###All primers even some duplicated
#Remove duplicate rows based on one or more column values (require dplyr): 
#primerInput <- primerInput %>% distinct(Primer_name, .keep_all = TRUE)
#write.csv(primerInput, "~/AA_Primer_evaluation/primerInputUnique.csv")

primerInput <- read.csv("/localstorage/victor/AA_Primer_evaluation/primerInputUnique.csv") 
primerInput$Primer_name <- gsub(pattern = " ", replacement = "", x = primerInput$Primer_name)
primerInput$X <- NULL
primerInput$Primer_name <- gsub(pattern = " ", replacement = "", x = primerInput$Primer_name)
##Primer characteristics (Tm and GC)
primerInput$Seq_F <- as.vector(primerInput$Seq_F)
primerInput$Seq_R <- as.vector(primerInput$Seq_R)

primerInput <- primerInput[-c(20,21,82),] ###Eliminating primers with inosine ... To make the loop work :S 

output <- data.frame(matrix(nrow = nrow(primerInput), ncol = 5))
colnames(output) <- c("Primer_name","Tm_F", "Tm_R", "GC_F", "GC_R")
output[,1]<-primerInput$Primer_name

for (row in 1: nrow(primerInput))
{
  Fwd <- primerInput[row, "Seq_F"]
  Rev <- primerInput[row, "Seq_R"]
  
  TmF <- Tm_NN(Fwd, ambiguous = T, imm_table = "DNA_IMM1")
  TmR <- Tm_NN(Rev, ambiguous = T, imm_table = "DNA_IMM1")
  GCF <- GC(Fwd, ambiguous = T)
  GCR <- GC(Rev, ambiguous = T)
  
  output[row, 2] <- TmF
  output[row, 3] <- TmR
  output[row, 4] <- GCF
  output[row, 5] <- GCR
  
}

primerInput <- read.csv("~/AA_Primer_evaluation/primerInputUnique.csv")
primerInput$Primer_name <- gsub(pattern = " ", replacement = "", x = primerInput$Primer_name)
primerInput <- merge(primerInput, output, by= "Primer_name", all= T)
primerInput$X <- NULL
primerInput$Primer_name <- gsub(pattern = " ", replacement = "", x = primerInput$Primer_name)
#write.csv(primerInput, "~/AA_Primer_evaluation/primerOutputUnique.csv")

n16S <- subset(primerInput, primerInput$Gen== "16S")
n18S <- subset(primerInput, primerInput$Gen== "18S")
n28S <- subset(primerInput, primerInput$Gen== "28S")
COI <-  subset(primerInput, primerInput$Gen== "COI")
plants <- subset(primerInput, primerInput$Target== "Viridiplantae")

###### Alpha diversity by marker in Fox dataset#####

richness <- data.frame()

for (i in 1:length(PS.l)){
  
  a <- data.frame()
  
  b <- data.frame(estimate_richness(PS.l[[i]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson")))
  
  a<-colMeans(b)
  
  richness<- rbind(richness, a)
}

richness[,7] <- names(STNC[keep])

colnames(richness) <- c("Observed","Chao1","se_Chao1", "Shannon", "Simpson", "InvSimpson", "Primer_name")

nsamp <- data.frame()

nsamp[1:41,1] <-names(STNC[keep])

for (i in 1:length(PS.l))
  nsamp[i,2]<- nrow(sample_data(PS.l[[i]]))


if(GetRaccoon){
  for (i in 1:length(PS.l))
    nsamp[i,2]<- nrow(data.frame(estimate_richness(PS.l[[i]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson"))))
}

colnames(nsamp) <- c("Primer_name", "Number_samples")

richness$Primer_name <- gsub(pattern = " ", replacement = "", x = richness$Primer_name) ### check if the primer names have extra spaces
nsamp$Primer_name <- gsub(pattern = " ", replacement = "", x= nsamp$Primer_name)

#primers$Primer_name <- gsub(pattern = " ", replacement = "", x= primers$Primer_name)

#setdiff(richness$Primer_name, nsamp$Primer_name)

richness <- merge(richness, nsamp, by= "Primer_name")

richness <- richness %>% distinct(Primer_name, .keep_all = TRUE)

richness$Primer_name <- gsub(pattern = "-", replacement = "_", x = richness$Primer_name)

richness <- merge(richness, rawcounts, "Primer_name", all= T) %>% distinct(Primer_name, .keep_all = TRUE)

richness <- merge(richness, primerInput, "Primer_name")

richness <- richness[which(richness$Raw_counts>8),] 

Div <- richness$Shannon
Rich <- log(richness$Observed + 1)
evenness = Div/Rich
richness$Evenness = evenness


ob <- ggplot(data = richness, aes(x = Gen, y = na.omit(Observed), color= Target)) +
  geom_point() +
  #facet_grid(.~Gen) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_y_continuous(name = "Observed richness")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") 

ch1<- ggplot(data = richness, aes(x = Gen, y = na.omit(Chao1), color= Target)) +
  geom_point() +
  #facet_grid(.~Gen) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_y_continuous(name = "Chao-1 richness")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") 

sh <- ggplot(data = richness, aes(x = Gen, y = na.omit(Shannon), color= Target)) +
  geom_point() +
  #facet_grid(.~Gen) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_y_continuous(name = "Diversity (Shannon Index)")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") 

sim<- ggplot(data = richness, aes(x = Gen, y = na.omit(Simpson), color= Target)) +
  geom_point() +
  #facet_grid(.~Gen) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_y_continuous(name = "Diversity (Simpson Index)")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") 

isim<-ggplot(data = richness, aes(x = Gen, y = na.omit(InvSimpson), color= Target)) +
  geom_point() +
  #facet_grid(.~Gen) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_y_continuous(name = "Diversity (Inverse Simpson Index)")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") 

eve<- ggplot(data = richness, aes(x = Gen, y = na.omit(Evenness), color= Target)) +
  geom_point() +
  #facet_grid(.~Gen) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_y_continuous(name = "Evenness (Shannon index/ log(Observed richness +1))")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") 


grid.arrange(ob, ch1, sh, sim, isim, eve, nrow= 2, ncol= 3)

ggplot(data = richness, aes(x = log10(Raw_counts), y = Observed, color= Target)) +
  geom_point() +
  #facet_grid(.~Gen) +
  theme_classic() +
  geom_smooth(method = "lm", col = "black") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_x_continuous(name = "log10 (Sequencing reads)")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") 

ggplot(data = richness, aes(x = log10(Raw_counts), y = Shannon, color= Target)) +
  geom_point() +
  #facet_grid(.~Gen) +
  theme_classic() +
  geom_smooth(method = "lm", col = "black") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_x_continuous(name = "log10 (Sequencing reads)")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") 


ggplot(data = richness, aes(x = log10(Raw_counts), y = Evenness, color= Target)) +
  geom_point() +
  #facet_grid(.~Gen) +
  theme_classic() +
  geom_smooth(method = "lm", col = "black") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_x_continuous(name = "log10 (Sequencing reads)")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") 

#summary(lm(log10(richness$Raw_counts)~ richness$Observed))$adj.r.squared

###Plot raw counts by primer pair 
ggplot(richness, aes(x= reorder(Primer_name, -Raw_counts), y=Raw_counts, color= Gen, fill= Gen)) + 
  geom_bar(stat = "identity")+
  scale_y_log10(name = "log10 Total sequencing reads")+
  theme_minimal() +
  labs(x= "Primer name")+
  geom_hline(yintercept=100, linetype="dashed", color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1))


###Amplicon size vs raw counts
ggplot(richness, aes(x = Expected, y = Raw_counts, color= Target), geom=c("point", "smooth")) +
  scale_x_continuous(name = "Expected amplicon size") +
  #scale_y_continuous(name = "Sequencing reads") + #, limits=c(0.95, 1.60)) + ggtitle("Oocysts L/W ratio by group")+
  scale_y_log10(name = "log10 Total sequencing reads")+ 
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.5, aes(size= 25))+
  theme_bw() +
  theme(legend.text=element_text(size=20)) +
  theme(legend.key.size = unit(3,"line")) +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  geom_hline(yintercept=10000, linetype="dashed", color = "black") + 
  theme(text = element_text(size=20))

#adjrsq <- summary(lm(pinfo$Amplicon_expected ~ pinfo$Raw_counts))$adj.r.squared

#######cool function########

names(PS.l)<- names(STNC[keep])


setdiff(names(PS.l), as.vector(primerInput$Primer_name))


Index<- "Shannon" ###Change according to the index is required 
vector_of_primers<- as.vector(primerInput$Primer_name)

Pdiversity <- function(vector_of_primers, PS.l, Index)
{
  finalPrime <- list(list())
  for(i in vector_of_primers)
  {
    if(i %in% names(PS.l))
    {
      id <- which(i == names(PS.l))
      finalPrime[[i]][["richness_index"]] <- data.frame(estimate_richness(PS.l[[id]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson")))
      Div <- finalPrime[[i]][["richness_index"]]$Shannon
      Rich <- log(finalPrime[[i]][["richness_index"]]$Observed + 1)
      evenness = Div/Rich
      finalPrime[[i]][["richness_index"]]$Evenness = evenness
      finalPrime[[i]][["number_of_samples"]] <- nrow(sample_data(PS.l[[id]])) ##For raccoon change 
      mean_ind <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = mean)
      std_dev <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = sd)
      ci <- qnorm(0.95)*std_dev/sqrt(finalPrime[[i]][["number_of_samples"]])
      upper <- mean_ind+ci
      lower <- mean_ind-ci
      finalPrime[[i]][["central_tendency"]] <- data.frame(Index_name = c("Observed","Chao1", "se.chao1", "Shannon", "Simpson", "InvSimpson", "Evenness"), 
                                                          Mean = mean_ind, 
                                                          SD = std_dev,
                                                          Upper = upper, 
                                                          Lower = lower)
      finalPrime[[i]][["primer_info"]] <- data.frame(primerInput[ primerInput$Primer_name %in% i, c(6:ncol(primerInput))])
    }
  }
  finalPrime[[1]] <- NULL
  # plots
  plot_df <- data.frame()
  x <- data.frame()
  for(j in names(finalPrime))
  {
    x <- finalPrime[[j]][["richness_index"]]
    x[1:nrow(x),8] <- j
    x[1:nrow(x),9] <- finalPrime[[j]][["primer_info"]]$Gen
    x[1:nrow(x),10] <- finalPrime[[j]][["primer_info"]]$Target
    plot_df <- rbind(plot_df, x)
  }
  colnames(plot_df)[8] <- c("primer_name")
  colnames(plot_df)[9] <- "gene"
  colnames(plot_df)[10] <- "target"
  
  ggplot(plot_df, aes(x= primer_name, y = plot_df[,Index], fill = gene))+
    geom_boxplot(outlier.colour="black", alpha = 0.5)+
    geom_jitter(shape=1)+
    xlab("Primer name")+
    ylab(paste0(Index, " Diversity Index", collapse = ''))+ ###Change according to the index that is ploted 
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 1))
  
  
  # upsetr
  require(UpSetR)
  # We need the list of sample names for a given vector of primers
  upset_list <- list()
  for(k in vector_of_primers)
  {
    upset_list[[k]] <- rownames(finalPrime[[k]][["richness_index"]])
  }
  upset(fromList(upset_list), order.by = "freq", mainbar.y.label = "Number of samples", sets.x.label = "Sample size", nintersects = NA, nsets = 41)
  
}


###############Table of taxa ##################
####Phylum
AbPhyF <- data.frame() ###Create the data frame 

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
  colnames(phyla) <- c("Reads", "Phyla", "Primer_name") ### change the names for the columns 
  AbPhyF <- rbind(AbPhyF, phyla) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

rownames(AbPhyF) <- c(1:nrow(AbPhyF)) ### change the rownames to consecutive numbers 
AbPhyF <- data.frame(Primer_name = AbPhyF$Primer_name, Phyla = AbPhyF$Phyla, Reads = AbPhyF$Reads) ###change the order of the columns

AbPhyF %>%
  group_by(Primer_name) %>% 
  mutate(Total_reads = sum(Reads)) -> AbPhyF ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

Relative_abundance = AbPhyF$Reads/AbPhyF$Total_reads ### create a vector with the result of the operation 

AbPhyF[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbPhyF)[5] <- "Relative_abundance" ### Change the name of the column 

AbPhyF$Primer_name <- gsub(pattern = " ", replacement = "", x = AbPhyF$Primer_name) ### check if the primer names have extra spaces

AbPhyF$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbPhyF$Primer_name)


#intersect(AbPhyF$Primer_name, primers$Primer_name)
#setdiff(AbPhyF$Primer_name, primers$Primer_name)
#setdiff(primers$Primer_name, AbPhyF$Primer_name)

AbPhyF <- merge(AbPhyF, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

#write.csv(AbPhyF, file = "/localstorage/victor/AA_Primer_evaluation/Phyloseq_Raccoon/Phyla_list_Racc_all.csv")


###Table of Genus 

AbGenF <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByGenus)) ### Start a loop: fro every element in the list ...
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
  colnames(genus) <- c("Reads", "Genus", "Primer_name") ### change the names for the columns 
  AbGenF <- rbind(AbGenF, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

rownames(AbGenF) <- c(1:nrow(AbGenF)) ### change the rownames to consecutive numbers 
AbGenF <- data.frame(Primer_name = AbGenF$Primer_name, genus = AbGenF$genus, Reads = AbGenF$Reads) ###change the order of the columns

AbGenF %>%
  group_by(Primer_name) %>% 
  mutate(Total_reads = sum(Reads)) -> AbGenF ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

Relative_abundance = AbGenF$Reads/AbGenF$Total_reads ### create a vector with the result of the operation 

AbGenF[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbGenF)[5] <- "Relative_abundance" ### Change the name of the column 

#unique(AbPhy$Primer_name)

AbGenF$Primer_name <- gsub(pattern = " ", replacement = "", x = AbGenF$Primer_name) ### check if the primer names have extra spaces
AbGenF$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbGenF$Primer_name)


#intersect(AbGenF$Primer_name, primers$Primer_name)
#setdiff(AbGenF$Primer_name, primers$Primer_name)
#setdiff(primers$Primer_name, AbGenF$Primer_name)

AbGenF <- merge(AbGenF, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

#write.csv(AbGenF, file = "/localstorage/victor/AA_Primer_evaluation/Phyloseq_Raccoon/Genus_list_Racc_all.csv")

#### Table by family

AbFamF <- data.frame() ###Create the data frame 

for (i in 1: length(readNumbyFam)) ### Start a loop: fro every element in the list ...
{ 
  family <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumbyFam[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    family[1,1] <- 0    ### Add a zero in the first column
    family[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    family <- as.data.frame((readNumbyFam[[i]]))  ###Make a data frame with the data included in each element of the list 
    family[,2] <- rownames(family) ### And use the rownames as information of the second column 
  }
  
  family[,3] <- names(readNumbyFam)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(family) <- c("Reads", "Family", "Primer_name") ### change the names for the columns 
  AbFamF <- rbind(AbFamF, family) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

rownames(AbFamF) <- c(1:nrow(AbFamF)) ### change the rownames to consecutive numbers 
AbFamF <- data.frame(Primer_name = AbFamF$Primer_name, Family = AbFamF$Family, Reads = AbFamF$Reads) ###change the order of the columns

AbFamF %>%
  group_by(Primer_name) %>% 
  mutate(Total_reads = sum(Reads)) -> AbFamF ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

Relative_abundance = AbFamF$Reads/AbFamF$Total_reads ### create a vector with the result of the operation 

AbFamF[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbFamF)[5] <- "Relative_abundance" ### Change the name of the column 

#unique(AbPhy$Primer_name)

AbFamF$Primer_name <- gsub(pattern = " ", replacement = "", x = AbFamF$Primer_name) ### check if the primer names have extra spaces
AbFamF$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbFamF$Primer_name)


#intersect(AbGenF$Primer_name, primers$Primer_name)
#setdiff(AbGenF$Primer_name, primers$Primer_name)
#setdiff(primers$Primer_name, AbGenF$Primer_name)

AbFamF <- merge(AbFamF, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

#write.csv(AbFamF, file = "/localstorage/victor/AA_Primer_evaluation/Phyloseq_Raccoon/Family_list_Racc_all.csv")


##############Plots for primers#############
###Plot raw counts by primer pair 

##Phylum
ggplot(data=AbPhyF, aes(x= reorder(Primer_name, Gen), y= Relative_abundance, fill= Phyla)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(.~Experiment) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ AbyPp$Target_gene) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=20)) 


###Genus
ggplot(data=AbGenF, aes(x= Primer_name, y= Relative_abundance, fill= Genus)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_colour_brewer(palette = "Set1") +
  #facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ AbyPp$Target_gene) +
  theme(legend.position="down") + guides(fill=guide_legend(nrow=5)) 

###Family
ggplot(data=AbFamF, aes(x= Primer_name, y= Relative_abundance, fill= Family)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_colour_brewer(palette = "Set1") +
  #facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ AbyPp$Target_gene) +
  theme(legend.position="down") + guides(fill=guide_legend(nrow=5)) 


###Bacterial biome

ggplot(data = subset(AbPhyF, Gen == "16S" & Target == "Bacteria") , aes(x = Primer_name, y = Relative_abundance, fill = Phyla)) +
  geom_bar(stat="identity", position="stack") +
  #facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Genl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) 

###18S

ggplot(data = subset(AbPhyF, Gen == "18S"), aes(x = Primer_name, y = Relative_abundance, fill = Phyla)) +
  geom_bar(stat="identity", position="stack") +
  facet_grid(.~Experiment) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Genl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=10)) 


###Plant
ggplot(data = subset(AbFamF, Target == "Viridiplantae"), aes(x = Primer_name, y = Relative_abundance, fill = Family)) +
  geom_bar(stat="identity", position="stack") +
  #facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Phyl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=25)) 

###Diet 16S
ggplot(data = subset(AbFamF,  Gen == "16S" & Target != "Bacteria"), aes(x = Primer_name, y = Relative_abundance, fill = Family)) +
  geom_bar(stat="identity", position="stack") +
  #facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Phyl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=20)) 


###COI

ggplot(data = subset(AbPhyF,  Gen == "COI"), aes(x = Primer_name, y = Relative_abundance, fill = Phyla)) +
  geom_bar(stat="identity", position="stack") +
  #facet_grid(.~Access_array_Status) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  # geom_text(aes(label = unique(Phyl16S$Primer_pair_name)), angle = 90)# +
  scale_colour_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ Phyl16S$Status) +
  theme(legend.position="right") + guides(fill=guide_legend(nrow=8)) 


######### Tabulate by specific Taxa####

annot.list <- lapply(seq_along(annot.list), function (i){
  an <- as.matrix(annot.list[[i]])
  rownames(an) <- colnames(STNC[[i]])
  an
})

names(annot.list) <- names(STNC)

phylalist <- lapply(annot.list, function (x) {
  if(nrow(x)>0){
    table(x[, "phylum"])
  }
})

lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Cestoda"))
lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Nematoda"))
lapply(annot.list, function (x) tabulate.taxa(x,  "genus", "Apicomplexa"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus",  "Platyhelminthes"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus",  "Dinoflagellate"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Streptophyta"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Ascomycota"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Chordata"))
lapply(annot.list, function (x) tabulate.taxa(x, "family", "Chordata"))
lapply(annot.list, function (x) tabulate.taxa(x, "phylum", "Ascomycota"))
lapply(annot.list, function (x) tabulate.taxa(x, "genus", "Amphibia"))



#################Phyloseq richness####################

PS.l[["mlCOIintF_XT_104_F.JgHCO2198_104_R"]]

for (z in 1:length(PS.l)){
  
  z<- "mlCOIintF_XT_104_F.JgHCO2198_104_R" 
  plot_richness(PS.l[[z]], c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson"))
  c <- transform_sample_counts(PS.l[[z]], function(otu) otu/sum(otu))
  o <- ordinate(c, method="NMDS", distance="bray")
  plot_ordination(c, o, title="Bray NMDS")
  t <- names(sort(taxa_sums(PS.l[[z]]), decreasing=TRUE))[1:60]
  p <- transform_sample_counts(PS.l[[z]], function(OTU) OTU/sum(OTU))
  p <- prune_taxa(t, p)
  plot_bar(p, fill= "genus")
}

names(PS.l)

estimate_richness(PS.lf[[1]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson"))
transform_sample_counts(PS.l[[z]], function(otu) otu/sum(otu))



lapply(PS.l, function (x) subset_taxa(x, phylum%in%c("Nematoda", "Apicomplexa","Platyhelminthes")))

for (z in 1:length(PS.l)){
  
  subset_taxa(PS.l[[2]], phylum%in%c("Nematoda", "Apicomplexa","Platyhelminthes"))
}



Parasites <- list(list())
for(i in vector_of_primers)
{
  if(i %in% names(PS.l))
  {
    id <- which(i == names(PS.l))
    finalPrime[[i]][["richness_index"]] <- data.frame(estimate_richness(PS.l[[id]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson")))
    Div <- finalPrime[[i]][["richness_index"]]$Shannon
    Rich <- log(finalPrime[[i]][["richness_index"]]$Observed + 1)
    evenness = Div/Rich
    finalPrime[[i]][["richness_index"]]$Evenness = evenness
    finalPrime[[i]][["number_of_samples"]] <- nrow(sample_data(PS.l[[id]])) ##For raccoon change 
    mean_ind <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = mean)
    std_dev <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = sd)
    ci <- qnorm(0.95)*std_dev/sqrt(finalPrime[[i]][["number_of_samples"]])
    upper <- mean_ind+ci
    lower <- mean_ind-ci
    finalPrime[[i]][["central_tendency"]] <- data.frame(Index_name = c("Observed","Chao1", "se.chao1", "Shannon", "Simpson", "InvSimpson", "Evenness"), 
                                                        Mean = mean_ind, 
                                                        SD = std_dev,
                                                        Upper = upper, 
                                                        Lower = lower)
    finalPrime[[i]][["primer_info"]] <- data.frame(primerInput[ primerInput$Primer_name %in% i, c(6:ncol(primerInput))])
  }
}
finalPrime[[1]] <- NULL


######Cumulative graph####
## All

num.taxa <-sapply(c("genus", "family", "order", "class", "phylum", "superkingdom"), function (rank){
  lapply(PS.l, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads <- unlist(lapply(PS.l, function (x) 
  sum(otu_table(x))))

PrimTax <- as.data.frame(cbind(num.taxa, num.reads))

PrimTax <- as.data.frame(apply(PrimTax, 2, unlist))
PrimTax[,8] <- row.names(PrimTax)  
colnames(PrimTax)[8] <- "Primer_name"


unique.taxa.cumsum <- function(ps.numord.l, rank){
  taxonNa <- lapply(ps.numord.l, function (x) 
    get_taxa_unique(x, taxonomic.rank = rank))
  names(taxonNa) <- names(ps.numord.l)
  unique.taxonNa.cumsum <- sapply(0:length(taxonNa),
                                  function (i) length(unique(unlist(taxonNa[0:i]))))
  return(unique.taxonNa.cumsum)
}

ps.numord.l <- PS.l[order(PrimTax$num.reads, decreasing=TRUE)]

cum.tax <- sapply(c("genus", "family", "class", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l, rank)
})

colnames(cum.tax) <- paste0("cum.", colnames(cum.tax))


cum.tax <- as.data.frame(cbind(1:nrow(cum.tax), cum.tax))
prnames <- names(ps.numord.l) ## get names from the primers
cum.tax[2:42,1] <- prnames ## add primer names 
colnames(cum.tax)[1] <- "Primer_name"

cum.plot <- ggplot(cum.tax) +
  geom_step(aes(Primer_name, cum.genus), color="red") +
  geom_text(aes(20, 1100, label="Genera"), color="red") +
  geom_step(aes(Primer_name, cum.family), color="purple") +
  geom_text(aes(20, 610, label="Families"), color="purple") +
  geom_step(aes(Primer_name, cum.order), color="blue") +
  geom_text(aes(20, 310, label="Orders"), color="blue") +
  geom_step(aes(Primer_name, cum.class), color="seagreen") +
  geom_text(aes(20, 135, label="Classes"), color="seagreen") +
  geom_step(aes(Primer_name, cum.phylum), color="green") +
  geom_text(aes(20, 63, label="Phyla"), color="green") +
  scale_y_log10("Cummulative count of taxa") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

####By marker 
###First merge cum.tax with primerInput
PrimTax$Primer_name <- gsub(pattern = " ", replacement = "", x = PrimTax$Primer_name) ### check if the primer names have extra spaces
intersect(cum.tax$Primer_name, primerInput$Primer_name) 
#setdiff(cum.tax$Primer_name, primerInput$Primer_name)

library(plyr)
PrimTax <-join(PrimTax, primerInput) ###merge the selected information with the origial data frame created 

#####18S####

Primtax.18S <- subset(PrimTax, Gen=="18S")

ps.numord.l.18S <- PS.l[c(Primtax.18S$Primer_name)][order(Primtax.18S$num.reads, decreasing=TRUE)]


cum.tax.18S <- sapply(c("genus", "family", "class", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.18S, rank)
})

colnames(cum.tax.18S) <- paste0("cum.", colnames(cum.tax.18S))


cum.tax.18S <- as.data.frame(cbind(1:nrow(cum.tax.18S), cum.tax.18S))
prnames18s <- names(ps.numord.l.18S) ## get names from the primers
#cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
colnames(cum.tax.18S)[1] <- "Primer_name"

cum.plot.18 <- ggplot(cum.tax.18S) +
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 1500, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 670, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 310, label="Orders"), color="#CC4678FF") +
  geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  geom_text(aes(20, 130, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
  geom_text(aes(20, 45, label="Phyla"), color="#F0F921FF") +
  scale_y_log10("Cummulative count of taxa") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  coord_cartesian(ylim = c(20, 2100))

ggsave("~/AA_Primer_evaluation/Cummulative_18S.pdf")

####28S####

Primtax.28S <- subset(PrimTax, Gen=="28S")

ps.numord.l.28S <- PS.l[c(Primtax.28S$Primer_name)][order(Primtax.28S$num.reads, decreasing=TRUE)]


cum.tax.28S <- sapply(c("genus", "family", "class", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.28S, rank)
})

colnames(cum.tax.28S) <- paste0("cum.", colnames(cum.tax.28S))


cum.tax.28S <- as.data.frame(cbind(1:nrow(cum.tax.28S), cum.tax.28S))
prnames28s <- names(ps.numord.l.28S) ## get names from the primers
#cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
colnames(cum.tax.28S)[1] <- "Primer_name"

cum.plot.28 <- ggplot(cum.tax.28S) +
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 1100, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 610, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 310, label="Orders"), color="#CC4678FF") +
  geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  geom_text(aes(20, 135, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
  geom_text(aes(20, 63, label="Phyla"), color="#F0F921FF") +
  scale_y_log10("Cummulative count of taxa") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

####28S and 18S####

Primtax.comb <- subset(PrimTax, Gen=="28S" | Gen== "18S")

ps.numord.l.comb <- PS.l[c(Primtax.comb$Primer_name)][order(Primtax.comb$num.reads, decreasing=TRUE)]


cum.tax.comb <- sapply(c("genus", "family", "class", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.comb, rank)
})

colnames(cum.tax.comb) <- paste0("cum.", colnames(cum.tax.comb))


cum.tax.comb <- as.data.frame(cbind(1:nrow(cum.tax.comb), cum.tax.comb))
prnames28s <- names(ps.numord.l.28S) ## get names from the primers
#cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
colnames(cum.tax.comb)[1] <- "Primer_name"

cum.plot.comb <- ggplot(cum.tax.comb) +
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 2000, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 810, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 410, label="Orders"), color="#CC4678FF") +
  geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  geom_text(aes(20, 135, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
  geom_text(aes(20, 45, label="Phyla"), color="#F0F921FF") +
  scale_y_log10("Cummulative count of taxa") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  coord_cartesian(ylim = c(20, 2100))

####28S, 18S, COI####

Primtax.comb2 <- subset(PrimTax, Gen=="28S" | Gen== "18S" | Gen== "COI")

ps.numord.l.comb2 <- PS.l[c(Primtax.comb2$Primer_name)][order(Primtax.comb2$num.reads, decreasing=TRUE)]


cum.tax.comb2 <- sapply(c("genus", "family", "class", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.comb2, rank)
})

colnames(cum.tax.comb2) <- paste0("cum.", colnames(cum.tax.comb2))


cum.tax.comb2 <- as.data.frame(cbind(1:nrow(cum.tax.comb2), cum.tax.comb2))
#prnames28s <- names(ps.numord.l.28S) ## get names from the primers
#cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
colnames(cum.tax.comb2)[1] <- "Primer_name"

cum.plot.comb2 <- ggplot(cum.tax.comb2) +
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 2000, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 900, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 360, label="Orders"), color="#CC4678FF") +
  geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  geom_text(aes(20, 135, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
  geom_text(aes(20, 45, label="Phyla"), color="#F0F921FF") +
  scale_y_log10("Cummulative count of taxa") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())+
  coord_cartesian(ylim = c(20, 2100))

ggsave("~/AA_Primer_evaluation/Cummulative_multimarker.pdf")

grid.arrange(cum.plot.18, cum.plot.comb, cum.plot.comb2, nrow= 1, ncol= 3)

####Parasites 18S+28S+COI
ps.numord.l.comb2[["ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R"]]<- NULL ### Just artropods

PM.para.comb2 <- lapply(ps.numord.l.comb2, function (x) {
  subset_taxa(x, phylum%in%c("Nematoda",
                             "Apicomplexa",
                             "Platyhelminthes"))
})

PG.para.comb2 <- lapply(PM.para.comb2, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


lapply(PG.para.comb2, function (x) {
  table(tax_table(x)[,6])
}) 


mat.para.comb2 <- lapply(PG.para.comb2, function (x) {
  as.matrix(tax_table(x))
}) 


P.comb2 <- lapply(PG.para.comb2, function (y) {
  make.names(tax_table(y)[, 6])
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
  scale_y_continuous(name = "Number of Parasites genera amplified")+
  scale_x_discrete(name = "Primer combination")


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

col_groups <- Primtax.comb2 %>%
  select("Primer_name", "Gen")

col_groups<- col_groups[-28,]

row.names(col_groups)<- colnames(P.intersection)

colour_groups <- list( Gen= c("18S"= "#440154FF", "28S"= "#21908CFF", "COI"= "#FDE725FF"))


library("pheatmap")
library("viridisLite") 
pheatmap(P.intersection, 
         color = plasma(100),
         border_color = "black",
         annotation_col = col_groups,
         #annotation_row = col_groups,
         annotation_colors = colour_groups,
         main= "Redundant genus of parasites amplified")


####Cummulative curve for parasites
num.taxa.para <-sapply(c("genus", "family", "order", "class", "phylum", "superkingdom"), function (rank){
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


cum.tax.para <- sapply(c("genus", "family", "class", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.para, rank)
})

colnames(cum.tax.para) <- paste0("cum.", colnames(cum.tax.para))

cum.tax.para <- as.data.frame(cbind(1:nrow(cum.tax.para), cum.tax.para))
cum.tax.para[2:28,7] <- PrimTax.para$Gen ## add marker amplified 
colnames(cum.tax.para)[1] <- "Primer_name"
colnames(cum.tax.para)[7] <- "Gen"

cum.plot.para <- ggplot(cum.tax.para) +
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 160, label="Genera"), color="#0D0887FF") +
  geom_point(aes(Primer_name, cum.genus, colour = Gen, size= 2))+
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))+
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 95, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 20, label="Orders"), color="#CC4678FF") +
  geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  geom_text(aes(20, 10, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
  geom_text(aes(20, 4, label="Phyla"), color="#F0F921FF") +
  scale_y_log10("Cummulative count of taxa (Parasites)") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = 'none')

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
  table(tax_table(x)[,6])
}) 


mat.fungi.comb2 <- lapply(PG.fungi.comb2, function (x) {
  as.matrix(tax_table(x))
}) 


F.comb2 <- lapply(PG.fungi.comb2, function (y) {
  make.names(tax_table(y)[, 6])
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


fungiprimer <- ggplot(Fungi.comb2, aes(Primer_name, color=Gen, fill=Gen))+
  geom_bar()+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(name = "Number of Fungi genera amplified")+
  scale_x_discrete(name = "Primer combination")


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

ParFun <- merge(Count.fungi.comb2, Count.para.comb2, by= "Primer_name")

ParFun <- merge(ParFun, Primtax.comb2, by= "Primer_name")

###Plot
ggplot(data = ParFun, aes(x = num.reads, y = Count_fungi, color= Gen)) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_y_continuous(name = "Number of Fungi genera")+
  scale_x_log10(name = "Expected amplicon size")+
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position="right") #+
  #geom_smooth(method = "lm", col = "black")
  

grid.arrange(paraprimer, fungiprimer, nrow= 1, ncol= 2)

####compair all the primers for multiple intersections

F.intersection <- sapply(seq_len(length(F.comb2)), function(x)
  sapply(seq_len(length(F.comb2)), function(y) length(intersect(unlist(F.comb2[x]), unlist(F.comb2[y]))))
)
rownames(F.intersection)<- names(F.comb2)
colnames(F.intersection)<- names(F.comb2)

col_groups <- Primtax.comb2 %>%
  select("Primer_name", "Gen")

col_groups<- col_groups[-28,]

row.names(col_groups)<- colnames(F.intersection)

colour_groups <- list( Gen= c("18S"= "#440154FF", "28S"= "#21908CFF", "COI"= "#FDE725FF"))


library("pheatmap")
library("viridisLite") 
pheatmap(log10(F.intersection), 
         color = plasma(100),
         border_color = "black",
         annotation_col = col_groups,
         #annotation_row = col_groups,
         annotation_colors = colour_groups,
         main= "log10 Redundant genus of Fungi amplified")


####Cummulative curve for fungi
num.taxa.fungi <-sapply(c("genus", "family", "order", "class", "phylum", "superkingdom"), function (rank){
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


cum.tax.fungi <- sapply(c("genus", "family", "class", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.fungi, rank)
})

colnames(cum.tax.fungi) <- paste0("cum.", colnames(cum.tax.fungi))


cum.tax.fungi <- as.data.frame(cbind(1:nrow(cum.tax.fungi), cum.tax.fungi))
cum.tax.fungi[2:28,7] <- PrimTax.fungi$Gen ## add marker amplified 
colnames(cum.tax.fungi)[1] <- "Primer_name"
colnames(cum.tax.fungi)[7] <- "Gen"

cum.plot.fungi <- ggplot(cum.tax.fungi) +
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 770, label="Genera"), color="#0D0887FF") +
  geom_point(aes(Primer_name, cum.genus, colour = Gen, size= 2))+
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))+
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 310, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 120, label="Orders"), color="#CC4678FF") +
  geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  geom_text(aes(20, 50, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
  geom_text(aes(20, 8, label="Phyla"), color="#F0F921FF") +
  scale_y_log10("Cummulative count of taxa (Fungi)") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
  

grid.arrange(cum.plot.para, cum.plot.fungi, nrow= 1, ncol= 2) ###Put parasites and fungi together 

####Bacterial biome#####
Primtax.16S <- subset(PrimTax, Gen=="16S" & Target=="Bacteria")

ps.numord.l.16S <- PS.l[c(Primtax.16S$Primer_name)][order(Primtax.16S$num.reads, decreasing=TRUE)]


cum.tax.16S <- sapply(c("genus", "family", "class", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.16S, rank)
})

colnames(cum.tax.16S) <- paste0("cum.", colnames(cum.tax.16S))


cum.tax.16S <- as.data.frame(cbind(1:nrow(cum.tax.16S), cum.tax.16S))
prnames16s <- names(ps.numord.l.16S) ## get names from the primers
#cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
colnames(cum.tax.16S)[1] <- "Primer_name"

cum.plot.16 <- ggplot(cum.tax.16S) +
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(7, 540, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(7, 255, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(7, 135, label="Orders"), color="#CC4678FF") +
  geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  geom_text(aes(7, 70, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
  geom_text(aes(7, 34, label="Phyla"), color="#F0F921FF") +
  scale_y_log10("Cummulative count of taxa (Bacteria)") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

####Redundancy bacterial biome

PG.Bac <- lapply(ps.numord.l.16S , function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


lapply(PG.Bac, function (x) {
  table(tax_table(x)[,6])
}) 


mat.Bac <- lapply(PG.Bac, function (x) {
  as.matrix(tax_table(x))
}) 


Bac.comb <- lapply(PG.Bac, function (y) {
  make.names(tax_table(y)[, 6])
}) 

Bacterial.comb <- data.frame() ###Create the data frame 

for (i in 1: length(Bac.comb)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(Bac.comb[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((Bac.comb[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(Bac.comb)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Genus", "Count", "Primer_name") ### change the names for the columns 
  Bacterial.comb <- rbind(Bacterial.comb, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

Bacterial.comb <- unique(Bacterial.comb[c("Genus", "Primer_name")]) ##Take unique combinations 
Bacterial.comb <- join(Bacterial.comb, PrimTax) 


Bacteriaprimer <- ggplot(Bacterial.comb, aes(Primer_name, color=Gen, fill=Gen))+
  geom_bar()+
  coord_flip()+
  theme_classic()+
  scale_y_continuous(name = "Number of Bacteria genera amplified")+
  scale_x_discrete(name = "Primer combination")


Bacteria <- Bacterial.comb$Genus
Bacteria <- as.data.frame(table(Bacteria))
Bacteria <- Bacteria[order(-Bacteria$Freq),]
###Top ten Bacteria 
head(Bacteria, 10)

B.intersection <- sapply(seq_len(length(Bac.comb)), function(x)
  sapply(seq_len(length(Bac.comb)), function(y) length(intersect(unlist(Bac.comb[x]), unlist(Bac.comb[y]))))
)
rownames(B.intersection)<- names(Bac.comb)
colnames(B.intersection)<- names(Bac.comb)

col_groups.bac <- Primtax.16S %>%
  select("Primer_name", "Expected")

row.names(col_groups.bac)<- colnames(B.intersection)

library("pheatmap")
library("viridisLite") 
pheatmap(B.intersection, 
         color = plasma(100),
         border_color = "black",
         annotation_col = col_groups.bac,
         show_rownames = F,
         show_colnames = F,
         #annotation_row = Primer_name,
         #annotation_colors = colour_groups,
         main= "log10 Redundant genus of Bacteria amplified")




###Parasites 18S####

PM.para.18S <- lapply(ps.numord.l.18S, function (x) {
                subset_taxa(x, phylum%in%c("Nematoda",
                              "Apicomplexa",
                              "Platyhelminthes"))
})

PG.para.18S <- lapply(PM.para.18S, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})
  
  
lapply(PG.para.18S, function (x) {
  table(tax_table(x)[,6])
}) 


mat.18S <- lapply(PG.para.18S, function (x) {
            as.matrix(tax_table(x))
}) 


P.18S <- lapply(PG.para.18S, function (y) {
    make.names(tax_table(y)[, 6])
}) 


Parasites.18S <- data.frame() ###Create the data frame 

for (i in 1: length(P.18S)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(P.18S[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((P.18S[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(P.18S)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Genus", "Count", "Primer_name") ### change the names for the columns 
  Parasites.18S <- rbind(Parasites.18S, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop



Parasites <- Parasites.18S$Genus
Parasites <- as.data.frame(table(Parasites))
Parasites <- Parasites[order(-Parasites$Freq),]


###Genus of parasites in 18S

Count.18S <- data.frame() ###Create the data frame 

for (i in 1: length(P.18S)) ### Start a loop: fro every element in the list ...
{ 
  cnt <- data.frame() #### make an individual data frame ...
  
    cnt <- as.data.frame(length(unique(P.18S[[i]])))  ###Make a data frame with the data included in each element of the list 
  cnt[,2] <- names(P.18S[i]) ### Take the names of every list and use them to fill column 2 as many times the logitude of the column 1
  colnames(cnt) <- c("Count", "Primer_name") ### change the names for the columns 
  Count.18S <- rbind(Count.18S, cnt) ### Join all the "individual" data frames into the final data frame 
} 

###Fungi 18S

PM.fungi.18S <- lapply(ps.numord.l.18S, function (x) {
  subset_taxa(x, phylum%in%c("Ascomycota",
                             "Basidiomycota",
                             "Zygomycota"))
})

PG.fungi.18S <- lapply(PM.fungi.18S, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


lapply(PG.fungi.18S, function (x) {
  table(tax_table(x)[,6])
}) 


mat.fungi.18S <- lapply(PG.fungi.18S, function (x) {
  as.matrix(tax_table(x))
}) 


F.18S <- lapply(PG.fungi.18S, function (y) {
  make.names(tax_table(y)[, 6])
}) 


Count.fungi.18s <- data.frame()

for (i in 1: length(P.18S)) ### Start a loop: fro every element in the list ...
{ 
  cnt <- data.frame() #### make an individual data frame ...
  
    cnt <- as.data.frame(length(unique(F.18S[[i]])))  ###Make a data frame with the data included in each element of the list 
  cnt[,2] <- names(F.18S[i]) ### Take the names of every list and use them to fill column 2 as many times the logitude of the column 1
  colnames(cnt) <- c("Count", "Primer_name") ### change the names for the columns 
  Count.fungi.18s <- rbind(Count.fungi.18s, cnt) ### Join all the "individual" data frames into the final data frame 
} 



####Parasites 18S+28S####

PM.para.comb <- lapply(ps.numord.l.comb, function (x) {
  subset_taxa(x, phylum%in%c("Nematoda",
                             "Apicomplexa",
                             "Platyhelminthes"))
})

PG.para.comb <- lapply(PM.para.comb, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


lapply(PG.para.comb, function (x) {
  table(tax_table(x)[,6])
}) 


mat.para.comb <- lapply(PG.para.comb, function (x) {
  as.matrix(tax_table(x))
}) 


P.comb <- lapply(PG.para.comb, function (y) {
  make.names(tax_table(y)[, 6])
}) 



Parasites.comb <- data.frame() ###Create the data frame 

for (i in 1: length(P.comb)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(P.comb[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((P.comb[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(P.comb)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Genus", "Count", "Primer_name") ### change the names for the columns 
  Parasites.comb <- rbind(Parasites.comb, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

Parasites.comb <- unique(Parasites.comb[c("Genus", "Primer_name")]) ##Take unique combinations 

Parasites <- Parasites.comb$Genus
Parasites <- as.data.frame(table(Parasites))
Parasites <- Parasites[order(-Parasites$Freq),]


###Add a bar plot (To improve)

ggplot(Parasites, aes(reorder(Parasites, Freq), Freq))+
  geom_point()+
  coord_flip()

###Genus of parasites in 18S and 28S

Count.para.comb <- data.frame() ###Create the data frame 

for (i in 1: length(P.comb)) ### Start a loop: fro every element in the list ...
{ 
  cnt <- data.frame() #### make an individual data frame ...
  
  cnt <- as.data.frame(length(unique(P.comb[[i]])))  ###Make a data frame with the data included in each element of the list 
  cnt[,2] <- names(P.comb[i]) ### Take the names of every list and use them to fill column 2 as many times the logitude of the column 1
  colnames(cnt) <- c("Count_para", "Primer_name") ### change the names for the columns 
  Count.para.comb <- rbind(Count.para.comb, cnt) ### Join all the "individual" data frames into the final data frame 
} 

###Fungi 18s + 28S

PM.fungi.comb <- lapply(ps.numord.l.comb, function (x) {
  subset_taxa(x, phylum%in%c("Ascomycota",
                             "Basidiomycota",
                             "Zygomycota", "Mucoromycota", "Zoopagomycota", "Blastocladiomycota", "Cryptomycota", "Microsporidia"))
})

PG.fungi.comb <- lapply(PM.fungi.comb, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


lapply(PG.fungi.comb, function (x) {
  table(tax_table(x)[,6])
}) 


mat.fungi.comb <- lapply(PG.fungi.comb, function (x) {
  as.matrix(tax_table(x))
}) 


F.comb <- lapply(PG.fungi.comb, function (y) {
  make.names(tax_table(y)[, 6])
}) 

Fungi.comb <- data.frame() ###Create the data frame 

for (i in 1: length(F.comb)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(F.comb[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((F.comb[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(F.comb)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Genus", "Count", "Primer_name") ### change the names for the columns 
  Fungi.comb <- rbind(Fungi.comb, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

Fungi.comb <- unique(Fungi.comb[c("Genus", "Primer_name")]) ##Take unique combinations 

Fungi <- Fungi.comb$Genus
Fungi <- as.data.frame(table(Fungi))
Fungi <- Fungi[order(-Fungi$Freq),]

###Add a bar plot (To improve)

ggplot(Fungi, aes(reorder(Fungi, Freq), Freq))+
  geom_point()+
  coord_flip()


Count.fungi.comb <- data.frame()

for (i in 1: length(F.comb)) ### Start a loop: fro every element in the list ...
{ 
  cnt <- data.frame() #### make an individual data frame ...
  
  cnt <- as.data.frame(length(unique(F.comb[[i]])))  ###Make a data frame with the data included in each element of the list 
  cnt[,2] <- names(F.comb[i]) ### Take the names of every list and use them to fill column 2 as many times the logitude of the column 1
  colnames(cnt) <- c("Count_fungi", "Primer_name") ### change the names for the columns 
  Count.fungi.comb <- rbind(Count.fungi.comb, cnt) ### Join all the "individual" data frames into the final data frame 
} 


ParFun <- merge(Count.fungi.comb, Count.para.comb, by= "Primer_name")

ParFun <- merge(ParFun, Primtax.comb, by= "Primer_name")


###What is not Parasite or Fungi???

PM.other.comb <- lapply(ps.numord.l.comb, function (x) {
  subset_taxa(x, !(phylum%in%c("Ascomycota",
                             "Basidiomycota",
                             "Zygomycota",
                             "Nematoda",
                             "Apicomplexa",
                             "Platyhelminthes")))
})

PG.other.comb <- lapply(PM.other.comb, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})

lapply(PG.other.comb, function (x) {
  table(tax_table(x)[,2])
}) 


######COI####
Primtax.coi <- subset(PrimTax, Gen== "COI")

ps.numord.l.coi <- PS.l[c(Primtax.coi$Primer_name)][order(Primtax.coi$num.reads, decreasing=TRUE)]


cum.tax.coi <- sapply(c("genus", "family", "class", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.coi, rank)
})

colnames(cum.tax.coi) <- paste0("cum.", colnames(cum.tax.coi))


cum.tax.coi <- as.data.frame(cbind(1:nrow(cum.tax.coi), cum.tax.coi))
#prnames28s <- names(ps.numord.l.28S) ## get names from the primers
#cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
colnames(cum.tax.coi)[1] <- "Primer_name"

cum.plot.comb2 <- ggplot(cum.tax.comb2) +
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 1100, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 610, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 310, label="Orders"), color="#CC4678FF") +
  geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  geom_text(aes(20, 135, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F0F921FF") +
  geom_text(aes(20, 63, label="Phyla"), color="#F0F921FF") +
  scale_y_log10("Cummulative count of taxa") +
  scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

####Parasites COI
ps.numord.l.coi[["ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R"]]<- NULL

PM.para.coi <- lapply(ps.numord.l.coi, function (x) {
  subset_taxa(x, phylum%in%c("Nematoda",
                             "Apicomplexa",
                             "Platyhelminthes"))
})

PG.para.coi <- lapply(PM.para.coi, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


lapply(PG.para.coi, function (x) {
  table(tax_table(x)[,6])
}) 


mat.para.coi <- lapply(PG.para.coi, function (x) {
  as.matrix(tax_table(x))
}) 


P.coi <- lapply(PG.para.coi, function (y) {
  make.names(tax_table(y)[, 6])
}) 



Parasites.coi <- data.frame() ###Create the data frame 

for (i in 1: length(P.coi)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(P.coi[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((P.coi[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(P.coi)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Genus", "Count", "Primer_name") ### change the names for the columns 
  Parasites.coi <- rbind(Parasites.coi, genus) ### Join all the "individual" data frames into the final data frame 
  
} 


Parasites.coi <- unique(Parasites.coi[c("Genus", "Primer_name")]) ##Take unique combinations 

Parasites <- Parasites.coi$Genus
Parasites <- as.data.frame(table(Parasites))
Parasites <- Parasites[order(-Parasites$Freq),]


###Add a bar plot (To improve)

ggplot(Parasites, aes(reorder(Parasites, Freq), Freq))+
  geom_point()+
  coord_flip()

###Genus of parasites in 18S and 28S

Count.para.comb <- data.frame() ###Create the data frame 

for (i in 1: length(P.comb)) ### Start a loop: fro every element in the list ...
{ 
  cnt <- data.frame() #### make an individual data frame ...
  
  cnt <- as.data.frame(length(unique(P.comb[[i]])))  ###Make a data frame with the data included in each element of the list 
  cnt[,2] <- names(P.comb[i]) ### Take the names of every list and use them to fill column 2 as many times the logitude of the column 1
  colnames(cnt) <- c("Count_para", "Primer_name") ### change the names for the columns 
  Count.para.comb <- rbind(Count.para.comb, cnt) ### Join all the "individual" data frames into the final data frame 
} 

###Fungi 18s + 28S

PM.fungi.comb <- lapply(ps.numord.l.comb, function (x) {
  subset_taxa(x, phylum%in%c("Ascomycota",
                             "Basidiomycota",
                             "Zygomycota", "Mucoromycota", "Zoopagomycota", "Blastocladiomycota", "Cryptomycota", "Microsporidia"))
})

PG.fungi.comb <- lapply(PM.fungi.comb, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


lapply(PG.fungi.comb, function (x) {
  table(tax_table(x)[,6])
}) 


mat.fungi.comb <- lapply(PG.fungi.comb, function (x) {
  as.matrix(tax_table(x))
}) 


F.comb <- lapply(PG.fungi.comb, function (y) {
  make.names(tax_table(y)[, 6])
}) 

Fungi.comb <- data.frame() ###Create the data frame 

for (i in 1: length(F.comb)) ### Start a loop: fro every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(F.comb[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((F.comb[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(F.comb)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("Genus", "Count", "Primer_name") ### change the names for the columns 
  Fungi.comb <- rbind(Fungi.comb, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

Fungi.comb <- unique(Fungi.comb[c("Genus", "Primer_name")]) ##Take unique combinations 

Fungi <- Fungi.comb$Genus
Fungi <- as.data.frame(table(Fungi))
Fungi <- Fungi[order(-Fungi$Freq),]

###Add a bar plot (To improve)

ggplot(Fungi, aes(reorder(Fungi, Freq), Freq))+
  geom_point()+
  coord_flip()


Count.fungi.comb <- data.frame()

for (i in 1: length(F.comb)) ### Start a loop: fro every element in the list ...
{ 
  cnt <- data.frame() #### make an individual data frame ...
  
  cnt <- as.data.frame(length(unique(F.comb[[i]])))  ###Make a data frame with the data included in each element of the list 
  cnt[,2] <- names(F.comb[i]) ### Take the names of every list and use them to fill column 2 as many times the logitude of the column 1
  colnames(cnt) <- c("Count_fungi", "Primer_name") ### change the names for the columns 
  Count.fungi.comb <- rbind(Count.fungi.comb, cnt) ### Join all the "individual" data frames into the final data frame 
} 
