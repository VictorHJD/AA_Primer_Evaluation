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
  join(primerInput, by= "Primer_comb_ID")->seqcounts

seqcounts[,1]<-NULL
seqcounts[,1]<-NULL

seqcounts%>%
  select(c("Primer_comb_ID", "Seq_F", "Seq_R", "Gen"))-> Primers

##PrimerTree pipeline
##Considering just the first 1000 alignments
##For degenerated primers uo to 25 combinations are tested for now
##18S primer pairs

Primers %>%
  filter(Gen== "18S")-> Primers_18S

Euk_18S_01<- search_primer_pair(name= as.character(Primers_18S[1,1]), as.character(Primers_18S[1,2]), as.character(Primers_18S[1,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_01, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_01.Rds")

Euk_18S_02<- search_primer_pair(name= as.character(Primers_18S[2,1]), as.character(Primers_18S[2,2]), as.character(Primers_18S[2,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_02, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_02.Rds")

Euk_18S_03<- search_primer_pair(name= as.character(Primers_18S[3,1]), as.character(Primers_18S[3,2]), as.character(Primers_18S[3,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_03, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_03.Rds")

Euk_18S_04<- search_primer_pair(name= as.character(Primers_18S[4,1]), as.character(Primers_18S[4,2]), as.character(Primers_18S[4,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_04, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_04.Rds")

Euk_18S_05<- search_primer_pair(name= as.character(Primers_18S[5,1]), as.character(Primers_18S[5,2]), as.character(Primers_18S[5,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_05, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_05.Rds")

Euk_18S_06<- search_primer_pair(name= as.character(Primers_18S[6,1]), as.character(Primers_18S[6,2]), as.character(Primers_18S[6,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_06, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_06.Rds")

Euk_18S_07<- search_primer_pair(name= as.character(Primers_18S[7,1]), as.character(Primers_18S[7,2]), as.character(Primers_18S[7,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_07, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_07.Rds")

Euk_18S_08<- search_primer_pair(name= as.character(Primers_18S[8,1]), as.character(Primers_18S[8,2]), as.character(Primers_18S[8,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_08, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_08.Rds")

Euk_18S_09<- search_primer_pair(name= as.character(Primers_18S[9,1]), as.character(Primers_18S[9,2]), as.character(Primers_18S[9,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_09, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_09.Rds")

Euk_18S_10<- search_primer_pair(name= as.character(Primers_18S[10,1]), as.character(Primers_18S[10,2]), as.character(Primers_18S[10,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_10, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_10.Rds")

Euk_18S_11<- search_primer_pair(name= as.character(Primers_18S[11,1]), as.character(Primers_18S[11,2]), as.character(Primers_18S[11,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_11, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_11.Rds")

Euk_18S_12<- search_primer_pair(name= as.character(Primers_18S[12,1]), as.character(Primers_18S[12,2]), as.character(Primers_18S[12,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_12, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_12.Rds")

Euk_18S_13<- search_primer_pair(name= as.character(Primers_18S[13,1]), as.character(Primers_18S[13,2]), as.character(Primers_18S[13,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_13, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_13.Rds")

Euk_18S_14<- search_primer_pair(name= as.character(Primers_18S[14,1]), as.character(Primers_18S[14,2]), as.character(Primers_18S[14,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_14, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_14.Rds")

Euk_18S_15<- search_primer_pair(name= as.character(Primers_18S[15,1]), as.character(Primers_18S[15,2]), as.character(Primers_18S[15,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_15, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_15.Rds")

Euk_18S_16<- search_primer_pair(name= as.character(Primers_18S[16,1]), as.character(Primers_18S[16,2]), as.character(Primers_18S[16,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_16, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_16.Rds")

Euk_18S_17<- search_primer_pair(name= as.character(Primers_18S[17,1]), as.character(Primers_18S[17,2]), as.character(Primers_18S[17,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_17, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_17.Rds")

Euk_18S_18<- search_primer_pair(name= as.character(Primers_18S[18,1]), as.character(Primers_18S[18,2]), as.character(Primers_18S[18,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_18, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_18.Rds")

Euk_18S_19<- search_primer_pair(name=as.character(Primers_18S[19,1]), as.character(Primers_18S[19,2]), as.character(Primers_18S[19,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_19, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_19.Rds")

Euk_18S_20<- search_primer_pair(name= as.character(Primers_18S[20,1]), as.character(Primers_18S[20,2]), as.character(Primers_18S[20,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_20, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_20.Rds")

Euk_18S_21<- search_primer_pair(name= as.character(Primers_18S[21,1]), as.character(Primers_18S[21,2]), as.character(Primers_18S[21,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_21, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_21.Rds")

#Make a list of primerTree objects saved
done18S <- list.files(path = "~/AA_Primer_evaluation/output/primerTreeObj/", pattern = "Euk_18S_")
done18S<- lapply(done18S, function(x) gsub(".Rds", "", x))

##General overview of the results
plot(Euk_18S_01, ranks= "phylum")

##Extract the relevant information
Euk_18S_01_Tax<- Euk_18S_01$taxonomy ##Taxonomy information 
Euk_18S_01_BLAST<- Euk_18S_01$BLAST_result #BLAST information

Euk_18S_01_Tax[is.na(Euk_18S_01_Tax)]<- "Unassigned" ##Change NA's into Unassigned 

Euk_18S_01_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_01_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_01_Results

hist(Euk_18S_01_Results$product_length)

Euk_18S_01_RA<- data.frame(count(Euk_18S_01_Results$phylum))

Euk_18S_01_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_01")-> Euk_18S_01_RA ##Data frame with all the relative abundance in silico per primer pair

colnames(Euk_18S_01_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_02
plot(Euk_18S_02, ranks= "phylum")
Euk_18S_02_Tax<- Euk_18S_02$taxonomy ##Taxonomy information 
Euk_18S_02_BLAST<- Euk_18S_02$BLAST_result #BLAST information
Euk_18S_02_Tax[is.na(Euk_18S_02_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_02_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_02_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_02_Results
hist(Euk_18S_02_Results$product_length)
Euk_18S_02_RA<- data.frame(count(Euk_18S_02_Results$phylum))
Euk_18S_02_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_02")->Euk_18S_02_RA

colnames(Euk_18S_02_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_03
plot(Euk_18S_03, ranks= "phylum")
Euk_18S_03_Tax<- Euk_18S_03$taxonomy ##Taxonomy information 
Euk_18S_03_BLAST<- Euk_18S_03$BLAST_result #BLAST information
Euk_18S_03_Tax[is.na(Euk_18S_03_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_03_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_03_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_03_Results
hist(Euk_18S_03_Results$product_length)
Euk_18S_03_RA<- data.frame(count(Euk_18S_03_Results$phylum))
Euk_18S_03_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_03")->Euk_18S_03_RA

colnames(Euk_18S_03_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_04
plot(Euk_18S_04, ranks= "phylum")
Euk_18S_04_Tax<- Euk_18S_04$taxonomy ##Taxonomy information 
Euk_18S_04_BLAST<- Euk_18S_04$BLAST_result #BLAST information
Euk_18S_04_Tax[is.na(Euk_18S_04_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_04_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_04_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_04_Results
hist(Euk_18S_04_Results$product_length)
Euk_18S_04_RA<- data.frame(count(Euk_18S_04_Results$phylum))
Euk_18S_04_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_04")->Euk_18S_04_RA

colnames(Euk_18S_04_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_05
plot(Euk_18S_05, ranks= "phylum")
Euk_18S_05_Tax<- Euk_18S_05$taxonomy ##Taxonomy information 
Euk_18S_05_BLAST<- Euk_18S_05$BLAST_result #BLAST information
Euk_18S_05_Tax[is.na(Euk_18S_05_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_05_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_05_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_05_Results
hist(Euk_18S_05_Results$product_length)
Euk_18S_05_RA<- data.frame(count(Euk_18S_05_Results$phylum))
Euk_18S_05_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_05")->Euk_18S_05_RA

colnames(Euk_18S_05_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_06
plot(Euk_18S_06, ranks= "phylum")
Euk_18S_06_Tax<- Euk_18S_06$taxonomy ##Taxonomy information 
Euk_18S_06_BLAST<- Euk_18S_06$BLAST_result #BLAST information
Euk_18S_06_Tax[is.na(Euk_18S_06_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_06_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_06_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_06_Results
hist(Euk_18S_06_Results$product_length)
Euk_18S_06_RA<- data.frame(count(Euk_18S_06_Results$phylum))
Euk_18S_06_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_06")->Euk_18S_06_RA

colnames(Euk_18S_06_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_07
plot(Euk_18S_07, ranks= "phylum")
Euk_18S_07_Tax<- Euk_18S_07$taxonomy ##Taxonomy information 
Euk_18S_07_BLAST<- Euk_18S_07$BLAST_result #BLAST information
Euk_18S_07_Tax[is.na(Euk_18S_07_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_07_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_07_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_07_Results
hist(Euk_18S_07_Results$product_length)
Euk_18S_07_RA<- data.frame(count(Euk_18S_07_Results$phylum))
Euk_18S_07_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_07")->Euk_18S_07_RA

colnames(Euk_18S_07_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_08
plot(Euk_18S_08, ranks= "phylum")
Euk_18S_08_Tax<- Euk_18S_08$taxonomy ##Taxonomy information 
Euk_18S_08_BLAST<- Euk_18S_08$BLAST_result #BLAST information
Euk_18S_08_Tax[is.na(Euk_18S_08_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_08_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_08_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_08_Results
hist(Euk_18S_08_Results$product_length)
Euk_18S_08_RA<- data.frame(count(Euk_18S_08_Results$phylum))
Euk_18S_08_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_08")->Euk_18S_08_RA

colnames(Euk_18S_08_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_09
plot(Euk_18S_09, ranks= "phylum")
Euk_18S_09_Tax<- Euk_18S_09$taxonomy ##Taxonomy information 
Euk_18S_09_BLAST<- Euk_18S_09$BLAST_result #BLAST information
Euk_18S_09_Tax[is.na(Euk_18S_09_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_09_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_09_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_09_Results
hist(Euk_18S_09_Results$product_length)
Euk_18S_09_RA<- data.frame(count(Euk_18S_09_Results$phylum))
Euk_18S_09_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_09")->Euk_18S_09_RA

colnames(Euk_18S_09_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_10
plot(Euk_18S_10, ranks= "phylum")
Euk_18S_10_Tax<- Euk_18S_10$taxonomy ##Taxonomy information 
Euk_18S_10_BLAST<- Euk_18S_10$BLAST_result #BLAST information
Euk_18S_10_Tax[is.na(Euk_18S_10_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_10_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_10_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_10_Results
hist(Euk_18S_10_Results$product_length)
Euk_18S_10_RA<- data.frame(count(Euk_18S_10_Results$phylum))
Euk_18S_10_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_10")->Euk_18S_10_RA

colnames(Euk_18S_10_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_11
plot(Euk_18S_11, ranks= "phylum")
Euk_18S_11_Tax<- Euk_18S_11$taxonomy ##Taxonomy information 
Euk_18S_11_BLAST<- Euk_18S_11$BLAST_result #BLAST information
Euk_18S_11_Tax[is.na(Euk_18S_11_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_11_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_11_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_11_Results
hist(Euk_18S_11_Results$product_length)
Euk_18S_11_RA<- data.frame(count(Euk_18S_11_Results$phylum))
Euk_18S_11_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_11")->Euk_18S_11_RA

colnames(Euk_18S_11_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_12
plot(Euk_18S_12, ranks= "phylum")
Euk_18S_12_Tax<- Euk_18S_12$taxonomy ##Taxonomy information 
Euk_18S_12_BLAST<- Euk_18S_12$BLAST_result #BLAST information
Euk_18S_12_Tax[is.na(Euk_18S_12_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_12_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_12_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_12_Results
hist(Euk_18S_12_Results$product_length)
Euk_18S_12_RA<- data.frame(count(Euk_18S_12_Results$phylum))
Euk_18S_12_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_12")->Euk_18S_12_RA

colnames(Euk_18S_12_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_13
plot(Euk_18S_13, ranks= "phylum")
Euk_18S_13_Tax<- Euk_18S_13$taxonomy ##Taxonomy information 
Euk_18S_13_BLAST<- Euk_18S_13$BLAST_result #BLAST information
Euk_18S_13_Tax[is.na(Euk_18S_13_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_13_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_13_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_13_Results
hist(Euk_18S_13_Results$product_length)
Euk_18S_13_RA<- data.frame(count(Euk_18S_13_Results$phylum))
Euk_18S_13_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_13")->Euk_18S_13_RA

colnames(Euk_18S_13_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_14
plot(Euk_18S_14, ranks= "phylum")
Euk_18S_14_Tax<- Euk_18S_14$taxonomy ##Taxonomy information 
Euk_18S_14_BLAST<- Euk_18S_14$BLAST_result #BLAST information
Euk_18S_14_Tax[is.na(Euk_18S_14_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_14_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_14_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_14_Results
hist(Euk_18S_14_Results$product_length)
Euk_18S_14_RA<- data.frame(count(Euk_18S_14_Results$phylum))
Euk_18S_14_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_14")->Euk_18S_14_RA

colnames(Euk_18S_14_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_15
plot(Euk_18S_15, ranks= "phylum")
Euk_18S_15_Tax<- Euk_18S_15$taxonomy ##Taxonomy information 
Euk_18S_15_BLAST<- Euk_18S_15$BLAST_result #BLAST information
Euk_18S_15_Tax[is.na(Euk_18S_15_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_15_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_15_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_15_Results
hist(Euk_18S_15_Results$product_length)
Euk_18S_15_RA<- data.frame(count(Euk_18S_15_Results$phylum))
Euk_18S_15_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_15")->Euk_18S_15_RA

colnames(Euk_18S_15_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_16
plot(Euk_18S_16, ranks= "phylum")
Euk_18S_16_Tax<- Euk_18S_16$taxonomy ##Taxonomy information 
Euk_18S_16_BLAST<- Euk_18S_16$BLAST_result #BLAST information
Euk_18S_16_Tax[is.na(Euk_18S_16_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_16_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_16_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_16_Results
hist(Euk_18S_16_Results$product_length)
Euk_18S_16_RA<- data.frame(count(Euk_18S_16_Results$phylum))
Euk_18S_16_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_16")->Euk_18S_16_RA

colnames(Euk_18S_16_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_17
plot(Euk_18S_17, ranks= "phylum")
Euk_18S_17_Tax<- Euk_18S_17$taxonomy ##Taxonomy information 
Euk_18S_17_BLAST<- Euk_18S_17$BLAST_result #BLAST information
Euk_18S_17_Tax[is.na(Euk_18S_17_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_17_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_17_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_17_Results
hist(Euk_18S_17_Results$product_length)
Euk_18S_17_RA<- data.frame(count(Euk_18S_17_Results$phylum))
Euk_18S_17_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_17")->Euk_18S_17_RA

colnames(Euk_18S_17_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_18
plot(Euk_18S_18, ranks= "phylum")
Euk_18S_18_Tax<- Euk_18S_18$taxonomy ##Taxonomy information 
Euk_18S_18_BLAST<- Euk_18S_18$BLAST_result #BLAST information
Euk_18S_18_Tax[is.na(Euk_18S_18_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_18_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_18_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_18_Results
hist(Euk_18S_18_Results$product_length)
Euk_18S_18_RA<- data.frame(count(Euk_18S_18_Results$phylum))
Euk_18S_18_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_18")->Euk_18S_18_RA

colnames(Euk_18S_18_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_19
plot(Euk_18S_19, ranks= "phylum")
Euk_18S_19_Tax<- Euk_18S_19$taxonomy ##Taxonomy information 
Euk_18S_19_BLAST<- Euk_18S_19$BLAST_result #BLAST information
Euk_18S_19_Tax[is.na(Euk_18S_19_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_19_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_19_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_19_Results
hist(Euk_18S_19_Results$product_length)
Euk_18S_19_RA<- data.frame(count(Euk_18S_19_Results$phylum))
Euk_18S_19_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_19")->Euk_18S_19_RA

colnames(Euk_18S_19_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_20
plot(Euk_18S_20, ranks= "phylum")
Euk_18S_20_Tax<- Euk_18S_20$taxonomy ##Taxonomy information 
Euk_18S_20_BLAST<- Euk_18S_20$BLAST_result #BLAST information
Euk_18S_20_Tax[is.na(Euk_18S_20_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_20_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_20_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_20_Results
hist(Euk_18S_20_Results$product_length)
Euk_18S_20_RA<- data.frame(count(Euk_18S_20_Results$phylum))
Euk_18S_20_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_20")->Euk_18S_20_RA

colnames(Euk_18S_20_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_18S_21
plot(Euk_18S_21, ranks= "phylum")
Euk_18S_21_Tax<- Euk_18S_21$taxonomy ##Taxonomy information 
Euk_18S_21_BLAST<- Euk_18S_21$BLAST_result #BLAST information
Euk_18S_21_Tax[is.na(Euk_18S_21_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_18S_21_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_18S_21_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_21_Results
hist(Euk_18S_21_Results$product_length)
Euk_18S_21_RA<- data.frame(count(Euk_18S_21_Results$phylum))
Euk_18S_21_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_21")->Euk_18S_21_RA

colnames(Euk_18S_21_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

##Relative abundance in silico comparison

Relative_abundance_18S<- bind_rows(Euk_18S_01_RA, Euk_18S_02_RA, Euk_18S_03_RA, Euk_18S_04_RA, Euk_18S_05_RA, Euk_18S_06_RA,
                               Euk_18S_07_RA, Euk_18S_08_RA, Euk_18S_09_RA, Euk_18S_10_RA, Euk_18S_11_RA, Euk_18S_12_RA,
                               Euk_18S_13_RA, Euk_18S_14_RA, Euk_18S_15_RA, Euk_18S_16_RA, Euk_18S_17_RA, Euk_18S_18_RA,
                               Euk_18S_19_RA, Euk_18S_20_RA, Euk_18S_21_RA)

Relative_abundance_18S%>%
  group_by(Primer_name)%>%
  mutate(Main_taxa= Rel_abund>= 0.05) %>%
  mutate(phylum= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(.$phylum))) %>%
  arrange(Primer_name, desc(phylum))->Relative_abundance_18S

#RA_18S <- 
ggplot(data=Relative_abundance_18S, aes(x= Primer_name,y= Rel_abund, fill= phylum)) +
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061", "#c00000","darkolivegreen","#FFDB77FF","lightgrey","#004CFFFF","#969696"))+
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Primer combination ID", y= "Relative abundance (In silico)", tag = "A)")+
  guides(fill= guide_legend(nrow = 8))

#pdf(file = "~/AA_Primer_evaluation/Figures/In_silico_evaluation/Preliminary/Figure_1.pdf", width = 10, height = 8)
#RA_18S
#dev.off()

##28S primer pairs
##Select Primer sequeces and Primer combination ID

Primers %>%
  filter(Gen== "28S")-> Primers_28S

Euk_28S_01<- search_primer_pair(name= as.character(Primers_28S[1,1]), forward = as.character(Primers_28S[1,2]), reverse = as.character(Primers_28S[1,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_28S_01, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_28S_01.Rds")

Euk_28S_02<- search_primer_pair(name= as.character(Primers_28S[2,1]), forward = as.character(Primers_28S[2,2]), reverse = as.character(Primers_28S[2,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_28S_02, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_28S_02.Rds")

Euk_28S_03<- search_primer_pair(name= as.character(Primers_28S[2,1]), forward = as.character(Primers_28S[3,2]), reverse = as.character(Primers_28S[3,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_28S_03, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_28S_03.Rds")

#Make a list of primerTree objects saved
done28S <- list.files(path = "~/AA_Primer_evaluation/output/primerTreeObj/", pattern = "Euk_28S_")
done28S<- lapply(done28S, function(x) gsub(".Rds", "_RA", x))

###Euk_28S_01
plot(Euk_28S_01, ranks= "phylum")
Euk_28S_01_Tax<- Euk_28S_01$taxonomy ##Taxonomy information 
Euk_28S_01_BLAST<- Euk_28S_01$BLAST_result #BLAST information
Euk_28S_01_Tax[is.na(Euk_28S_01_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_28S_01_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_28S_01_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_28S_01_Results
hist(Euk_28S_01_Results$product_length)
Euk_28S_01_RA<- data.frame(count(Euk_28S_01_Results$phylum))
Euk_28S_01_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_28S_01")->Euk_28S_01_RA
colnames(Euk_28S_01_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_28S_02
plot(Euk_28S_02, ranks= "phylum")
Euk_28S_02_Tax<- Euk_28S_02$taxonomy ##Taxonomy information 
Euk_28S_02_BLAST<- Euk_28S_02$BLAST_result #BLAST information
Euk_28S_02_Tax[is.na(Euk_28S_02_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_28S_02_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_28S_02_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_28S_02_Results
hist(Euk_28S_02_Results$product_length)
Euk_28S_02_RA<- data.frame(count(Euk_28S_02_Results$phylum))
Euk_28S_02_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_28S_02")->Euk_28S_02_RA
colnames(Euk_28S_02_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_28S_03
plot(Euk_28S_03, ranks= "phylum")
Euk_28S_03_Tax<- Euk_28S_03$taxonomy ##Taxonomy information 
Euk_28S_03_BLAST<- Euk_28S_03$BLAST_result #BLAST information
Euk_28S_03_Tax[is.na(Euk_28S_03_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_28S_03_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_28S_03_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_28S_03_Results
hist(Euk_28S_03_Results$product_length)
Euk_28S_03_RA<- data.frame(count(Euk_28S_03_Results$phylum))
Euk_28S_03_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_28S_03")->Euk_28S_03_RA
colnames(Euk_28S_03_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

Relative_abundance_28S<- bind_rows(Euk_28S_01_RA, Euk_28S_02_RA, Euk_28S_03_RA) ##Change everytime to add next primer information
  
Relative_abundance_28S%>%
  group_by(Primer_name)%>%
  mutate(Main_taxa= Rel_abund>= 0.05) %>%
  mutate(phylum= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(.$phylum))) %>%
  arrange(Primer_name, desc(phylum))->Relative_abundance_28S
  
#RA_28S <- 
ggplot(data=Relative_abundance_28S, aes(x= Primer_name,y= Rel_abund, fill= phylum)) +
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
                              # "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061", "#c00000",
                              # "aquamarine","darkolivegreen","#FFDB77FF",
                              "lightgrey","#004CFFFF","#969696"))+
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Primer combination ID", y= "Relative abundance (In silico)", tag = "B)")+
  guides(fill= guide_legend(nrow = 8))

#pdf(file = "~/AA_Primer_evaluation/Figures/In_silico_evaluation/Preliminary/Figure_2.pdf", width = 10, height = 8)
#RA_28S
#dev.off()

##Plant primer pairs
##Select Primer sequeces and Primer combination ID

Primers %>%
  filter(Gen%in% c("ITS", "rbcL"))-> Primers_plants

Euk_rbcL_01<- search_primer_pair(name= as.character(Primers_plants[1,1]), forward = as.character(Primers_plants[1,2]), reverse = as.character(Primers_plants[1,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_rbcL_01, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_rbcL_01.Rds")

Euk_ITS_01<- search_primer_pair(name= as.character(Primers_plants[2,1]), forward = as.character(Primers_plants[2,2]), reverse = as.character(Primers_plants[2,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_ITS_01, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_ITS_01.Rds")

###Euk_rbcL_01
plot(Euk_rbcL_01, ranks= "phylum")
Euk_rbcL_01_Tax<- Euk_rbcL_01$taxonomy ##Taxonomy information 
Euk_rbcL_01_BLAST<- Euk_rbcL_01$BLAST_result #BLAST information
Euk_rbcL_01_Tax[is.na(Euk_rbcL_01_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_rbcL_01_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_rbcL_01_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_rbcL_01_Results
hist(Euk_rbcL_01_Results$product_length)
Euk_rbcL_01_RA<- data.frame(count(Euk_rbcL_01_Results$phylum))
Euk_rbcL_01_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_rbcL_01")->Euk_rbcL_01_RA
colnames(Euk_rbcL_01_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_ITS_01
plot(Euk_ITS_01, ranks= "phylum")
Euk_ITS_01_Tax<- Euk_ITS_01$taxonomy ##Taxonomy information 
Euk_ITS_01_BLAST<- Euk_ITS_01$BLAST_result #BLAST information
Euk_ITS_01_Tax[is.na(Euk_ITS_01_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_ITS_01_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_ITS_01_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_ITS_01_Results
hist(Euk_ITS_01_Results$product_length)
Euk_ITS_01_RA<- data.frame(count(Euk_ITS_01_Results$phylum))
Euk_ITS_01_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_ITS_01")->Euk_ITS_01_RA
colnames(Euk_ITS_01_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

Relative_abundance_plants<- bind_rows(Euk_rbcL_01_RA, Euk_ITS_01_RA) ##Change everytime to add next primer information

Relative_abundance_plants%>%
  group_by(Primer_name)%>%
  mutate(Main_taxa= Rel_abund>= 0.05) %>%
  mutate(phylum= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(.$phylum))) %>%
  arrange(Primer_name, desc(phylum))->Relative_abundance_plants

#RA_plants <- 
ggplot(data=Relative_abundance_plants, aes(x= Primer_name,y= Rel_abund, fill= phylum)) +
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
                               # "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061", "#c00000",
                               # "aquamarine","darkolivegreen","#FFDB77FF",
                               "lightgrey","#004CFFFF","#969696"))+
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Primer combination ID", y= "Relative abundance (In silico)", tag = "B)")+
  guides(fill= guide_legend(nrow = 8))

##COI primer pairs
##Select Primer sequeces and Primer combination ID

Primers %>%
  filter(Gen== "COI")-> Primers_COI

Euk_COI_01<- search_primer_pair(name= as.character(Primers_COI[1,1]), forward = as.character(Primers_COI[1,2]), reverse = as.character(Primers_COI[1,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_COI_01, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_COI_01.Rds")

Euk_COI_02<- search_primer_pair(name= as.character(Primers_COI[2,1]), forward = as.character(Primers_COI[2,2]), reverse = as.character(Primers_COI[2,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_COI_02, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_COI_02.Rds")

Euk_COI_03<- search_primer_pair(name= as.character(Primers_COI[3,1]), forward = as.character(Primers_COI[3,2]), reverse = as.character(Primers_COI[3,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_COI_03, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_COI_03.Rds")

Euk_COI_04<- search_primer_pair(name= as.character(Primers_COI[4,1]), forward = as.character(Primers_COI[4,2]), reverse = as.character(Primers_COI[4,3]), .parallel = T, num_aligns = 1000)
#saveRDS(Euk_COI_04, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_COI_04.Rds")

###Euk_COI_01
plot(Euk_COI_01, ranks= "phylum")
Euk_COI_01_Tax<- Euk_COI_01$taxonomy ##Taxonomy information 
Euk_COI_01_BLAST<- Euk_COI_01$BLAST_result #BLAST information
Euk_COI_01_Tax[is.na(Euk_COI_01_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_COI_01_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_COI_01_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_COI_01_Results
hist(Euk_COI_01_Results$product_length)
Euk_COI_01_RA<- data.frame(count(Euk_COI_01_Results$phylum))
Euk_COI_01_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_COI_01")->Euk_COI_01_RA
colnames(Euk_COI_01_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_COI_02
plot(Euk_COI_02, ranks= "phylum")
Euk_COI_02_Tax<- Euk_COI_02$taxonomy ##Taxonomy information 
Euk_COI_02_BLAST<- Euk_COI_02$BLAST_result #BLAST information
Euk_COI_02_Tax[is.na(Euk_COI_02_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_COI_02_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_COI_02_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_COI_02_Results
hist(Euk_COI_02_Results$product_length)
Euk_COI_02_RA<- data.frame(count(Euk_COI_02_Results$phylum))
Euk_COI_02_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_COI_02")->Euk_COI_02_RA
colnames(Euk_COI_02_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_COI_03
plot(Euk_COI_03, ranks= "phylum")
Euk_COI_03_Tax<- Euk_COI_03$taxonomy ##Taxonomy information 
Euk_COI_03_BLAST<- Euk_COI_03$BLAST_result #BLAST information
Euk_COI_03_Tax[is.na(Euk_COI_03_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_COI_03_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_COI_03_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_COI_03_Results
hist(Euk_COI_03_Results$product_length)
Euk_COI_03_RA<- data.frame(count(Euk_COI_03_Results$phylum))
Euk_COI_03_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_COI_03")->Euk_COI_03_RA
colnames(Euk_COI_03_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

###Euk_COI_04
plot(Euk_COI_04, ranks= "phylum")
Euk_COI_04_Tax<- Euk_COI_04$taxonomy ##Taxonomy information 
Euk_COI_04_BLAST<- Euk_COI_04$BLAST_result #BLAST information
Euk_COI_04_Tax[is.na(Euk_COI_04_Tax)]<- "Unassigned" ##Change NA's into Unassigned 
Euk_COI_04_Tax%>%
  select(c("taxId", "gi", "species", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"))%>%
  join(Euk_COI_04_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_COI_04_Results
hist(Euk_COI_04_Results$product_length)
Euk_COI_04_RA<- data.frame(count(Euk_COI_04_Results$phylum))
Euk_COI_04_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_COI_04")->Euk_COI_04_RA
colnames(Euk_COI_04_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

Relative_abundance_COI<- bind_rows(Euk_COI_01_RA, Euk_COI_02_RA, Euk_COI_03_RA, Euk_COI_04_RA) ##Change everytime to add next primer information

Relative_abundance_COI%>%
  group_by(Primer_name)%>%
  mutate(Main_taxa= Rel_abund>= 0.05) %>%
  mutate(phylum= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(.$phylum))) %>%
  arrange(Primer_name, desc(phylum))->Relative_abundance_COI

#RA_COI <- 
ggplot(data=Relative_abundance_COI, aes(x= Primer_name,y= Rel_abund, fill= phylum)) +
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C",#"#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
                               #"#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061", "#c00000",
                               #"aquamarine","darkolivegreen","#FFDB77FF",
                               "lightgrey","#004CFFFF","#969696"))+
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Primer combination ID", y= "Relative abundance (In silico)", tag = "B)")+
  guides(fill= guide_legend(nrow = 8))


##Merge all
Relative_abundance<- bind_rows(Relative_abundance_18S, Relative_abundance_28S, 
                               Relative_abundance_COI, Relative_abundance_plants)

#RA_Euk<-
ggplot(data=Relative_abundance, aes(x= Primer_name,y= Rel_abund, fill= phylum)) +
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C", "#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061", "#c00000",
                               "aquamarine","darkolivegreen", "burlywood1","#FFDB77FF",
                               "lightgrey","#004CFFFF","#969696"
                               ))+
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Primer combination ID", y= "Relative abundance (In silico)", tag = "B)")+
  guides(fill= guide_legend(nrow = 8))

#pdf(file = "~/AA_Primer_evaluation/Figures/In_silico_evaluation/Preliminary/Figure_3.pdf", width = 10, height = 8)
#RA_Euk
#dev.off()
