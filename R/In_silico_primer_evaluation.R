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

##PrimerTree pipeline
##Considering just the first 1000 alignments
##For degenerated primers uo to 25 combinations are tested for now
##Try with one primer pair

Euk_18S_01<- search_primer_pair(name='Euk_18S_01', 'CAGCAGCCGCGGTAATTCC', 'TCCGTCAATTTCTTTAAGTTTCAG', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_01, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_01.Rds")

Euk_18S_02<- search_primer_pair(name='Euk_18S_02', 'GCTTGTCTCAAAGATTAAGCCN', 'CAGGCYCSCTCTCCGGA', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_02, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_02.Rds")

Euk_18S_03<- search_primer_pair(name='Euk_18S_03', 'GCCGCGGTAATTCCAGCT', 'CCCGTGTTGAGTCAAATTAAGC', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_03, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_03.Rds")

Euk_18S_04<- search_primer_pair(name='Euk_18S_04', 'GCCGCGGTAATTCCAGCT', 'TCCGTCAATTTCTTTAAGTTTCAG', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_04, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_04.Rds")

Euk_18S_05<- search_primer_pair(name='Euk_18S_05', 'GCTTGTCTCAAAGATTAAGCCN', 'TCTCAKGCKCCYTCTCCG', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_05, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_05.Rds")

Euk_18S_06<- search_primer_pair(name='Euk_18S_06', 'TGGTGCCAGCAGCCGC', 'TCCGTCAATTTCTTTAAGTTTCAG', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_06, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_06.Rds")

Euk_18S_07<- search_primer_pair(name='Euk_18S_07', 'CCGCGGTAATTCCAGCTC', 'TKGTGGTRCCMYTCCGTCAA', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_07, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_07.Rds")

Euk_18S_08<- search_primer_pair(name='Euk_18S_08', 'CCGCGGTAATTCCAGCTC', 'GGTGCCCTTCCGTCAATTC', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_08, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_08.Rds")

Euk_18S_09<- search_primer_pair(name='Euk_18S_09', 'ATAACAGGTCTGTGATGCCCTT', 'TGATCCTTCTGCAGGTTCACC', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_09, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_09.Rds")

Euk_18S_10<- search_primer_pair(name='Euk_18S_10', 'GCCGCGGTAATTCCAGCT', 'ARACATTCTTGGCAAATGCYTTC', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_10, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_10.Rds")

Euk_18S_11<- search_primer_pair(name='Euk_18S_11', 'TGCCAGTAGTCATATGCTTGTYT', 'CAGGCYCSCTCTCCGGA', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_11, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_11.Rds")

Euk_18S_12<- search_primer_pair(name='Euk_18S_12', 'GTCTGGTGCCAGCASCC', 'GGTGCCCTTCCGTCAATTC', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_12, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_12.Rds")

Euk_18S_13<- search_primer_pair(name='Euk_18S_13', 'YGGTGRTGCATGGCCGYT', 'GGGCGGTGTGTACAAAGG', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_13, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_13.Rds")

Euk_18S_14<- search_primer_pair(name='Euk_18S_14', 'ATCCTGCCAGTAGTCATATGCTTG', 'TACGCTWYTGGAGCTGGAGTTACCG', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_14, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_14.Rds")

Euk_18S_15<- search_primer_pair(name='Euk_18S_15', 'CAGCAGCCGCGGTAATTCC', 'ARACATTCTTGGCAAATGCYTTC', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_15, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_15.Rds")

Euk_18S_16<- search_primer_pair(name='Euk_18S_16', 'YGGTGRTGCATGGCCGYT', 'TGATCCTTCTGCAGGTTCACC', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_16, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_16.Rds")

Euk_18S_17<- search_primer_pair(name='Euk_18S_17', 'TTGACGGARKGGYACCACMA', 'GGGCGGTGTGTACAAAGG', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_17, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_17.Rds")

Euk_18S_18<- search_primer_pair(name='Euk_18S_18', 'GCCATGCATGYCTAAGTATMA', 'AGGCTCCYTCTCCGG', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_18, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_18.Rds")

Euk_18S_19<- search_primer_pair(name='Euk_18S_19', 'YGGTGRTGCATGGCCGYT', 'TCCTTCTGCAGGTTCACCTAC', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_19, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_19.Rds")

Euk_18S_20<- search_primer_pair(name='Euk_18S_20', 'AAGCCATGCATGYCTAAGTATM', 'TCTCAGGCTCCYTCTCC', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_20, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_20.Rds")

Euk_18S_21<- search_primer_pair(name='Euk_18S_21', 'GAATTGACGGAAGGGCACC', 'AAGGGCATCACAGACCTGTTAT', .parallel = T, num_aligns = 1000)
#saveRDS(Euk_18S_21, file = "~/AA_Primer_evaluation/output/primerTreeObj/Euk_18S_21.Rds")

##General overview of the results
plot(Euk_18S_01, ranks= "phylum")

##Extract the relevant information
Euk_18S_01_Tax<- Euk_18S_01$taxonomy ##Taxonomy information 
Euk_18S_01_BLAST<- Euk_18S_01$BLAST_result #BLAST information

Euk_18S_01_Tax[is.na(Euk_18S_01_Tax)]<- "Unassigned" ##Change NA's into Unassigned 

Euk_18S_01_Tax%>%
  select(1:5,7,9,11,12,13)%>%
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
  select(1:5,7,9,11,12,13)%>%
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
  select(1:10)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:5,7,9,11:13)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:9,11)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:5,7,9,11:13)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:6,8,12,16,17)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:5,7,9:12)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:6,8:10,12)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:6,9,11:13)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:6,8:11)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:9,11)%>% ##Check for each primer pair because not always the same taxa are retrived :S
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
  select(1:9,11)%>% ##Check for each primer pair because not always the same taxa are retrived :S
  join(Euk_18S_13_BLAST, by= "gi")%>%
  distinct(species, .keep_all= T)%>%
  select(-(X1))-> Euk_18S_13_Results
hist(Euk_18S_13_Results$product_length)
Euk_18S_13_RA<- data.frame(count(Euk_18S_13_Results$phylum))
Euk_18S_13_RA%>%
  mutate(Rel_abund= freq/sum(freq))%>%
  mutate(Primer_name= "Euk_18S_13")->Euk_18S_13_RA

colnames(Euk_18S_13_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")
##Relative abundance in silico comparison

#Relative_abundance<- bind_rows(Euk_18S_01_RA, Euk_18S_02_RA)

Relative_abundance<- bind_rows(Relative_abundance, Euk_18S_13_RA) ##Change everytime to add next primer information

Relative_abundance%>%
  group_by(Primer_name)%>%
  mutate(Main_taxa= Rel_abund>= 0.05) %>%
  mutate(phylum= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(.$phylum))) %>%
  arrange(Primer_name, desc(phylum))->Relative_abundance

ggplot(data=Relative_abundance, aes(x= Primer_name,y= Rel_abund, fill= phylum)) +
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061", "lightpink","pink","lightgrey","#969696"))+
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Primer combination ID", y= "Relative abundance", tag = "A)")+
  guides(fill= guide_legend(nrow = 8))

