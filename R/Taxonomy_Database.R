##Get all database content information
##Load Libraries
library("doMC")
library("ggplot2")
library("reshape2")
library("plyr")
library("ape")
library("TmCalculator")
library("dplyr")
library("PrimerMiner")
library("Biostrings")
library("DECIPHER")
library("annotate")
library("taxonomizr")

##Load data base 18S 
Seq_18S_db<- readDNAStringSet("/SAN/db/blastdb/18S_ENA/18S_All_ENA.fasta.gzcleaned.fasta", format = "fasta")
##Align seqs
#AlignSeqs(Seq_18S_db)

##Extract names from the sequences 
Seq_18S_names<- Seq_18S_db@ranges@NAMES
accession_18S <- as.data.frame(str_extract(Seq_18S_names, "KY\\d+.1")) ##Extract accession numbers from them 
colnames(accession_18S)<- "Accesion_number"

##convert accession numbers to taxonomic IDs 
##SQLite database is already in Harriet
taxID_18S<- accessionToTaxa(as.character(accession_18S$Accesion_number), sqlFile = "/SAN/db/taxonomy/taxonomizr.sql")
accession_18S$Taxonomic_ID<- taxID_18S
rm(taxID_18S)

##convert taxonomic IDs to taxonomy
taxonomy_18S_DB<- getTaxonomy(as.character(accession_18S$Taxonomic_ID), sqlFile = "/SAN/db/taxonomy/taxonomizr.sql")

temporal<- as.data.frame(taxonomy_18S_DB)
rm(taxonomy_18S_DB)
rownames(temporal) <- c(1:nrow(temporal))
temporal$Taxonomic_ID<- accession_18S$Taxonomic_ID
temporal$Accession_18S<- accession_18S$Accesion_number

Taxonomy_18S<- temporal
rm(temporal)

##Get counts of unique taxa per taxonomic level
Unique_taxa<- as.data.frame(sapply(Taxonomy_18S, function(x) n_distinct(x, na.rm = T)))
colnames(Unique_taxa) <- "Database"

##Get counts of unique taxa per primer pair using in silico predictions on PrimerTree

#test<- list.files(path = "~/AA_Primer_evaluation/output/primerTreeObj", pattern = ".Rds")
#insilico <- ls(envir = .GlobalEnv, pattern = "_Results$")

y<- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "taxId", "accession") ##Vector with desired order

##Euk_18S_01
tmp<- as.data.frame(sapply(Euk_18S_01_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_01"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_02
tmp<- as.data.frame(sapply(Euk_18S_02_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_02"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_03
tmp<- as.data.frame(sapply(Euk_18S_03_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_03"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_04
tmp<- as.data.frame(sapply(Euk_18S_04_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_04"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_05
tmp<- as.data.frame(sapply(Euk_18S_05_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_05"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_06
tmp<- as.data.frame(sapply(Euk_18S_06_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_06"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_07
tmp<- as.data.frame(sapply(Euk_18S_07_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_07"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_08
tmp<- as.data.frame(sapply(Euk_18S_08_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_08"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_09
tmp<- as.data.frame(sapply(Euk_18S_09_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_09"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_10
tmp<- as.data.frame(sapply(Euk_18S_10_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_10"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_11
tmp<- as.data.frame(sapply(Euk_18S_11_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_11"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_12
tmp<- as.data.frame(sapply(Euk_18S_12_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_12"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_13
tmp<- as.data.frame(sapply(Euk_18S_13_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_13"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_14
tmp<- as.data.frame(sapply(Euk_18S_14_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_14"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_15
tmp<- as.data.frame(sapply(Euk_18S_15_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_15"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_16
tmp<- as.data.frame(sapply(Euk_18S_16_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_16"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_17
tmp<- as.data.frame(sapply(Euk_18S_17_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_17"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_18
tmp<- as.data.frame(sapply(Euk_18S_18_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_18"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_19
tmp<- as.data.frame(sapply(Euk_18S_19_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_19"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_20
tmp<- as.data.frame(sapply(Euk_18S_20_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_20"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)

##Euk_18S_21
tmp<- as.data.frame(sapply(Euk_18S_21_Results, function(x) n_distinct(x, na.rm = T)))
colnames(tmp) <- "Euk_18S_21"
tmp[,2]<-rownames(tmp)
tmp%>%
  slice(match(y, V2))%>%
  dplyr::select(-V2)%>%
  mutate(V2= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S"))%>%
  column_to_rownames("V2")%>%
  cbind(Unique_taxa)-> Unique_taxa
rm(tmp)


