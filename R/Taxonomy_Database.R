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

##Functions
##Convert objects with in silico taxnomoy for 18S results in the global environment into a list
env2list <- function(env){
  names = ls(env, pattern = "18S_\\d+_Results$")
  mget(names, env)
}

##
source("~/GitProjects/AA_Primer_Evaluation/R/In_silico_primer_evaluation.R")

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
insilico_18S <- env2list(.GlobalEnv)

names_18S <- gsub("_Results", "\\1", names(insilico_18S))

##Data frame with total unique taxa by primer combinations
final<- list() ###Start with an empty list 

for (i in 1: length(insilico_18S)) ### Start a loop: for every element in the list ...
{ 
  tmp <- data.frame() #### make an individual data frame ...
  
  {
    tmp <- as.data.frame(sapply(insilico_18S[[i]], function(x) n_distinct(x, na.rm = T)))  ###Make a data frame with the data included in each element of the list 
    colnames(tmp)<- i
    tmp[,2]<-rownames(tmp) ### And use the rownames as information of the second column
    tmp<- slice(tmp, match(c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "taxId", "accession"), V2))
    tmp<- dplyr::mutate(tmp, V3= c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S")) ##Vector with desired order
    tmp<-  column_to_rownames(tmp, "V3") 
    tmp<- dplyr::select(tmp, -V2)
    
  }
  final[[i]] <- tmp ### Join all the "individual" data frames into the list 
}

final<- data.frame(final)

colnames(final)<- names_18S

Unique_taxa<- cbind(Unique_taxa, final)

rm(tmp)

write.csv(Unique_taxa, "~/AA_Primer_evaluation/Unique_taxa_by_primer_combination.csv")

final_trans <- data.table::transpose(final)
rownames(final_trans)<- colnames(final)
colnames(final_trans)<- rownames(final)

##Data frame with cummulative new taxa by primer combination
names(insilico_18S)<- names_18S

alltaxa<- data.frame() ###Start with an empty data frame

for (i in 1: length(insilico_18S)) ### Start a loop: for every element in the list ...
{ 
  tmp <- data.frame() #### make an individual data frame ...
  {
    tmp <- as.data.frame(distinct(insilico_18S[[i]]))  ###Make a data frame with the data included in each element of the list 
    tmp<- dplyr::select(tmp, c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "taxId", "accession"))
    colnames(tmp)<- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "Taxonomic_ID", "Accession_18S")
    tmp<- dplyr::mutate(tmp, Primer_comb_ID= names(insilico_18S[i])) ##Add a colum with Primer combination name
  }
  alltaxa<- rbind(alltaxa, tmp) ### Join all the "individual" data frames into the final data frame 
}

sapply(alltaxa, function(x) n_distinct(x, na.rm = T))

unique(alltaxa$species)
