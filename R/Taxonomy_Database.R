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
n_distinct(Taxonomy_18S$superkingdom, na.rm = T)
n_distinct(Taxonomy_18S$phylum, na.rm = T)
n_distinct(Taxonomy_18S$class, na.rm = T)
n_distinct(Taxonomy_18S$order, na.rm = T)
n_distinct(Taxonomy_18S$family, na.rm = T)
n_distinct(Taxonomy_18S$genus, na.rm = T)
n_distinct(Taxonomy_18S$species, na.rm = T)

