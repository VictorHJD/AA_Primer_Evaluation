##Get all database content information
##Load Libraries
library("doMC")
library("ggplot2")
library("reshape2")
library("plyr")
library("ape")
library("TmCalculator")
library("dplyr")
library("Biostrings")
library("DECIPHER")
library("annotate")
library("taxonomizr")

##
db18S<- T
db28S<- F
dbCOI<- F
db16S<- F
dbITS<- F
dbrbcL<- F

##
if(db18S){
##Load data base 18S 
#Seq_18S_db<- Biostrings::readDNAStringSet("/SAN/db/ENA_marker/18S_All_ENA.fasta.gz1700cleaned.fasta", format = "fasta") ##Seq names un formated
Seq_18S_db<- Biostrings::readDNAStringSet("/SAN/db/blastdb/18S_ENA/Full1700_18S.fasta", format = "fasta") ##Seq names are already taxID

##Align seqs
#AlignSeqs(Seq_18S_db)

##Extract names from the sequences 
Seq_18S_names<- Seq_18S_db@ranges@NAMES
#accession_18S <- as.data.frame(str_extract(Seq_18S_names, "[:alpha:][:alpha:]\\d+.1")) ##Extract accession numbers from them 
#colnames(accession_18S)<- "Accesion_number"

##convert accession numbers to taxonomic IDs 
##SQLite database is already in Harriet
#taxID_18S<- accessionToTaxa(as.character(accession_18S$Accesion_number), sqlFile = "/SAN/db/taxonomy/taxonomizr.sql") ##Not working :S
#accession_18S$Taxonomic_ID<- taxID_18S
taxID_18S<- as.data.frame(str_extract(Seq_18S_names, "\\d+"))
colnames(taxID_18S)<- "Taxonomic_ID"
#accession_18S<- taxonomizr::getAccessions(as.character(taxID_18S$Taxonomic_ID), sqlFile = "/SAN/db/taxonomy/taxonomizr.sql")
#accession_18S$Taxonomic_ID<- taxID_18S
#rm(taxID_18S)

##Create a file for mapping taxid 
ENA_18S_db<-as.data.frame(Seq_18S_names)
ENA_18S_db<- bind_cols(ENA_18S_db, taxID_18S)
#write.table(ENA_18S_db, file = "/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/input/ENA_18S_db.txt", sep = "\t",
#            row.names = FALSE, col.names = FALSE, quote = FALSE)

##convert taxonomic IDs to taxonomy
Taxonomy_18S<- taxonomizr::getTaxonomy(as.character(taxID_18S$Taxonomic_ID), sqlFile = "/SAN/db/taxonomy/taxonomizr.sql")

Taxonomy_18S<- as.data.frame(Taxonomy_18S)
rownames(Taxonomy_18S) <- c(1:nrow(Taxonomy_18S))
Taxonomy_18S$Taxonomic_ID<- taxID_18S$Taxonomic_ID
Taxonomy_18S$Accession_18S<- NA

##Get counts of unique taxa per taxonomic level
Unique_taxa<- as.data.frame(sapply(Taxonomy_18S, function(x) n_distinct(x, na.rm = T)))
colnames(Unique_taxa) <- "Database"
}

rm(Taxonomy_18S, ENA_18S_db, taxID_18S, Seq_18S_names, seq_18S_db)

