##Analysis In silico PCR results (primerSearch EMBOSS tool) 
##Load packages
library("Biostrings")
library("DECIPHER")
library("annotate")
library("taxonomizr")
library("metacoder")
library("plyr")
library("dplyr")
library("tidyverse")

##Functions
#Modification from metabarcoder::parse_primersearch to reshape primersearch output into data frame containing:
#Primer name (not included in the original function)
#Sequence name (input) --> for our db is directly the taxid  
#Product amplicon (not included in the original function) 

parse_primersearch <- function(raw_output) {
  # Split output into chunks for each primer
  primer_indexes <- grep("Primer name ", raw_output, fixed = TRUE, value = FALSE)
  primer_chunk_id <- findInterval(seq_along(raw_output), primer_indexes)
  primer_chunks <- vapply(split(raw_output, primer_chunk_id)[-1],
                          paste, character(1), collapse = "\n")
  names(primer_chunks) <- stringr::str_match(primer_chunks, "Primer name ([^\n]*)")[,2]
  # Extract amplicon data from each chunk and combine
  pattern <- paste("Amplimer ([0-9]+)",
                   "\tSequence: ([^\n]*)",
                   "\t([^\n]*)",
                   "\t([^\n]+) hits forward strand at ([0-9]+) with ([0-9]+) mismatches",
                   "\t([^\n]+) hits reverse strand at \\[([0-9]+)\\] with ([0-9]+) mismatches",
                   "\tAmplimer length: ([0-9]+) bp", sep = '\n')
  primer_data <- stringr::str_match_all(primer_chunks, pattern)
  primer_data <- as.data.frame(cbind(rep(names(primer_chunks), vapply(primer_data, nrow, numeric(1))),
                                     do.call(rbind, primer_data)[, -1]), stringsAsFactors = FALSE)
  # Reformat amplicon data
  colnames(primer_data) <- c("pair_name", "amplimer", "input", "name", "f_primer", "f_start",
                             "f_mismatch",  "r_primer", "r_start", "r_mismatch", "product_length")
  primer_data <- primer_data[, c("pair_name","input",  "f_primer", "f_start", "f_mismatch",
                                 "r_primer", "r_start", "r_mismatch", "product_length")]
  numeric_cols <- c("f_start", "f_mismatch", "r_start", "r_mismatch", "input", "product_length")
  for (col in numeric_cols) primer_data[, col] <- as.numeric(primer_data[, col]) 
  return(primer_data)
} 

###Get taxonomy: retrieve taxonomic information based on taxID values from in silico PCR
#input: reformated raw primersearch results for a given primer combination
taxResults <- function(input){
  input_Tax<- taxonomizr::getTaxonomy(as.character(input$taxID), sqlFile = "/SAN/db/taxonomy/taxonomizr.sql")
  input_Tax<- as.data.frame(input_Tax)
  rownames(input_Tax) <- c(1:nrow(input_Tax))
  input_Tax$taxID<- input$taxID
  return(input_Tax)
}

###Compile results for the rest of analysis: Extract unique taxa within a range of product size 
#taxdata: result from taxResults()
#insilico: reformated raw primersearch results for a given primer combination

compInsilico<- function(taxdata, insilico){
  taxdata%>%
    join(insilico, by= "taxID")%>%
    dplyr::filter(between(product_length, 200, 700))%>% ##Defined by library preparation
    distinct(species, .keep_all= T)-> insilico_Results ##Keep a single ocurrence of each species 
  hist(insilico_Results$product_length)
  return(insilico_Results)
}

##Convert objects with in silico taxonomy results for 18S primer combinations from global environment into a list
env2list18S <- function(env){
  names = ls(env, pattern = "18S_\\d+_Results$")
  mget(names, env)
}

##Activate genetic marker 
Primers18S<- T
Primers28S<- F
Primers16S<- F
PrimersCOI<- F
PrimersOther<- F ##ITS, rbcL, 12S 

##Load raw results 
if(Primers18S){
##In silico PCR against ENA full 18S database
raw_output <- readLines("/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/output/primersearch_18S_ENA")
}

##Analysis
if(Primers18S){
##Transform the raw outpunt into a data frame with all the results
Primer_search_18S_results<- parse_primersearch(raw_output = raw_output)
Primer_search_18S_results$taxID<- stringr::str_extract(Primer_search_18S_results$input, "\\d+") ##create taxID colum for retrieve taxonomic information (seq name contains .#)

##Create a list to split information into each primer combination
primers18S_list<- split(Primer_search_18S_results, Primer_search_18S_results$pair_name)
#saveRDS(primers18S_list, "/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/output/primersearch_18S_list.rds")

##Send each to the global environment
list2env(primers18S_list, envir = .GlobalEnv)

##Work with the results from individual primer combinations
##Euk_18S_01
##convert taxonomic IDs to taxonomy
Euk_18S_01_Tax<- taxResults(input = Euk_18S_01)
##Summarize results for later compilation 
Euk_18S_01_Results<- compInsilico(taxdata = Euk_18S_01_Tax, insilico = Euk_18S_01)
Euk_18S_01_RA<- data.frame(dplyr::count(Euk_18S_01_Results, phylum))
Euk_18S_01_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_01")-> Euk_18S_01_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_01_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_01, Euk_18S_01_Tax)

##Euk_18S_02
##convert taxonomic IDs to taxonomy
Euk_18S_02_Tax<- taxResults(input = Euk_18S_02)
##Summarize results for later compilation 
Euk_18S_02_Results<- compInsilico(taxdata = Euk_18S_02_Tax, insilico = Euk_18S_02)
Euk_18S_02_RA<- data.frame(dplyr::count(Euk_18S_02_Results, phylum))
Euk_18S_02_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_02")-> Euk_18S_02_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_02_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_02, Euk_18S_02_Tax)

##Euk_18S_03
##convert taxonomic IDs to taxonomy
Euk_18S_03_Tax<- Euk_18S_03_Tax<- taxResults(input = Euk_18S_03)
##Summarize results for later compilation 
Euk_18S_03_Results<- compInsilico(taxdata = Euk_18S_03_Tax, insilico = Euk_18S_03)
Euk_18S_03_RA<- data.frame(dplyr::count(Euk_18S_03_Results, phylum))
Euk_18S_03_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_03")-> Euk_18S_03_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_03_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_03, Euk_18S_03_Tax)

##Euk_18S_04
Euk_18S_04_Tax<- taxResults(input = Euk_18S_04)
##Summarize results for later compilation 
Euk_18S_04_Results<- compInsilico(taxdata = Euk_18S_04_Tax, insilico = Euk_18S_04)
Euk_18S_04_RA<- data.frame(dplyr::count(Euk_18S_04_Results, phylum))
Euk_18S_04_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_04")-> Euk_18S_04_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_04_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_04, Euk_18S_04_Tax)

##Euk_18S_05
Euk_18S_05_Tax<- taxResults(input = Euk_18S_05)
##Summarize results for later compilation 
Euk_18S_05_Results<- compInsilico(taxdata = Euk_18S_05_Tax, insilico = Euk_18S_05)
Euk_18S_05_RA<- data.frame(dplyr::count(Euk_18S_05_Results, phylum))
Euk_18S_05_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_05")-> Euk_18S_05_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_05_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_05, Euk_18S_05_Tax)

##Euk_18S_06
Euk_18S_06_Tax<- taxResults(input = Euk_18S_06)
##Summarize results for later compilation 
Euk_18S_06_Results<- compInsilico(taxdata = Euk_18S_06_Tax, insilico = Euk_18S_06)
Euk_18S_06_RA<- data.frame(dplyr::count(Euk_18S_06_Results, phylum))
Euk_18S_06_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_06")-> Euk_18S_06_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_06_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_06, Euk_18S_06_Tax)

##Euk_18S_07
Euk_18S_07_Tax<- taxResults(input = Euk_18S_07)
##Summarize results for later compilation 
Euk_18S_07_Results<- compInsilico(taxdata = Euk_18S_07_Tax, insilico = Euk_18S_07)
Euk_18S_07_RA<- data.frame(dplyr::count(Euk_18S_07_Results, phylum))
Euk_18S_07_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_07")-> Euk_18S_07_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_07_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_07, Euk_18S_07_Tax)

##Euk_18S_08
Euk_18S_08_Tax<- taxResults(input = Euk_18S_08)
##Summarize results for later compilation 
Euk_18S_08_Results<- compInsilico(taxdata = Euk_18S_08_Tax, insilico = Euk_18S_08)
Euk_18S_08_RA<- data.frame(dplyr::count(Euk_18S_08_Results, phylum))
Euk_18S_08_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_08")-> Euk_18S_08_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_08_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_08, Euk_18S_08_Tax)

##Euk_18S_09
Euk_18S_09_Tax<- taxResults(input = Euk_18S_09)
##Summarize results for later compilation 
Euk_18S_09_Results<- compInsilico(taxdata = Euk_18S_09_Tax, insilico = Euk_18S_09)
Euk_18S_09_RA<- data.frame(dplyr::count(Euk_18S_09_Results, phylum))
Euk_18S_09_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_09")-> Euk_18S_09_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_09_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_09, Euk_18S_09_Tax)

##Euk_18S_10
Euk_18S_10_Tax<- taxResults(input = Euk_18S_10)
##Summarize results for later compilation 
Euk_18S_10_Results<- compInsilico(taxdata = Euk_18S_10_Tax, insilico = Euk_18S_10)
Euk_18S_10_RA<- data.frame(dplyr::count(Euk_18S_10_Results, phylum))
Euk_18S_10_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_10")-> Euk_18S_10_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_10_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_10, Euk_18S_10_Tax)

##Euk_18S_11
Euk_18S_11_Tax<- taxResults(input = Euk_18S_11)
##Summarize results for later compilation 
Euk_18S_11_Results<- compInsilico(taxdata = Euk_18S_11_Tax, insilico = Euk_18S_11)
Euk_18S_11_RA<- data.frame(dplyr::count(Euk_18S_11_Results, phylum))
Euk_18S_11_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_11")-> Euk_18S_11_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_11_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_11, Euk_18S_11_Tax)

##Euk_18S_12
Euk_18S_12_Tax<- taxResults(input = Euk_18S_12)
##Summarize results for later compilation 
Euk_18S_12_Results<- compInsilico(taxdata = Euk_18S_12_Tax, insilico = Euk_18S_12)
Euk_18S_12_RA<- data.frame(dplyr::count(Euk_18S_12_Results, phylum))
Euk_18S_12_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_12")-> Euk_18S_12_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_12_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_12, Euk_18S_12_Tax)

##Euk_18S_13
Euk_18S_13_Tax<- taxResults(input = Euk_18S_13)
##Summarize results for later compilation 
Euk_18S_13_Results<- compInsilico(taxdata = Euk_18S_13_Tax, insilico = Euk_18S_13)
Euk_18S_13_RA<- data.frame(dplyr::count(Euk_18S_13_Results, phylum))
Euk_18S_13_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_13")-> Euk_18S_13_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_13_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_13, Euk_18S_13_Tax)

##Euk_18S_14
Euk_18S_14_Tax<- taxResults(input = Euk_18S_14)
##Summarize results for later compilation 
Euk_18S_14_Results<- compInsilico(taxdata = Euk_18S_14_Tax, insilico = Euk_18S_14)
Euk_18S_14_RA<- data.frame(dplyr::count(Euk_18S_14_Results, phylum))
Euk_18S_14_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_14")-> Euk_18S_14_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_14_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_14, Euk_18S_14_Tax)

##Euk_18S_15
Euk_18S_15_Tax<- taxResults(input = Euk_18S_15)
##Summarize results for later compilation 
Euk_18S_15_Results<- compInsilico(taxdata = Euk_18S_15_Tax, insilico = Euk_18S_15)
Euk_18S_15_RA<- data.frame(dplyr::count(Euk_18S_15_Results, phylum))
Euk_18S_15_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_15")-> Euk_18S_15_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_15_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_15, Euk_18S_15_Tax)

##Euk_18S_16
Euk_18S_16_Tax<- taxResults(input = Euk_18S_16)
##Summarize results for later compilation 
Euk_18S_16_Results<- compInsilico(taxdata = Euk_18S_16_Tax, insilico = Euk_18S_16)
Euk_18S_16_RA<- data.frame(dplyr::count(Euk_18S_16_Results, phylum))
Euk_18S_16_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_16")-> Euk_18S_16_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_16_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_16, Euk_18S_16_Tax)

##Euk_18S_17
Euk_18S_17_Tax<- taxResults(input = Euk_18S_17)
##Summarize results for later compilation 
Euk_18S_17_Results<- compInsilico(taxdata = Euk_18S_17_Tax, insilico = Euk_18S_17)
Euk_18S_17_RA<- data.frame(dplyr::count(Euk_18S_17_Results, phylum))
Euk_18S_17_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_17")-> Euk_18S_17_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_17_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_17, Euk_18S_17_Tax)

##Euk_18S_18
Euk_18S_18_Tax<- taxResults(input = Euk_18S_18)
##Summarize results for later compilation 
Euk_18S_18_Results<- compInsilico(taxdata = Euk_18S_18_Tax, insilico = Euk_18S_18)
Euk_18S_18_RA<- data.frame(dplyr::count(Euk_18S_18_Results, phylum))
Euk_18S_18_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_18")-> Euk_18S_18_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_18_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_18, Euk_18S_18_Tax)

##Euk_18S_19
Euk_18S_19_Tax<- taxResults(input = Euk_18S_19)
##Summarize results for later compilation 
Euk_18S_19_Results<- compInsilico(taxdata = Euk_18S_19_Tax, insilico = Euk_18S_19)
Euk_18S_19_RA<- data.frame(dplyr::count(Euk_18S_19_Results, phylum))
Euk_18S_19_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_19")-> Euk_18S_19_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_19_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_19, Euk_18S_19_Tax)

##Euk_18S_20
Euk_18S_20_Tax<- taxResults(input = Euk_18S_20)
##Summarize results for later compilation 
Euk_18S_20_Results<- compInsilico(taxdata = Euk_18S_20_Tax, insilico = Euk_18S_20)
Euk_18S_20_RA<- data.frame(dplyr::count(Euk_18S_20_Results, phylum))
Euk_18S_20_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_20")-> Euk_18S_20_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_20_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_20, Euk_18S_20_Tax)

##Euk_18S_21
Euk_18S_21_Tax<- taxResults(input = Euk_18S_21)
##Summarize results for later compilation 
Euk_18S_21_Results<- compInsilico(taxdata = Euk_18S_21_Tax, insilico = Euk_18S_21)
Euk_18S_21_RA<- data.frame(dplyr::count(Euk_18S_21_Results, phylum))
Euk_18S_21_RA%>%
  mutate(Rel_abund= n/sum(n))%>%
  mutate(Primer_name= "Euk_18S_21")-> Euk_18S_21_RA ##Data frame with all the relative abundance in silico per primer pair
colnames(Euk_18S_21_RA)<- c("phylum", "Freq", "Rel_abund", "Primer_name")

rm(Euk_18S_21, Euk_18S_21_Tax)

##Relative abundance in silico comparison
##Phylum level
Relative_abundance_18S<- bind_rows(Euk_18S_01_RA, Euk_18S_02_RA, Euk_18S_03_RA, Euk_18S_04_RA, Euk_18S_05_RA, Euk_18S_06_RA,
                                   Euk_18S_07_RA, Euk_18S_08_RA, Euk_18S_09_RA, Euk_18S_10_RA, Euk_18S_11_RA, Euk_18S_12_RA,
                                   Euk_18S_13_RA, Euk_18S_14_RA, Euk_18S_15_RA, Euk_18S_16_RA, Euk_18S_17_RA, Euk_18S_18_RA,
                                   Euk_18S_19_RA, Euk_18S_20_RA, Euk_18S_21_RA)

Relative_abundance_18S%>%
  group_by(Primer_name)%>%
  mutate(Main_taxa= Rel_abund>= 0.05) %>%
  dplyr::mutate(phylum= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(phylum)))%>%
  arrange(Primer_name, desc(phylum))->Relative_abundance_18S

rm(Euk_18S_01_RA, Euk_18S_02_RA, Euk_18S_03_RA, Euk_18S_04_RA, Euk_18S_05_RA, Euk_18S_06_RA,
   Euk_18S_07_RA, Euk_18S_08_RA, Euk_18S_09_RA, Euk_18S_10_RA, Euk_18S_11_RA, Euk_18S_12_RA,
   Euk_18S_13_RA, Euk_18S_14_RA, Euk_18S_15_RA, Euk_18S_16_RA, Euk_18S_17_RA, Euk_18S_18_RA,
   Euk_18S_19_RA, Euk_18S_20_RA, Euk_18S_21_RA)

#RA_18S <- 
ggplot(data=Relative_abundance_18S, aes(x= Primer_name,y= Rel_abund, fill= phylum)) +
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061", "#c00000","darkolivegreen","#FFDB77FF","lightgrey","#004CFFFF","#969696"))+
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Primer combination ID", y= "Relative abundance (In silico)", tag = "A)")+
  guides(fill= guide_legend(nrow = 8))

#pdf(file = "~/AA_Primer_evaluation/Figures/In_silico_evaluation/Figure_1_primersearch18S.pdf", width = 10, height = 8)
#RA_18S
#dev.off()

}
