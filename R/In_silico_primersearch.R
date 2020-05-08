##Analysis In silico PCR results (primerSearch EMBOSS tool) 
##Load packages
library("Biostrings")
library("DECIPHER")
library("annotate")
library("taxonomizr")
library("metacoder")

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

##Load raw results 
##In silico PCR against ENA full 18S database
raw_output <- readLines("/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/output/primersearch_18S_ENA")


Primer_search_18S_results<- parse_primersearch(raw_output = raw_output)
