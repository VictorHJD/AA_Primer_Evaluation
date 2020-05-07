##Analysis In silico PCR results (primerSearch EMBOSS tool) 
##Load packages
library("Biostrings")
library("DECIPHER")
library("annotate")
library("taxonomizr")
library("metacoder")

#Using as reference the function from metacoder package transform primersearch output into an R readable file  

raw_output <- readLines("/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/output/primersearch_18S_ENA")

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
                             "f_mismatch",  "r_primer", "r_start", "r_mismatch")
  primer_data <- primer_data[, c("input",  "f_primer", "f_start", "f_mismatch",
                                 "r_primer", "r_start", "r_mismatch")]
  numeric_cols <- c("f_start", "f_mismatch", "r_start", "r_mismatch", "input")
  for (col in numeric_cols) primer_data[, col] <- as.numeric(primer_data[, col]) 
  return(primer_data)
} 