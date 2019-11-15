###Information from all the primer pairs evaluated 
library("TmCalculator")

####Primer information
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

rm(GCF, GCR, Fwd, Rev, row, TmF, TmR)
#write.csv(primerInput, "~/AA_Primer_evaluation/primerOutputUnique.csv")

#n16S <- subset(primerInput, primerInput$Gen== "16S")
#n18S <- subset(primerInput, primerInput$Gen== "18S")
#n28S <- subset(primerInput, primerInput$Gen== "28S")
#COI <-  subset(primerInput, primerInput$Gen== "COI")
#plants <- subset(primerInput, primerInput$Target== "Viridiplantae")
