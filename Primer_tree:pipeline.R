library("doMC")
registerDoMC(10)
library("primerTree")
library("ggplot2")
library("reshape2")
library("plyr")
library("ape")
library("TmCalculator")
library("dplyr")



###############Some try outs############
#Dice11_14_COI <- search_primer_pair(name='Dice11F_53Mod-94_F.Dice14R_53Mod-93_R', 
#                                   'GATCGDAAATTTDGTWCWGCG', 'CCHACMRTAAACATATGATGCC')
#plot(Dice11_14_COI, "genus")
#get_taxonomy(Dice11_14_COI)

#plot_tree(Dice11_14_COI$tree, type = "unrooted", main = "Genus", guide_size = NULL,
#          rank = "genus", taxonomy = Dice11_14_COI$taxonomy, size = 2, legend_cutoff = 60)

#plot_tree_ranks(Dice11_14_COI$tree, Dice11_14_COI$taxonomy ,ranks=c('phylum','family','genus'), size=3)

# Show the counts for each length
#seq_lengths(Dice11_14_COI)
# Plot the distribution of lengths
#seqL_Dice11_14 <- seq_lengths(Dice11_14_COI)
#barplot(seqL_Dice11_14,
#        main = "Frequency of sequence lengths for Dice11F_14R primers",
#        xlab="Amplicon length (in bp)",
#        ylab=("Frequency"))


tufA <- search_primer_pair(name='tufA', 
                                             'TGAAACAGAAMAWCGTCATTATGC', 'CCTTCNCGAATMGCRAAWCGC')

plot(Mamm1F_Mamm1R)
#plot_tree_ranks(Mamm1F_Mamm1R$tree, Mamm1F_Mamm1R$taxonomy ,ranks=c('phylum','family','genus'), size=2)
#seq_lengths(Mamm1F_Mamm1R)

Amoebozoa<- search_primer_pair(name='Amoebozoa', 'GAATTGACGGAAGGGCACAC', 'GCCCYRTCTAAGGGCATCAC', num_aligns = 1000)

Entamoeba<- search_primer_pair(name='Entamoeba', 'GTAGTGACGACAAATAACTCTTG', 'TTTCGTTCTTGATTAATGAACG')

?search_primer_pair

get_taxonomy(Amoebozoa)

Amoebozoa$taxonomy
plot_tree_ranks(Amoebozoa$tree, Amoebozoa$taxonomy ,ranks=c('phylum','family','genus'), size=5)
seq_lengths(Amoebozoa)

barplot(seq_lengths(Amoebozoa)) 

plot(Amoebozoa)

?plot()

#################Evaluate all primer pairs#################
primerInput <- read.csv("~/AA_Primer_evaluation/primerInput.csv")

#Remove duplicate rows based on one or more column values (require dplyr): 
primerInput <- primerInput %>% distinct(Primer_name, .keep_all = TRUE)
write.csv(primerInput, "~/AA_Primer_evaluation/primerInputUnique.csv")

#primers <- read.csv("~/AA_Primer_evaluation/primerInputUnique.csv")
#primers$X <- NULL

##Primer characteristics (Tm and GC)
primerInput$Seq_F <- as.vector(primerInput$Seq_F)
primerInput$Seq_R <- as.vector(primerInput$Seq_R)

primerInput <- primerInput[-c(97,103,106),] ###Eliminating primers with inosine ... To make the loop work :S 

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
primerInput <- merge(primerInput, output, by= "Primer_name", all= T)
primerInput$X <- NULL

write.csv(primerInput, "~/AA_Primer_evaluation/primerOutputUnique.csv")


###PrimerTree pipeline

doneFiles <- list.files(path = "~/AA_Primer_evaluation/output/primerTreeObj/", pattern = "_primerTree.RData")
doneFiles <- as.character(lapply(doneFiles, function(x) gsub("_primerTree.RData", "", x)))

for (i in 1:nrow(primers))
{
  
  Fwd <- as.character(primers[i, 4])
  Rev <- as.character(primers[i, 5])
  min_len <- primers[i, 6]
  max_len <- primers[i, 7]
  
  lab <- as.character(primers[i, 1], sep = "") ##labels 
  
  print(lab)
  
  if(lab%in% doneFiles) {
    
    load(paste("~/AA_Primer_evaluation/output/primerTreeObj/", lab, "_primerTree.RData", sep = ""))
  }else{
    
    primTree <- search_primer_pair(forward = Fwd,
                                   reverse = Rev,
                                   num_aligns = 10000,
                                   num_permutations = 150)
    
    assign(lab, primTree)
    save(list = (as.character(lab)), file = paste("~/AA_Primer_evaluation/output/primerTreeObj/", lab, "_primerTree.RData", sep = ""))
    
    if(length(names(get(lab))) == 9) {
      taxonomy <- get(lab)$taxonomy
      
      write.table(taxonomy, file = paste("~/AA_Primer_evaluation/output/taxonomy/", lab, "_taxonomy.txt", sep = ""),
                  quote = F,
                  sep = ",",
                  col.names = T,
                  row.names = F)
      
      sequences <- get(lab)$sequence
      write.dna(sequences, file = paste("~/AA_Primer_evaluation/output/sequences/", lab, "_sequences.fasta", sep = ""),
                format = "fasta")
      
      if(length(names(get(lab)$sequence))>5) {
        png(filename = paste("~/AA_Primer_evaluation/output/figures/", lab, "_primerTree.png", sep = ""),
            width = 4000,
            height = 4000,
            res = 300)
        
        print(plot(get(lab), size= 1, main= paste(lab, " primers", sep = "")))
        
        dev.off()
        
        seqLength <- as.data.frame(seq_lengths(get(lab)))
        colnames(seqLength) <- c("Lenght", "Count")
        seqLength$Lenght <- as.numeric(as.character(seqLength$Lenght))
        
        png(filename = paste("~/AA_Primer_evaluation/output/figures/", lab, "_SeqLens.png", sep = ""), 
            width = 4000, 
            height = 4000, 
            res = 300)
        print(
          ggplot(seqLength, aes(x = Length, y = Count)) + 
            geom_bar(stat = "identity") +
            theme(text = element_text(size = 20)) +
            ggtitle(paste(lab, " primers amplicon length\n", sep = ""))
        )
        dev.off()
        
      }
    }
  }
  
  rm(list = ls(pattern = lab))
  
}



