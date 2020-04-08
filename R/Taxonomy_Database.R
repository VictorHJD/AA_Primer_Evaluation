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

##Get number of cummulative unique by 18S primer combination in silico 

pre_cum_in_silico<- list()

for(i in names(alltaxa)){
  tmp <- data.frame() #### make an individual data frame ...
  {
    tmp<-as.data.frame(dplyr::distinct(alltaxa, !!as.name(i), .keep_all= T)) 
    tmp<- dplyr::select(tmp, c(!!as.name(i), "Primer_comb_ID"))
    tmp<- dplyr::group_by(tmp, Primer_comb_ID)
    tmp<- dplyr::summarise(tmp, n())
    if(colnames(tmp) == i) ### A tiny condition for the list "Primer_comb_ID"
    {
      colnames(tmp)<- c("Primer_comb_ID", "Primer_ID") ### Add a zero in the first column
    }else   
    colnames(tmp)<- c("Primer_comb_ID", i)
  }
  pre_cum_in_silico[[i]] <- tmp
}

pre_cum_in_silico<- lapply(pre_cum_in_silico, data.frame, row.names=NULL)

pre_cum_in_silico<- Reduce(function(x,y) merge(x, y, by = "Primer_comb_ID", all.x = TRUE, all.y = TRUE),
                           pre_cum_in_silico)

pre_cum_in_silico[is.na(pre_cum_in_silico)]<- 0 ##Transform false NAs to zeros

pre_cum_in_silico$Primer_ID<- NULL

#write.csv(pre_cum_in_silico, "~/AA_Primer_evaluation/Pre_cum_tax_comb18_in_silico.csv") ##save again just in case of change

##Create a real cummulative count of taxa by primer combination for cummulative curves
pre_cum_in_silico%>%
  mutate(cum.species= cumsum(species))%>%
  mutate(cum.genus= cumsum(genus))%>%
  mutate(cum.family= cumsum(family))%>%
  mutate(cum.order= cumsum(order))%>%
  mutate(cum.class= cumsum(class))%>%
  mutate(cum.phylum= cumsum(phylum))%>%
  mutate(Primer_name= as.integer(rownames(pre_cum_in_silico)))%>%
  dplyr::select(-(1:10))%>%
  dplyr::select(Primer_name, cum.species, cum.genus, cum.family, cum.order, cum.class, cum.phylum)-> cum.tax.comb18.insilico

rm(pre_cum_in_silico, tmp)

baseline.zero <- data.frame(0, 0, 0, 0, 0, 0, 0)
names(baseline.zero)<- c("Primer_name", "cum.species", "cum.genus", "cum.family", "cum.order", "cum.class", "cum.phylum")

cum.tax.comb18.insilico<- rbind(baseline.zero, cum.tax.comb18.insilico)

cum.plot.comb18.insilico <- ggplot(cum.tax.comb18.insilico) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(24, 4300, label="Species"), color="red")+ 
  geom_hline(yintercept = 4098, color="red", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(24, 2100, label="Genera"), color="#0D0887FF") +
  geom_hline(yintercept = 1988, color="#0D0887FF", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(24, 950, label="Families"), color="#7E03A8FF") +
  geom_hline(yintercept = 850, color="#7E03A8FF", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(24, 360, label="Orders"), color="#CC4678FF") +
  geom_hline(yintercept = 321, color="#CC4678FF", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  geom_text(aes(24, 115, label="Classes"), color="#F0F921FF") +
  geom_hline(yintercept = 102, color="#F0F921FF", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(24, 38, label="Phyla"), color="#F89441FF") +
  geom_hline(yintercept = 34, color="#F89441FF", linetype= "dashed")+
  scale_y_log10("Cummulative count of taxa") +
  scale_x_continuous("Number of primers considered \n (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), text = element_text(size=20)#,
        #axis.title.x = element_blank()
  )+
  labs(tag = "A)")+
  coord_cartesian(ylim = c(20, 4300), xlim = c(0,25))

#pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Pre_Figure_4.pdf", width = 8, height = 10)
#grid.arrange(cum.plot.comb18.insilico, cum.plot.comb18) 
#dev.off()

###Matrix with presence/absence (1/0) data for the in silico analysis
##Prepair matrix for PCA analysis at genus level
alltaxa%>%
  dplyr::select(6,10)%>%
  group_by(Primer_comb_ID)%>%
  distinct(genus, .keep_all = TRUE)%>%
  mutate(Abundance= 1)-> foo

foo<- as.data.frame(foo)
foo<- reshape(foo, v.names = "Abundance",timevar= "genus", idvar= "Primer_comb_ID", direction= "wide")
colnames(foo) <- gsub("Abundance.", "\\1", colnames(foo)) ##Remove "Rel_abund_amp"
foo[is.na(foo)]<- 0
rownames(foo)<-foo[,1]
foo[,1]<- NULL

foo.matrix<- as.matrix(foo)
foo.braycurt<- vegdist(foo.matrix, method = "bray")
as.matrix(foo.braycurt)

###Using pheatmap to include annotations 
foo.clust <- hclust(dist(foo.braycurt), method = "complete") ##Dendogram
require(dendextend)
as.dendrogram(foo.clust) %>%
  plot(horiz = TRUE)

foo.col <- cutree(tree = foo.clust, k = 2)
foo.col  <- data.frame(cluster = ifelse(test = foo.col  == 1, yes = "cluster 1", no = "cluster 2"))
foo.col$Primer_comb_ID <- rownames(foo.col)

foo.col <- merge(foo.col, primerInput, by="Primer_comb_ID", sort= F)

col_groups <- foo.col %>%
  dplyr::select("Primer_comb_ID", "Gen", "Region") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Primer_comb_ID

col_groups$Primer_comb_ID<- NULL

colour_groups <- list( Gen= c("18S"= "#440154FF"),
                       Region= c("V1-V2"= "#8DD3C7", "V3-V4"= "#FFFFB3", "V4"= "#BEBADA", "V4-V5"= "#FB8072", "V6-V7"= "#80B1D3",
                                 "V6-V8"="#FDB462", "V7-V8"= "#B3DE69", "V8-V9"= "#FCCDE5", "V9"= "#D9D9D9"))
require(pheatmap)
require(viridis)
BCheatmap.insilico <- pheatmap(foo.braycurt, 
                      color = plasma(100),
                      border_color = NA,
                      annotation_col = col_groups, 
                      #annotation_row = col_groups,
                      annotation_colors = colour_groups,
                      #cutree_rows = 2,
                      #cutree_cols = 2,
                      show_rownames = F,
                      show_colnames = F,
                      main= "Bray-Curtis dissimilarity among 18S primers in silico")

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_2_insilico.pdf", width = 10, height = 8)
BCheatmap.insilico
dev.off()

