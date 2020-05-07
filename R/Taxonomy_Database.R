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

##Load PrimerTree results 
primerTreeResults<- list.files(path = "~/AA_Primer_evaluation/output/taxonomy", pattern = ".csv", full.names = T)

lapply(primerTreeResults, function(x) {
  objname<- gsub("/home/victor/AA_Primer_evaluation/output/taxonomy/", "", x)
  objname<- gsub(".csv", "", objname)
  tmp <- read_csv(x)
  assign(objname, tmp, envir = .GlobalEnv)
})

#source("~/GitProjects/AA_Primer_Evaluation/R/In_silico_primer_evaluation.R")

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
write.table(ENA_18S_db, file = "/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/input/ENA_18S_db.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

##convert taxonomic IDs to taxonomy
Taxonomy_18S<- taxonomizr::getTaxonomy(as.character(taxID_18S$Taxonomic_ID), sqlFile = "/SAN/db/taxonomy/taxonomizr.sql")

Taxonomy_18S<- as.data.frame(Taxonomy_18S)
rownames(Taxonomy_18S) <- c(1:nrow(Taxonomy_18S))
Taxonomy_18S$Taxonomic_ID<- taxID_18S$Taxonomic_ID
Taxonomy_18S$Accession_18S<- NA

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
    tmp<- dplyr::slice(tmp, match(c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "taxId", "accession"), V2))
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
  geom_text(aes(24, cum.tax.comb18.insilico[25,2]+3000, label="Species"), color="red")+ 
  geom_hline(yintercept = Unique_taxa[7,1], color="red", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(24, cum.tax.comb18.insilico[25,3]+1000, label="Genera"), color="#0D0887FF") +
  geom_hline(yintercept = Unique_taxa[6,1], color="#0D0887FF", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(24, cum.tax.comb18.insilico[25,4]+400, label="Families"), color="#7E03A8FF") +
  geom_hline(yintercept = Unique_taxa[5,1], color="#7E03A8FF", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(24, cum.tax.comb18.insilico[25,5]+100, label="Orders"), color="#CC4678FF") +
  geom_hline(yintercept = Unique_taxa[4,1], color="#CC4678FF", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  geom_text(aes(24, cum.tax.comb18.insilico[25,6]+50, label="Classes"), color="#F0F921FF") +
  geom_hline(yintercept = Unique_taxa[3,1], color="#F0F921FF", linetype= "dashed")+
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(24, cum.tax.comb18.insilico[25,7]+10, label="Phyla"), color="#F89441FF") +
  geom_hline(yintercept = Unique_taxa[2,1], color="#F89441FF", linetype= "dashed")+
  scale_y_log10("Cummulative count of taxa") +
  scale_x_continuous("Number of primers considered \n (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), text = element_text(size=20)#,
        #axis.title.x = element_blank()
  )+
  labs(tag = "A)")#+
  #coord_cartesian(ylim = c(20, 4300), xlim = c(0,25))

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

##Prepair matrix for PCA analysis at phylum level (just for primers in manuscript)
alltaxa%>%
  dplyr::select(2,10)%>%
  group_by(Primer_comb_ID)%>%
  distinct(phylum, .keep_all = TRUE)%>%
  mutate(Abundance= 1)%>%
  dplyr::filter(!(Primer_comb_ID %in% c("Euk_18S_22", "Euk_18S_23", "Euk_18S_24")))-> foo

foo<- as.data.frame(foo)
foo<- reshape(foo, v.names = "Abundance",timevar= "phylum", idvar= "Primer_comb_ID", direction= "wide")
colnames(foo) <- gsub("Abundance.", "\\1", colnames(foo)) ##Remove "Rel_abund_amp"
foo[is.na(foo)]<- 0
rownames(foo)<-foo[,1]
foo[,1]<- NULL

require("vegan")
foo.matrix<- as.matrix(foo)
foo.braycurt<- vegdist(foo.matrix, method = "bray")
as.matrix(foo.braycurt)

##Subset matrix for parasites in silico
parasites.18S2.insilico<- foo.matrix[,c("Apicomplexa", "Nematoda", "Microsporidia", "Platyhelminthes")]

##Load matrix for parasites in vitro
parasites.18S2.invitro<- readRDS(file = "~/AA_Primer_evaluation/In_vitro_parasites_18S.RDS")
parasites.18S2.invitro<- as.matrix(parasites.18S2.invitro)

##Transform relative abundance data into presence/absence data
for(i in 1:ncol(parasites.18S2.invitro)){
  for(j in 1:nrow(parasites.18S2.invitro)){
if(parasites.18S2.invitro[j,i]>0){
  parasites.18S2.invitro[j,i]<-1
    } 
  }
}

##Estimate BC dissimilarity in vitro and in silico for parasites 
parasites.invitro.BC<- vegdist(parasites.18S2.invitro, method = "bray")
parasites.insilico.BC<- vegdist(parasites.18S2.insilico, method = "bray")

parasites.insilico.BC<- as.matrix(parasites.insilico.BC)
parasites.invitro.BC<- as.matrix(parasites.invitro.BC)

mantel(parasites.insilico.BC, parasites.invitro.BC, permutations = 1000)

##Plot BC dissimiliraty and geographic distance 
combi.dist<- data.frame(insilico= as.vector(parasites.insilico.BC),
                        invitro= as.vector(parasites.invitro.BC))

geo.plot<- ggplot(combi.dist, aes(insilico, invitro)) +
  geom_jitter(alpha=1, width=0.3, height=0, color= "lightblue") + 
  stat_smooth(se=T, method="lm") +
  scale_x_continuous("In silico Bray-Curtis dissimilarity between primers") +
  scale_y_continuous("In vitro Bray-Curtis dissimilarity between primers") +
  labs(tag = "B)")+
  theme_classic(base_size = 15, base_family = "Helvetica")



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
                       Region= c("V1-V2"= "#8DD3C7", "V1-V3"= "#009999","V3-V4"= "#FFFFB3", "V4"= "#BEBADA", "V4-V5"= "#FB8072", "V6-V7"= "#80B1D3",
                                 "V6-V8"="#FDB462", "V7-V8"= "#B3DE69", "V7-V9"= "#FC4E07","V8-V9"= "#FCCDE5", "V9"= "#D9D9D9"))
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

#pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_2_insilico.pdf", width = 10, height = 8)
pdf(file = "~/AA_Primer_evaluation/Figures/In_silico_evaluation/Preliminary/Figure_4_All.pdf", width = 10, height = 8)
BCheatmap.insilico
dev.off()

###PCA with Primer tree data 
foo$Primer_comb_ID<- rownames(foo)
foo.primer<-join(foo, primerInput, by= "Primer_comb_ID")
rownames(foo.primer)<- foo.primer$Primer_comb_ID
foo$Primer_comb_ID<- NULL

require("FactoMineR") 
foo.18S_Is.pca<- PCA(foo, graph = T)
foo.18S_Is.eig<- get_eigenvalue(foo.18S_Is.pca)

a18Is<- fviz_pca_ind(foo.18S_Is.pca,
                     pointsize = "cos2", ##It is possible to size it by number of reads pointsize = log10(foo.primer$Total_reads)
                     pointshape = 21,
                     geom.ind = "point", # show points only (but not "text")
                     col.ind = "black", # color by region
                     fill.ind =  foo.primer$Region,
                     palette = c("#8DD3C7", "#009999","#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FC4E07","#FCCDE5", "#D9D9D9"),
                     addEllipses = F, # Concentration ellipses
                     legend.title = "Region")+
  labs(tag = "A)")

b18IS<-fviz_pca_ind(foo.18S_Is.pca,
                    pointsize = log10(foo.primer$Expected), 
                    geom.ind = "point", # show points only (but not "text")
                    col.ind = foo.primer$Expected, # color by groups
                    gradient.cols = c("blue", "red"),
                    alpha.ind = 0.7,
                    addEllipses = F, # Concentration ellipses
                    legend.title = "Expected \namplicon size")+
  labs(tag = "B)")

c18Is<- fviz_contrib(foo.18S_Is.pca, choice = "ind", axes = 1:2, fill = "#440154FF", color = "#440154FF")+
  labs(tag = "C)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1))

d18Is<- fviz_pca_var(foo.18S_Is.pca, col.var = "cos2", select.var = list(contrib = 10),
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                     repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

##Composition bar plot is made based on data from In_silico_primer_evaluation.R (for now later will be necessary to re-structure codes :S)
#e18Is<-ggplot(data=Relative_abundance_18S, aes(x= Primer_name,y= Rel_abund, fill= phylum)) +
#  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
#                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061", #"#c00000",
#                               "darkolivegreen",#"#FFDB77FF",
#                               "lightgrey",#"#004CFFFF",
#                               "#969696"))+
#  geom_bar(aes(), stat="identity", position="fill") +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  labs(x = "Primer combination ID", y= "Relative abundance (In silico)", tag = "E)")+
  guides(fill= guide_legend(nrow = 8))

pdf(file = "~/AA_Primer_evaluation/Figures/In_silico_evaluation/Preliminary/Figure_PrimerTree.pdf", width = 10, height = 15)
grid.arrange(a18Is,b18IS, c18Is, d18Is, e18Is, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)))
dev.off()

