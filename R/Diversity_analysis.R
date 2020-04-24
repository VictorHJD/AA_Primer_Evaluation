##Phyloseq pipeline for diversity analysis by amplicon

library("ggplot2")
library("reshape")
library("phyloseq")
library("data.table")
library("parallel")
library("microbiome")
library("tidyverse")
library("plyr")
library("gridExtra")
library("grid")
library("vegan")
library("FactoMineR")
library("factoextra")
library("corrplot")

##Functions
sumSeqByTax <- function(PS.l, rank) {
  lapply(PS.l, function(x){ 
    counts <- data.frame(cbind(asvCount=colSums(otu_table(x)), tax_table(x)))
    counts$asvCount <- as.numeric(as.character(counts$asvCount))
    tapply(counts$asvCount, counts[, rank], sum)})
}

if(!exists("primerInput")){
  source("~/GitProjects/AA_Primer_Evaluation/R/General_Primer_information.R") ##   
}

if(!exists("PS.l")){
  PS.l<- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/MergedPhyloSeqList.Rds") ##   
}

Distribution <- F ##Ploting distribution by amplicon 
Rarefaction_curves <- F ##Ploting rarefaction curves by amplicon 
Beta_diversity_plots<- F ##Ploting PCoA for each marke
Phylum<- T ##Beta diversity analysis to phylum level
Genus<- F ##Beta diversity analysis to genus level
Figures<- F ##Activate to save new versions of the figures 

###Non rarefied raw counts by amplicon
rawcounts <- as.data.frame(unlist(lapply(PS.l, function(x){
  data.frame(cbind(asvCount=colSums(data.frame(rowSums(otu_table(x))))))
})))

rawcounts[,2] <- rownames(rawcounts)
rownames(rawcounts) <- c(1:nrow(rawcounts))
colnames(rawcounts) <- c("Reads_non_rare","Primer_comb_ID")
rawcounts <- data.frame(Primer_comb_ID = rawcounts$Primer_comb_ID, Reads_non_rare = rawcounts$Reads_non_rare) 
rawcounts$Primer_comb_ID <- gsub(".asvCount", "\\1", rawcounts$Primer_comb_ID)

##Create a table with total number of ASVs per primer pair
#write.csv(rawcounts, file = "~/AA_Primer_evaluation/ASVs_per_Primer_Pair.csv")

###Rarefaction curves by amplicon 
###Summary by amplicon
s.list <- lapply(PS.l, function (x) {
  summary(sample_sums(x))
})

###Distribution of reads by amplicon
if(Distribution){
  pdf("~/AA_Primer_evaluation/Figures/Manuscript/Supplementary_1.pdf", 
    width=8, height=10, onefile=T)
for(i in 1:length(PS.l)) {
  y<- names(PS.l)
  hist<- histogram(rowSums(otu_table(PS.l[[i]])),
          main = paste("Histogram of", y[i]),
          xlab = "Read count",
          ylab= "Percent of samples",
          freq = T)
  plot(hist)
}
dev.off()
}
###Rarefaction curves 
if(Rarefaction_curves){
  rarecurv.l <- lapply(PS.l, function (x) {
  vegan::rarecurve(otu_table(x),
                   label = F)
})

#pdf("~/AA_HMHZ/Primers/Supplementary_2.pdf", 
#    width=8, height=10, onefile=T)
#rarecurv.l
#dev.off()
}
## Eliminate samples with no counts for each primer
##Check when each marker reached the species saturation (from rarecurve!) for now just eliminate 0
PS.l.High <- lapply(PS.l, function (x) {
  prune_samples(sample_sums(x)>0, x)
})

###Non rarefied high counts by amplicon
highcounts <- as.data.frame(unlist(lapply(PS.l.High, function(x){
  data.frame(cbind(asvCount=colSums(data.frame(rowSums(otu_table(x))))))
})))

highcounts[,2] <- rownames(highcounts)
rownames(highcounts) <- c(1:nrow(highcounts))
colnames(highcounts) <- c("Reads_high","Primer_comb_ID")
highcounts <- data.frame(Primer_comb_ID = highcounts$Primer_comb_ID, Reads_high = highcounts$Reads_high) 
highcounts$Primer_comb_ID <- gsub(".asvCount", "\\1", highcounts$Primer_comb_ID)


###Rarefied data
##Each amplicon has their own depth, so rarefy to the total sum of reads by ASV by amplicon 
PS.l.rar <- lapply(PS.l.High, function (x) {
  rarefy_even_depth(x, rngseed=1234, sample.size= rowSums(otu_table(x)), replace=F) #mean(sample_sums(x))
})

### Rarefied sequence counts
rarecounts <- as.data.frame(unlist(lapply(PS.l.rar, function(x){
  data.frame(cbind(asvCount=colSums(data.frame(rowSums(otu_table(x))))))
})))

rarecounts[,2] <- rownames(rarecounts)
rownames(rarecounts) <- c(1:nrow(rarecounts))
colnames(rarecounts) <- c("Reads_rare","Primer_comb_ID")
rarecounts <- data.frame("Primer_comb_ID" = rarecounts$Primer_comb_ID, Reads_rare = rarecounts$Reads_rare) 
rarecounts$Primer_comb_ID <- gsub(".asvCount", "\\1", rarecounts$Primer_comb_ID)

##Alpha diversity 
alphaDiv <- data.frame()
for (i in 1:length(PS.l.rar)){
  a <- data.frame()
  b <- data.frame(estimate_richness(PS.l.rar[[i]], measures = c("Observed","Chao1", "Shannon")))
  a<-colMeans(b)
  alphaDiv<- rbind(alphaDiv, a)
}
alphaDiv[,5] <- names(PS.l.rar)
colnames(alphaDiv) <- c("Observed","Chao1","se_Chao1", "Shannon", "Primer_comb_ID")

##Beta diversity by marker considering each sample individually
ord.l <-lapply(PS.l.rar, function (x) {
  ordinate(x, method="PCoA", distance="bray")
})

if(Beta_diversity_plots){
  pdf("~/AA_Primer_evaluation/Figures/Manuscript/Supplementary_3.pdf", 
    width=10, height=10, onefile=T)
for (i in 1:length(PS.l.rar)) {
  x<- names(PS.l.rar)
  pca <- plot_ordination(PS.l.rar[[i]], ord.l[[i]],  type="taxa", color="phylum")+
    theme_bw()+
    geom_point(size=5, alpha= 0.75)+
    labs(title = x[i])
  print(pca)
}
dev.off()
}
###Beta diverisity without considering individual samples but considering each primer combination as an idividual community
##Variables= Taxa 
##Individuals= Primer combinations 
## Put all the infromation in the right format
if(Phylum){
readNumByPhylum <- sumSeqByTax(PS.l = PS.l.rar, rank = "phylum")

AbPhy <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByPhylum)) ### Start a loop: fro every element in the list ...
{ 
  phyla <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByPhylum[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    phyla[1,1] <- 0    ### Add a zero in the first column
    phyla[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    phyla <- as.data.frame((readNumByPhylum[[i]]))  ###Make a data frame with the data included in each element of the list 
    phyla[,2] <- rownames(phyla) ### And use the rownames as information of the second column 
  }
  
  phyla[,3] <- names(readNumByPhylum)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(phyla) <- c("ASV_count", "Phyla", "Primer_comb_ID") ### change the names for the columns 
  AbPhy <- rbind(AbPhy, phyla) ### Join all the "individual" data frames into the final data frame 
  
}   ### close loop

rownames(AbPhy) <- c(1:nrow(AbPhy)) ### change the rownames to consecutive numbers 
AbPhy <- data.frame(Primer_comb_ID = AbPhy$Primer_comb_ID, Phyla = AbPhy$Phyla, ASV_count = AbPhy$ASV_count) ###change the order of the columns
AbPhy$ASV_count <- as.numeric(AbPhy$ASV_count)

AbPhy %>%
  group_by(Primer_comb_ID) %>% 
  mutate(Total_count = sum(ASV_count))%>%
  mutate(Relative_abundance= ASV_count/Total_count)-> AbPhy ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

##This value represent the relative abundance respect to the total amount of reads
#Relative_abundance = AbPhy$ASV/AbPhy$Total_ASV ### create a vector with the result of the operation 

#AbPhy[,5] <- Relative_abundance ### And put it in the same order in the colum 5

#colnames(AbPhy)[5] <- "Relative_abundance" ### Change the name of the column 

AbPhy$Primer_comb_ID <- gsub(pattern = " ", replacement = "", x = AbPhy$Primer_comb_ID) ### check if the primer names have extra spaces

AbPhy$Primer_comb_ID <- gsub(pattern = "-", replacement = "_", x = AbPhy$Primer_comb_ID)

AbPhy <- merge(AbPhy, primerInput, by= "Primer_comb_ID") ###merge the selected information with the origial data frame created 

AbPhy <- plyr::join(AbPhy, rarecounts, by= "Primer_comb_ID")

##Create a REAL relative abundance value respect the total amount of reads per amplicon
Rel_abund_amp = AbPhy$ASV/AbPhy$Reads ### create a vector with the result of the operation 

AbPhy[,25] <- Rel_abund_amp ### And put it in the same order in the colum 25

colnames(AbPhy)[25] <- "Rel_abund_amp" ### Change the name of the column 

##Prepair matrix for PCA analysis
foo<- AbPhy%>%select(1,2,25)
foo<- reshape(foo, v.names = "Rel_abund_amp", timevar = "Phyla", idvar = "Primer_comb_ID",direction = "wide")
colnames(foo) <- gsub("Rel_abund_amp.", "\\1", colnames(foo)) ##Remove "Rel_abund_amp"
##Transform NA to 0
foo[is.na(foo)]<- 0
rownames(foo)<-foo[,1]
}

if(Genus){
  readNumByGenus <- sumSeqByTax(PS.l = PS.l.rar, rank = "genus")
  
  AbGen <- data.frame() ###Create the data frame 
  
  for (i in 1: length(readNumByGenus)) ### Start a loop: fro every element in the list ...
  { 
    genus <- data.frame() #### make an individual data frame ...
    
    if(nrow(as.data.frame(readNumByGenus[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
    {
      genus[1,1] <- 0    ### Add a zero in the first column
      genus[1,2] <- "NA" ### And a NA in the second column.
    }else               ### For the rest of the elements: 
    {
      genus <- as.data.frame((readNumByGenus[[i]]))  ###Make a data frame with the data included in each element of the list 
      genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
    }
    
    genus[,3] <- names(readNumByGenus)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
    colnames(genus) <- c("ASV", "Genus", "Primer_comb_ID") ### change the names for the columns 
    AbGen <- rbind(AbGen, genus) ### Join all the "individual" data frames into the final data frame 
    
  }   ### close loop
  
  rownames(AbGen) <- c(1:nrow(AbGen)) ### change the rownames to consecutive numbers 
  AbGen <- data.frame(Primer_comb_ID = AbGen$Primer_comb_ID, Genus = AbGen$Genus, ASV = AbGen$ASV) ###change the order of the columns
  AbGen$ASV <- as.numeric(AbGen$ASV)
  
  AbGen %>%
    group_by(Primer_comb_ID) %>% 
    mutate(Total_ASV = sum(ASV))%>%
    mutate(Relative_abundance= ASV/Total_ASV)-> AbGen ### Add a new variable that will contain the sum of all the sequencing reads by primer pair
  
  ##This value represent the relative abundance respect to the total amount of reads
  #Relative_abundance = AbGen$ASV/AbGen$Total_ASV ### create a vector with the result of the operation 
  
  #AbGen[,5] <- Relative_abundance ### And put it in the same order in the colum 5
  
  #colnames(AbGen)[5] <- "Relative_abundance" ### Change the name of the column 
  
  AbGen$Primer_comb_ID <- gsub(pattern = " ", replacement = "", x = AbGen$Primer_comb_ID) ### check if the primer names have extra spaces
  
  AbGen$Primer_comb_ID <- gsub(pattern = "-", replacement = "_", x = AbGen$Primer_comb_ID)
  
  AbGen <- merge(AbGen, primerInput, by= "Primer_comb_ID") ###merge the selected information with the origial data frame created 
  
  AbGen <- plyr::join(AbGen, rarecounts, by= "Primer_comb_ID")
  
  ##Create a REAL relative abundance value respect the total amount of reads per amplicon
  Rel_abund_amp = AbGen$ASV/AbGen$Reads ### create a vector with the result of the operation 
  
  AbGen[,25] <- Rel_abund_amp ### And put it in the same order in the colum 25
  
  colnames(AbGen)[25] <- "Rel_abund_amp" ### Change the name of the column 
  
  ##Prepair matrix for PCA analysis
  foo<- AbGen%>%select(1,2,25)
  foo<- reshape(foo, v.names = "Rel_abund_amp", timevar = "Genus", idvar = "Primer_comb_ID",direction = "wide")
  colnames(foo) <- gsub("Rel_abund_amp.", "\\1", colnames(foo)) ##Remove "Rel_abund_amp"
  ##Transform NA to 0
  foo[is.na(foo)]<- 0
  rownames(foo)<-foo[,1]
}

##Merge with primer information 
foo.primer<-join(foo, primerInput, by= "Primer_comb_ID")
rownames(foo.primer)<-foo.primer[,1]

##Add rarecounts
rarecounts<- as.data.frame(rarecounts)
rownames(rarecounts)<-rarecounts[,1]

foo.primer<-join(foo.primer, rarecounts, by= "Primer_comb_ID")
rownames(foo.primer)<- rownames(foo)

foo[,1]<- NULL
foo.primer[,1]<- NULL

######### Bray-Curtis dissimilarity estimation (Figure 2) ########
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
  select("Primer_comb_ID", "Gen", "Region") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Primer_comb_ID

col_groups$Primer_comb_ID<- NULL

colour_groups <- list( Gen= c("12S"= "#E3DAC9", "16S"= "pink","18S"= "#440154FF", "28S"= "#21908CFF", 
                              "COI"= "#FDE725FF", "ITS"= "#C46210", "rbcL"= "#D0FF14"),
                       Region= c("V1-V2"= "#8DD3C7", "V3-V4"= "#FFFFB3", "V4"= "#BEBADA", "V4-V5"= "#FB8072", "V6-V7"= "#80B1D3",
                                 "V6-V8"="#FDB462", "V7-V8"= "#B3DE69", "V8-V9"= "#FCCDE5", "V9"= "#D9D9D9",
                                 "D2"="#D95F02", "D3"= "#7570B3", "Folmer"= "#E7298A", "ITS1"= "#A6761D", "Chloroplast"= "#66A61E", "Mitochondrial"= "#E6AB02"))
require(pheatmap)
require(viridis)
BCheatmap <- pheatmap(foo.braycurt, 
                      color = plasma(100),
                      border_color = NA,
                      annotation_col = col_groups, 
                      #annotation_row = col_groups,
                      annotation_colors = colour_groups,
                      #cutree_rows = 2,
                      #cutree_cols = 2,
                      show_rownames = F,
                      show_colnames = F,
                      main= "Bray-Curtis dissimilarity among primers")

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_2.pdf", width = 10, height = 8)
BCheatmap
dev.off()

####PCA and composition plots (Figure 3 and supplementary)####
##Let's do a PCA of individuals
foo.pca<- PCA(foo, graph = T)
foo.eig<- get_eigenvalue(foo.pca)
fviz_eig(foo.pca, addlabels = TRUE, ylim = c(0, 25))

fviz_pca_ind(foo.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE) # Avoid text overlapping (slow if many points)

###Color them by gene
a<- fviz_pca_ind(foo.pca,
             pointsize = "cos2", ##It is possible to size it by number of reads pointsize = log10(foo.primer$Total_reads)
             pointshape = 21,
             geom.ind = "point", # show points only (but not "text")
             col.ind = "black", # color by groups
             fill.ind =  foo.primer$Gen,
             palette = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"),
             alpha.ind = 0.9,
             addEllipses = F, # Concentration ellipses
             legend.title = "Gen")+
  labs(tag = "A)")

###Color them by sequencing reads
b<- fviz_pca_ind(foo.pca,
             pointsize = log10(foo.primer$Reads), 
             geom.ind = "point", # show points only (but not "text")
             col.ind = foo.primer$Reads, # color by groups
             gradient.cols = c("blue", "red"),
             alpha.ind = 0.7,
             addEllipses = F, # Concentration ellipses
             legend.title = "Sequencing \nreads")+
  labs(tag = "B)")

##Color need to be changed manually according to righ order, annoying but let's find a solution later:S 
c<- fviz_contrib(foo.pca, choice = "ind", axes = 1:2, top = 10,
             fill = c("pink", "pink", "pink", "#21908CFF", "pink", "pink", 
                      "#440154FF", "#440154FF", "#440154FF", "#440154FF"), 
             color =  c("pink", "pink", "pink", "#21908CFF", "pink", "pink", 
                        "#440154FF", "#440154FF", "#440154FF", "#440154FF")) +
  labs(tag = "C)")+
  theme(axis.text.x = element_text(angle=90))

###PCA variables
d<- fviz_pca_var(foo.pca, col.var = "cos2", select.var =  list(contrib = 10), ##top ten taxa contributing
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

###Visualize quality of representation in each dimension
##Variable results
#foo.var <- get_pca_var(foo.pca)
#corrplot(foo.var$cos2, is.corr=FALSE) #Alternative

##Plot contribution
#fviz_contrib(foo.pca, choice = "var", axes = 1:2, top = 25)

###Clear difference between Bacterial and Eukaryotic... Let's split them

foo.euk<- foo.primer[foo.primer$Gen != "16S",]
foo.bac<- foo.primer[foo.primer$Gen == "16S",]

###Let's do a PCA for Eukaryotes 
if(Phylum){foo.euk2<- foo.euk[,1:45]}
if(Genus){foo.euk2<- foo.euk[,1:2067]}

foo.euk2.pca<- PCA(foo.euk2, graph = T)
foo.euk2.eig<- get_eigenvalue(foo.euk2.pca)
fviz_eig(foo.euk2.pca, addlabels = TRUE, ylim = c(0, 25))

fviz_pca_ind(foo.euk2.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE) # Avoid text overlapping (slow if many points)

aE<- fviz_pca_ind(foo.euk2.pca,
                 pointsize = "cos2", ##It is possible to size it by number of reads pointsize = log10(foo.primer$Total_reads)
                 pointshape = 21,
                 geom.ind = "point", # show points only (but not "text")
                 col.ind = "black", # color by groups
                 fill.ind =  foo.euk$Gen,
                 palette = c("#E3DAC9","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"),
                 addEllipses = F, # Concentration ellipses
                 legend.title = "Gen")+
  labs(tag = "A)")

bE<- fviz_pca_ind(foo.euk2.pca,
                 pointsize = log10(foo.euk$Reads), 
                 geom.ind = "point", # show points only (but not "text")
                 col.ind = foo.euk$Reads, # color by groups
                 gradient.cols = c("blue", "red"),
                 alpha.ind = 0.7,
                 addEllipses = F, # Concentration ellipses
                 legend.title = "Sequencing \nreads")+
  labs(tag = "B)")

cE<- fviz_contrib(foo.euk2.pca, choice = "ind", axes = 1:2, top = 10,
                 fill = c("#21908CFF", "#440154FF", "#440154FF", "#FDE725FF", "#440154FF",
                          "#E3DAC9", "#FDE725FF","#440154FF", "#440154FF",  "#440154FF"),
                 color = c("#21908CFF", "#440154FF", "#440154FF", "#FDE725FF", "#440154FF",
                           "#E3DAC9", "#FDE725FF","#440154FF", "#440154FF",  "#440154FF"))+
  labs(tag = "C)") +
  theme(axis.text.x = element_text(angle=90))

dE<- fviz_pca_var(foo.euk2.pca, col.var = "cos2", select.var = list(contrib = 15),
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                 repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

###Now do if for Prokaryotes 
if(Phylum){foo.bac2<- foo.bac[,1:45]}
if(Genus){foo.bac2<- foo.bac[,1:2067]}

foo.bac2.pca<- PCA(foo.bac2, graph = T)
foo.bac2.eig<- get_eigenvalue(foo.bac2.pca)
fviz_eig(foo.bac2.pca, addlabels = TRUE, ylim = c(0, 50))

fviz_pca_ind(foo.bac2.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE) # Avoid text overlapping (slow if many points)

aB<- fviz_pca_ind(foo.bac2.pca,
                  pointsize = "cos2", ##It is possible to size it by number of reads pointsize = log10(foo.primer$Total_reads)
                  pointshape = 21,
                  geom.ind = "point", # show points only (but not "text")
                  col.ind = "black", # color by region
                  fill.ind =  foo.bac$Region,
                  palette = c("#E6AB02", "#8DD3C7", "#FFFFB3", "#80B1D3"),
                  alpha.ind = 0.9,
                  addEllipses = F, # Concentration ellipses
                  legend.title = "Region")+
  labs(tag = "A)")

bB<- fviz_pca_ind(foo.bac2.pca,
                  pointsize = log10(foo.bac$Reads), 
                  geom.ind = "point", # show points only (but not "text")
                  col.ind = foo.bac$Reads, # color by groups
                  gradient.cols = c("blue", "red"),
                  alpha.ind = 0.7,
                  addEllipses = F, # Concentration ellipses
                  legend.title = "Sequencing \nreads")+
  labs(tag = "B)")

cB<- fviz_contrib(foo.bac2.pca, choice = "ind", axes = 1:2, fill = "pink", color = "pink")+
  labs(tag = "C)")+
  theme(axis.text.x = element_text(angle=90))

dB<- fviz_pca_var(foo.bac2.pca, col.var = "cos2", select.var = list(contrib = 10),
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                  repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

###Let's do it just for 18S primers

foo.18S<- foo.primer[foo.primer$Gen == "18S",]

if(Phylum){foo.18S2<- foo.18S[,1:45]}
if(Genus){foo.18S2<- foo.18S[,1:2067]}

foo.18S2.pca<- PCA(foo.18S2, graph = T)
foo.18S2.eig<- get_eigenvalue(foo.18S2.pca)
fviz_eig(foo.18S2.pca, addlabels = TRUE, ylim = c(0, 50))

fviz_pca_ind(foo.18S2.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE) # Avoid text overlapping (slow if many points)

a18<- fviz_pca_ind(foo.18S2.pca,
                   pointsize = "cos2", ##It is possible to size it by number of reads pointsize = log10(foo.primer$Total_reads)
                   pointshape = 21,
                   geom.ind = "point", # show points only (but not "text")
                   col.ind = "black", # color by region
                   fill.ind =  foo.18S$Region,
                   palette = c("#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9"),
                   addEllipses = F, # Concentration ellipses
                   legend.title = "Region")+
  labs(tag = "A)")

b18<- fviz_pca_ind(foo.18S2.pca,
                   pointsize = log10(foo.18S$Reads), 
                   geom.ind = "point", # show points only (but not "text")
                   col.ind = foo.18S$Reads, # color by groups
                   gradient.cols = c("blue", "red"),
                   alpha.ind = 0.7,
                   addEllipses = F, # Concentration ellipses
                   legend.title = "Sequencing \nreads")+
  labs(tag = "B)")

c18<- fviz_contrib(foo.18S2.pca, choice = "ind", axes = 1:2, fill = "#440154FF", color = "#440154FF")+
  labs(tag = "C)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1))

d18<- fviz_pca_var(foo.18S2.pca, col.var = "cos2", select.var = list(contrib = 10),
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                   repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

##Variable results
#foo.var18 <- get_pca_var(foo.18S2.pca)
#corrplot(na.omit(foo.var18$cos2), is.corr=FALSE) #Alternative

###Composition plots by primer pair 
if(Phylum){
  AbPhy%>%
  select(1,2,14,25)%>%
  group_by(Primer_comb_ID) %>%
  mutate(Main_taxa= Rel_abund_amp>= 0.05) %>%
  mutate(Phyla= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(.$Phyla))) %>%
  arrange(Primer_comb_ID, desc(Phyla))-> CompPhy

e <- ggplot(data=CompPhy, aes(x= Primer_comb_ID, y= Rel_abund_amp, fill= Phyla)) +
        scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
                                     "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061","#969696"))+
        geom_bar(aes(), stat="identity", position="fill") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
        labs(x = "Primer combination ID", y= "Relative abundance", tag = "E)")+
        guides(fill= guide_legend(nrow = 8))

eB <- ggplot(data=subset(CompPhy, Gen=="16S"), aes(x= Primer_comb_ID, y= Rel_abund_amp, fill= Phyla)) +
  scale_fill_manual(values = c("#33A02C","#CAB2D6","#6A3D9A","#01665e","#969696"))+
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Primer combination ID", y= "Relative abundance", tag = "E)")

eE <- ggplot(data=subset(CompPhy, Gen!="16S"), aes(x= Primer_comb_ID, y= Rel_abund_amp, fill= Phyla)) +
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00", 
                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061","#969696"))+
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Primer combination ID", y= "Relative abundance", tag = "E)")+
  guides(fill= guide_legend(nrow = 8))

e18 <- ggplot(data=subset(CompPhy, Gen=="18S"), aes(x= Primer_comb_ID, y= Rel_abund_amp, fill= Phyla)) +
        scale_fill_manual(values = c("#1F78B4","#B2DF8A","#FB9A99","#E31A1C","#FDBF6F",
                               "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#eddc12","#01665e","#053061","#969696"))+
        geom_bar(aes(), stat="identity", position="fill") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
        labs(x = "Primer combination ID", y= "Relative abundance", tag = "E)")+
        guides(fill= guide_legend(nrow = 8))
}

if(Genus){
  AbGen%>%
    select(1,2,14,25)%>%
    group_by(Primer_comb_ID) %>%
    mutate(Main_taxa= Rel_abund_amp>= 0.1) %>%
    mutate(Genus= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(.$Genus))) %>%
    arrange(Primer_comb_ID, desc(Genus))-> CompGen
  
  e <- ggplot(data=CompGen, aes(x= Primer_comb_ID, y= Rel_abund_amp, fill= Genus)) +
    scale_fill_manual(values = viridis(20))+
    geom_bar(aes(), stat="identity", position="fill") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
    labs(x = "Primer combination ID", y= "Relative abundance", tag = "E)")+
    guides(fill= guide_legend(nrow = 15))
  
  eB <- ggplot(data=subset(CompGen, Gen=="16S"), aes(x= Primer_comb_ID, y= Rel_abund_amp, fill= Genus)) +
    scale_fill_manual(values = c("#33A02C","#CAB2D6","#6A3D9A","#01665e","#969696"))+
    geom_bar(aes(), stat="identity", position="fill") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
    labs(x = "Primer combination ID", y= "Relative abundance", tag = "E)")
  
  eE <- ggplot(data=subset(CompGen, Gen!="16S"), aes(x= Primer_comb_ID, y= Rel_abund_amp, fill= Genus)) +
    scale_fill_manual(values = viridis(18))+
    geom_bar(aes(), stat="identity", position="fill") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
    labs(x = "Primer combination ID", y= "Relative abundance", tag = "E)")+
    guides(fill= guide_legend(nrow = 8))
  
  e18 <- ggplot(data=subset(CompGen, Gen=="18S"), aes(x= Primer_comb_ID, y= Rel_abund_amp, fill= Genus)) +
    scale_fill_manual(values = viridis(20))+
    geom_bar(aes(), stat="identity", position="fill") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
    labs(x = "Primer combination ID", y= "Relative abundance", tag = "E)")+
    guides(fill= guide_legend(nrow = 8))
}
###Compile and save all the figures

if(Figures){

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Supplementary_4.pdf", width = 10, height = 15)
grid.arrange(a, b, c, d, e, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)))
dev.off()

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Supplementary_5.pdf", width = 10, height = 15)
grid.arrange(aE, bE, cE, dE, eE, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)))
dev.off()

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Supplementary_6.pdf", width = 10, height = 15)
grid.arrange(aB, bB, cB, dB, eB, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)))
dev.off()

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_3.pdf", width = 10, height = 15)
grid.arrange(a18, b18, c18, d18, e18, widths = c(1, 1), layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 5)))
dev.off()
}