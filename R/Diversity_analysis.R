##Phyloseq pipeline for diversity analysis by amplicon

library(ggplot2)
library(reshape)
library(phyloseq)
library(data.table)
library(parallel)
library(microbiome)
library(tidyverse)
library(plyr)
library(gridExtra)
library(grid)
library(vegan)
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

###Raw counts by amplicon
rawcounts <- as.data.frame(unlist(lapply(PS.l, function(x){
  data.frame(cbind(asvCount=colSums(data.frame(rowSums(otu_table(x))))))
})))

rawcounts[,2] <- rownames(rawcounts)
rownames(rawcounts) <- c(1:nrow(rawcounts))
colnames(rawcounts) <- c("Reads","Primer_name")
rawcounts <- data.frame(Primer_name = rawcounts$Primer_name, Reads = rawcounts$Reads) 
rawcounts$Primer_name <- gsub(".asvCount", "\\1", rawcounts$Primer_name)

###Rarefaction curves by amplicon 
###Summary by amplicon
s.list <- lapply(PS.l, function (x) {
  summary(sample_sums(x))
})

###Distribution of reads by amplicon
pdf("~/AA_Primer_evaluation/Figures/Supplementary_3.pdf", 
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

###Rarefaction curves 
rarecurv.l <- lapply(PS.l, function (x) {
  vegan::rarecurve(otu_table(x),
                   label = F)
})

#pdf("~/AA_HMHZ/Primers/RareCurv_by_amplicon.pdf", 
#    width=8, height=10, onefile=T)
#rarecurv.l
#dev.off()

## Eliminate samples with no counts for each primer
##Check when each marker reached the species saturation (from rarecurve!) for now just eliminate 0
PS.l.High <- lapply(PS.l, function (x) {
  prune_samples(sample_sums(x)>0, x)
})

###Rarefied data
##Each amplicon has their own depth... so rarefy to the max count sum (option 1) or ... 
#PS.l.rar <- lapply(PS.l.High, function (x) {
#  rarefy_even_depth(x, rngseed=1234, sample.size= max(sample_sums(x)), replace=F) #
#})

##... to the total sum of reads by ASV by amplicon (option 2)
PS.l.rar <- lapply(PS.l.High, function (x) {
  rarefy_even_depth(x, rngseed=1234, sample.size= rowSums(otu_table(x)), replace=F) #max(sample_sums(x))
})

rawcounts <- as.data.frame(unlist(lapply(PS.l.rar, function(x){
  data.frame(cbind(asvCount=colSums(data.frame(rowSums(otu_table(x))))))
})))

rawcounts[,2] <- rownames(rawcounts)
rownames(rawcounts) <- c(1:nrow(rawcounts))
colnames(rawcounts) <- c("Reads","Primer_name")
rawcounts <- data.frame(Primer_name = rawcounts$Primer_name, Reads = rawcounts$Reads) 
rawcounts$Primer_name <- gsub(".asvCount", "\\1", rawcounts$Primer_name)

##Normalize data 
PS.l.Nor <- lapply(PS.l, function (x) {
  transform_sample_counts(x, function(i) i/sum(i))
})

##Alpha diversity 

alphaDiv <- data.frame()

for (i in 1:length(PS.l.rar)){
  
  a <- data.frame()
  
  b <- data.frame(estimate_richness(PS.l.rar[[i]], measures = c("Observed","Chao1", "Shannon")))
  
  a<-colMeans(b)
  
  alphaDiv<- rbind(alphaDiv, a)
}

alphaDiv[,5] <- names(PS.l.rar)

colnames(alphaDiv) <- c("Observed","Chao1","se_Chao1", "Shannon", "Primer_name")


##Beta diversity
ord.l <-lapply(PS.l.rar, function (x) {
  ordinate(x, method="PCoA", distance="bray")
})

pdf("~/AA_Primer_evaluation/Figures/Supplementary_4.pdf", 
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

###"Beta diverisity" without considering individual samples 
##Variables= Taxa 
##Individuals= Primer combinations 
## Put all the infromation in the right format

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
  colnames(phyla) <- c("ASV", "Phyla", "Primer_name") ### change the names for the columns 
  AbPhy <- rbind(AbPhy, phyla) ### Join all the "individual" data frames into the final data frame 
  
}   ### close loop

rownames(AbPhy) <- c(1:nrow(AbPhy)) ### change the rownames to consecutive numbers 
AbPhy <- data.frame(Primer_name = AbPhy$Primer_name, Phyla = AbPhy$Phyla, ASV = AbPhy$ASV) ###change the order of the columns
AbPhy$ASV <- as.numeric(AbPhy$ASV)

AbPhy %>%
  group_by(Primer_name) %>% 
  mutate(Total_ASV = sum(ASV)) -> AbPhy ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

##This value represent the relative abundance respect to the total amount of reads
Relative_abundance = AbPhy$ASV/AbPhy$Total_ASV ### create a vector with the result of the operation 

AbPhy[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbPhy)[5] <- "Relative_abundance" ### Change the name of the column 

AbPhy$Primer_name <- gsub(pattern = " ", replacement = "", x = AbPhy$Primer_name) ### check if the primer names have extra spaces

AbPhy$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbPhy$Primer_name)

AbPhy <- merge(AbPhy, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

AbPhy <- plyr::join(AbPhy, rawcounts, by= "Primer_name")

##Create a REAL relative abundance value respect the total amount of reads per amplicon
Rel_abund_amp = AbPhy$ASV/AbPhy$Reads ### create a vector with the result of the operation 

AbPhy[,24] <- Rel_abund_amp ### And put it in the same order in the colum 24

colnames(AbPhy)[24] <- "Rel_abund_amp" ### Change the name of the column 


##Prepair matrix for PCA analysis
foo<- AbPhy%>%select(1,2,24)
foo<- reshape(foo, v.names = "Rel_abund_amp", timevar = "Phyla", idvar = "Primer_name",direction = "wide")
colnames(foo) <- gsub("Rel_abund_amp.", "\\1", colnames(foo)) ##Remove "Rel_abund_amp"
##Transform NA to 0
foo[is.na(foo)]<- 0
rownames(foo)<-foo[,1]

##Merge with primer information 
foo.primer<-join(foo, primerInput, by= "Primer_name")
rownames(foo.primer)<-foo.primer[,1]

##Add rawcounts
rawcounts<- as.data.frame(rawcounts)
rownames(rawcounts)<-rawcounts[,1]

foo.primer<-join(foo.primer, rawcounts, by= "Primer_name")
rownames(foo.primer)<- rownames(foo)

foo[,1]<- NULL
foo.primer[,1]<- NULL

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
             col.ind = foo.primer$Gen, # color by groups
             fill.ind =  foo.primer$Gen,
             palette = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"),
             alpha.ind = 0.5,
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
             legend.title = "Sequencing reads")+
  labs(tag = "B)")

##Color need to be changed manually according to righ order, annoying but let's find a solution later:S 
c<- fviz_contrib(foo.pca, choice = "ind", axes = 1:2, top = 10,
             fill = c("pink", "pink", "pink", "#21908CFF", "pink", "pink", 
                      "#440154FF", "#440154FF", "#440154FF", "#440154FF"), 
             color =  c("pink", "pink", "pink", "#21908CFF", "pink", "pink", 
                        "#440154FF", "#440154FF", "#440154FF", "#440154FF")) +
  labs(tag = "C)")+
  theme(axis.text.x = element_text(angle=90))

##Variable results
foo.var <- get_pca_var(foo.pca)

###PCA variables
#fviz_pca_var(foo.pca, col.var = "black")
d<- fviz_pca_var(foo.pca, col.var = "cos2", select.var =  list(contrib = 10), ##top ten taxa contributing
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

###Visualize quality of representation 
corrplot(foo.var$cos2, is.corr=FALSE) #Alternative

##Plot contribution
e<- fviz_contrib(foo.pca, choice = "var", axes = 1:2, top = 25)+
  labs(tag = "E)")

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_6.pdf", width = 16, height = 8)
grid.arrange(a, b, c, d, e, widths = c(1, 1, 1), layout_matrix = rbind(c(1, 3, 5),
                                                                      c(2, 4, 5)))
dev.off()

###Clear difference between Bacterial and Eukaryotic... Let's split them

foo.euk<- foo.primer[foo.primer$Gen != "16S",]
foo.bac<- foo.primer[foo.primer$Gen == "16S",]

###Let's do a PCA for Eukaryotes 
#foo.euk2<- foo.euk[,1:21] ##Rarefaction option 1 
foo.euk2<- foo.euk[,1:45] ##Rarefaction option 2

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
                 col.ind = foo.euk$Gen, # color by groups
                 fill.ind =  foo.euk$Gen,
                 palette = c("#E3DAC9","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"),
                 alpha.ind = 0.5,
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
                 legend.title = "Sequencing reads")+
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

eE<- fviz_contrib(foo.euk2.pca, choice = "var", axes = 1:2, top = 25)+
  labs(tag = "E)")

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_7.pdf", width = 16, height = 8)
grid.arrange(aE, bE, cE, dE, eE, widths = c(1, 1, 1), layout_matrix = rbind(c(1, 3, 5),
                                                                       c(2, 4, 5)))
dev.off()

###Now do if for Prokaryotes 
#foo.bac2<- foo.bac[,1:21] ##Rarefaction option 1
foo.bac2<- foo.bac[,1:45] ##Rarefaction option 2

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
                  col.ind = foo.bac$Region, # color by region
                  fill.ind =  foo.bac$Region,
                  palette ="Set1",
                  alpha.ind = 0.5,
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
                  legend.title = "Sequencing reads")+
  labs(tag = "B)")

cB<- fviz_contrib(foo.bac2.pca, choice = "ind", axes = 1:2, fill = "pink", color = "pink")+
  labs(tag = "C)")+
  theme(axis.text.x = element_text(angle=90))

dB<- fviz_pca_var(foo.bac2.pca, col.var = "cos2", select.var = list(contrib = 10),
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                  repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

eB<- fviz_contrib(foo.bac2.pca, choice = "var", axes = 1:2, top = 25)+
  labs(tag = "E)")

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_8.pdf", width = 16, height = 8)
grid.arrange(aB, bB, cB, dB, eB, widths = c(1, 1, 1), layout_matrix = rbind(c(1, 3, 5),
                                                                            c(2, 4, 5)))
dev.off()

###Eliminating the highly dominant primers, this plot is usful just when non-rarefied data is used 
#foo2<- foo[-c(7,13,18),]
#foo.primer2<- foo.primer[-c(7,13,18),]

#foo.pca2<- PCA(foo2, graph = T)
#foo.eig2<- get_eigenvalue(foo.pca2)
#fviz_eig(foo.pca2, addlabels = TRUE, ylim = c(0, 25))

#fviz_pca_ind(foo.pca2, pointsize = "cos2", 
#             pointshape = 21, fill = "#E7B800",
#             repel = TRUE) # Avoid text overlapping (slow if many points)

###Color them by gene
#a2<- fviz_pca_ind(foo.pca2,
#                 pointsize = "cos2", ##It is possible to size it by number of reads pointsize = log10(foo.primer$Total_reads)
#                 pointshape = 21,
#                 geom.ind = "point", # show points only (but not "text")
#                 col.ind = foo.primer2$Gen, # color by groups
#                 fill.ind =  foo.primer2$Gen,
#                 palette = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"),
#                 alpha.ind = 0.5,
#                 addEllipses = F, # Concentration ellipses
#                 legend.title = "Gen")+
#  labs(tag = "A)")

###Color them by sequencing reads
#b2<- fviz_pca_ind(foo.pca2,
#                 pointsize = log10(foo.primer2$Total_reads), 
#                 geom.ind = "point", # show points only (but not "text")
#                 col.ind = foo.primer2$Total_reads, # color by groups
#                 gradient.cols = c("blue", "red"),
#                 alpha.ind = 0.7,
#                 addEllipses = F, # Concentration ellipses
#                 legend.title = "Sequencing reads")+
#  labs(tag = "B)")

#c2<- fviz_contrib(foo.pca2, choice = "ind", axes = 1:2, top = 10,
#                 fill = c("#440154FF", "pink", "pink", "pink", "pink",
#                          "#440154FF", "#440154FF", "#440154FF", "#440154FF", "#440154FF"),
#                 color = c("#440154FF", "pink", "pink", "pink", "pink",
#                           "#440154FF", "#440154FF", "#440154FF", "#440154FF", "#440154FF"))+
#  labs(tag = "C)")+
#  theme(axis.text.x = element_text(angle=90))

##Variable results
#foo.var2 <- get_pca_var(foo.pca2)

###PCA variables
#fviz_pca_var(foo.pca, col.var = "black")
#d2<- fviz_pca_var(foo.pca2, col.var = "cos2", select.var =  list(contrib = 10), ##top ten taxa contributing
#                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#                 repel = TRUE )+ # Avoid text overlapping
#  labs(tag = "D)")

###Visualize quality of representation 
#corrplot(na.omit(foo.var2$cos2), is.corr=FALSE) #Alternative

##Plot contribution
#e2<- fviz_contrib(foo.pca2, choice = "var", axes = 1:2, top = 25)+
#  labs(tag = "E)")

#pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_9.pdf", width = 16, height = 8)
#grid.arrange(a2, b2, c2, d2, e2, widths = c(1, 1, 1), heights= c(2,2), layout_matrix = rbind(c(1, 3, 5),
#                                                                       c(2, 4, 5)))
#dev.off()

###Let's do it just for 18S primers

foo.18S<- foo.primer[foo.primer$Gen == "18S",]

#foo.18S2<- foo.18S[,1:21] ##Rarefaction option 1
foo.18S2<- foo.18S[,1:45] ##Rarefaction option 2

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
                   col.ind = foo.18S$Region, # color by region
                   fill.ind =  foo.18S$Region,
                   palette ="Set1",
                   alpha.ind = 0.5,
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
                   legend.title = "Sequencing reads")+
  labs(tag = "B)")

c18<- fviz_contrib(foo.18S2.pca, choice = "ind", axes = 1:2, fill = "#440154FF", color = "#440154FF")+
  labs(tag = "C)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1))

d18<- fviz_pca_var(foo.18S2.pca, col.var = "cos2", select.var = list(contrib = 10),
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                   repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

e18<- fviz_contrib(foo.18S2.pca, choice = "var", axes = 1:2, top = 25)+
  labs(tag = "E)")

##Variable results
#foo.var18 <- get_pca_var(foo.18S2.pca)
#corrplot(na.omit(foo.var18$cos2), is.corr=FALSE) #Alternative

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_10.pdf", width = 16, height = 8)
grid.arrange(a18, b18, c18, d18, e18, widths = c(1, 1, 1), layout_matrix = rbind(c(1, 3, 5),
                                                                                 c(2, 4, 5)))
dev.off()

##Bray-Curtis dissimilarity estimation

foo.matrix<- as.matrix(foo)
foo.braycurt<- vegdist(foo.matrix, method = "bray")
as.matrix(foo.braycurt)

##plot heatmap of Bray-Curtis dissimilarity with factoextra
BC <- fviz_dist(foo.braycurt, lab_size = 8,  gradient = list(low = "#FC4E07", high = "#00AFBB")) +
  labs(tag = "A)", fill= "Bray-Curtis\ndissimilarity")

###Using pheatmap to include annotations 
foo.clust <- hclust(dist(foo.braycurt), method = "complete") ##Dendogram
require(dendextend)
as.dendrogram(foo.clust) %>%
  plot(horiz = TRUE)

foo.col <- cutree(tree = foo.clust, k = 2)
foo.col  <- data.frame(cluster = ifelse(test = foo.col  == 1, yes = "cluster 1", no = "cluster 2"))
foo.col$Primer_name <- rownames(foo.col)

foo.col <- merge(foo.col, primerInput, by="Primer_name", sort= F)

col_groups <- foo.col %>%
  select("Primer_name", "Gen", "Region") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Primer_name

col_groups$Primer_name<- NULL

colour_groups <- list( Gen= c("12S"= "#E3DAC9", "16S"= "pink","18S"= "#440154FF", "28S"= "#21908CFF", 
                              "COI"= "#FDE725FF", "ITS"= "#C46210", "rbcL"= "#D0FF14"))
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

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_12.pdf", width = 10, height = 8)
BCheatmap
dev.off()

###Composition plots by primer pair 
Comp <- ggplot(data=AbPhy, aes(x= Primer_name, y= Relative_abundance, fill= Phyla)) +
        geom_bar(aes(), stat="identity", position="stack") +
        scale_colour_brewer(palette = "Set1") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
        labs(tag = "B)") #+
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ AbyPp$Target_gene) +
        #theme(legend.position="none") + guides(fill=guide_legend(nrow=50)) +

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_11.pdf", width = 12, height = 16)
grid.arrange(BC, Comp, nrow= 2, ncol= 1)
dev.off()  
