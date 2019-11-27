##Phyloseq pipeline for diversity analysis by amplicon

library(ggplot2)
library(reshape)
library(phyloseq)
library(data.table)
library(parallel)
library(microbiome)
library(tidyverse)
library(gridExtra)
library(grid)
library(vegan)
library("FactoMineR")
library("factoextra")
library("corrplot")

if(!exists("primerInput")){
  source("~/GitProjects/AA_Primer_Evaluation/R/General_Primer_information.R") ##   
}

if(!exists("PS.l")){
  PS.l<- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/MergedPhyloSeqList.Rds") ##   
}

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

###Rarefied data
PS.l.rar <- lapply(PS.l, function (x) {
  rarefy_even_depth(x, rngseed=1, sample.size=0.3*mean(sample_sums(x)), replace=F)
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
foo<- AbPhy[,1:3]
foo<- reshape(foo, v.names = "ASV", timevar = "Phyla", idvar = "Primer_name",direction = "wide")
colnames(foo) <- gsub("ASV.", "\\1", colnames(foo)) ##Remove "ASV"
##Transform NA to 0
foo[is.na(foo)]<- 0
rownames(foo)<-foo[,1]

##Merge with primer information 
foo.primer<-join(foo, primerInput, by= "Primer_name")
rownames(foo.primer)<-foo.primer[,1]

foo[,1]<- NULL
foo.primer[,1]<- NULL

##Let's do a PCA of individuals
foo.pca<- PCA(foo, graph = T)
foo.eig<- get_eigenvalue(foo.pca)
fviz_eig(foo.pca, addlabels = TRUE, ylim = c(0, 25))

fviz_pca_ind(foo.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE)  # Avoid text overlapping (slow if many points)

fviz_contrib(foo.pca, choice = "ind", axes = 1:2)

###Group them by gene
fviz_pca_ind(foo.pca,
             geom.ind = "point", # show points only (but not "text")
             col.ind = foo.primer$Gen, # color by groups
             palette = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"),
             alpha.ind = 0.5,
             addEllipses = F, # Concentration ellipses
             legend.title = "Gen")


##Variable results
foo.var <- get_pca_var(foo.pca)
fviz_pca_var(foo.pca, col.var = "black")
corrplot(foo.var$cos2, is.corr=FALSE)
fviz_cos2(foo.pca, choice = "var", axes = 1:2)

fviz_pca_var(foo.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE ) # Avoid text overlapping


###Composition plots by primer pair 
ggplot(data=AbPhy, aes(x= Primer_name, y= Relative_abundance, fill= Phyla)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ AbyPp$Target_gene) +
  theme(legend.position="none") + guides(fill=guide_legend(nrow=50)) 