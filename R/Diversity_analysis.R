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
             pointsize = log10(foo.primer$Total_reads), 
             geom.ind = "point", # show points only (but not "text")
             col.ind = foo.primer$Total_reads, # color by groups
             gradient.cols = c("blue", "red"),
             alpha.ind = 0.7,
             addEllipses = F, # Concentration ellipses
             legend.title = "Sequencing reads")+
  labs(tag = "B)")

c<- fviz_contrib(foo.pca, choice = "ind", axes = 1:2, top = 10,
             fill = c("pink", "#21908CFF", "#21908CFF", "#440154FF", "pink",
                      "pink", "pink", "pink", "#440154FF", "#440154FF"),
             color = c("pink", "#21908CFF", "#21908CFF", "#440154FF", "pink",
                       "pink", "pink", "pink", "#440154FF", "#440154FF"))+
  labs(tag = "C)")

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
foo.euk2<- foo.euk[,1:48]

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
                 pointsize = log10(foo.euk$Total_reads), 
                 geom.ind = "point", # show points only (but not "text")
                 col.ind = foo.euk$Total_reads, # color by groups
                 gradient.cols = c("blue", "red"),
                 alpha.ind = 0.7,
                 addEllipses = F, # Concentration ellipses
                 legend.title = "Sequencing reads")+
  labs(tag = "B)")

cE<- fviz_contrib(foo.euk2.pca, choice = "ind", axes = 1:2, top = 10,
                 fill = c("#21908CFF", "#21908CFF", "#440154FF", "#440154FF", "#440154FF",
                          "#440154FF", "#440154FF","#440154FF", "#440154FF", "#D0FF14"),
                 color = c("#21908CFF", "#21908CFF", "#440154FF", "#440154FF", "#440154FF",
                           "#440154FF", "#440154FF","#440154FF", "#440154FF", "#D0FF14"))+
  labs(tag = "C)") +
  theme(axis.text.x = element_text(angle=90))

dE<- fviz_pca_var(foo.euk2.pca, col.var = "cos2", select.var = list(contrib = 10),
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
foo.bac2<- foo.bac[,1:48]

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
                  pointsize = log10(foo.bac$Total_reads), 
                  geom.ind = "point", # show points only (but not "text")
                  col.ind = foo.bac$Total_reads, # color by groups
                  gradient.cols = c("blue", "red"),
                  alpha.ind = 0.7,
                  addEllipses = F, # Concentration ellipses
                  legend.title = "Sequencing reads")+
  labs(tag = "B)")

cB<- fviz_contrib(foo.bac2.pca, choice = "ind", axes = 1:2, fill = "pink", color = "pink")+
  labs(tag = "C)")

dB<- fviz_pca_var(foo.bac2.pca, col.var = "cos2", select.var = list(contrib = 15),
                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                  repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

eB<- fviz_contrib(foo.bac2.pca, choice = "var", axes = 1:2, top = 25)+
  labs(tag = "E)")

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_8.pdf", width = 16, height = 8)
grid.arrange(aB, bB, cB, dB, eB, widths = c(1, 1, 1), layout_matrix = rbind(c(1, 3, 5),
                                                                            c(2, 4, 5)))
dev.off()

###Eliminating the highly dominant primers 
foo2<- foo[-c(7,13,18),]
foo.primer2<- foo.primer[-c(7,13,18),]

foo.pca2<- PCA(foo2, graph = T)
foo.eig2<- get_eigenvalue(foo.pca2)
fviz_eig(foo.pca2, addlabels = TRUE, ylim = c(0, 25))

fviz_pca_ind(foo.pca2, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE) # Avoid text overlapping (slow if many points)

###Color them by gene
a2<- fviz_pca_ind(foo.pca2,
                 pointsize = "cos2", ##It is possible to size it by number of reads pointsize = log10(foo.primer$Total_reads)
                 pointshape = 21,
                 geom.ind = "point", # show points only (but not "text")
                 col.ind = foo.primer2$Gen, # color by groups
                 fill.ind =  foo.primer2$Gen,
                 palette = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"),
                 alpha.ind = 0.5,
                 addEllipses = F, # Concentration ellipses
                 legend.title = "Gen")+
  labs(tag = "A)")

###Color them by sequencing reads
b2<- fviz_pca_ind(foo.pca2,
                 pointsize = log10(foo.primer2$Total_reads), 
                 geom.ind = "point", # show points only (but not "text")
                 col.ind = foo.primer2$Total_reads, # color by groups
                 gradient.cols = c("blue", "red"),
                 alpha.ind = 0.7,
                 addEllipses = F, # Concentration ellipses
                 legend.title = "Sequencing reads")+
  labs(tag = "B)")

c2<- fviz_contrib(foo.pca2, choice = "ind", axes = 1:2, top = 10,
                 fill = c("#440154FF", "pink", "pink", "pink", "pink",
                          "#440154FF", "#440154FF", "#440154FF", "#440154FF", "#440154FF"),
                 color = c("#440154FF", "pink", "pink", "pink", "pink",
                           "#440154FF", "#440154FF", "#440154FF", "#440154FF", "#440154FF"))+
  labs(tag = "C)")+
  theme(axis.text.x = element_text(angle=90))

##Variable results
foo.var2 <- get_pca_var(foo.pca2)

###PCA variables
#fviz_pca_var(foo.pca, col.var = "black")
d2<- fviz_pca_var(foo.pca2, col.var = "cos2", select.var =  list(contrib = 10), ##top ten taxa contributing
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                 repel = TRUE )+ # Avoid text overlapping
  labs(tag = "D)")

###Visualize quality of representation 
corrplot(na.omit(foo.var2$cos2), is.corr=FALSE) #Alternative

##Plot contribution
e2<- fviz_contrib(foo.pca2, choice = "var", axes = 1:2, top = 25)+
  labs(tag = "E)")

pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_9.pdf", width = 16, height = 8)
grid.arrange(a2, b2, c2, d2, e2, widths = c(1, 1, 1), heights= c(2,2), layout_matrix = rbind(c(1, 3, 5),
                                                                       c(2, 4, 5)))
dev.off()

###Composition plots by primer pair 
ggplot(data=AbPhy, aes(x= Primer_name, y= Relative_abundance, fill= Phyla)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ AbyPp$Target_gene) +
  theme(legend.position="none") + guides(fill=guide_legend(nrow=50)) 