###Pipeline to generate the cummulative curves (Fig. 1)
###Libraries
library(ggplot2)
library(data.table)
library(phyloseq)
library(dplyr)
library(plyr)
library(vegan)
library(gridExtra)
library(grid)
library(lattice)
library("pheatmap")
library("viridisLite")

###Function
unique.taxa.cumsum <- function(ps.numord.l, rank){
  taxonNa <- lapply(ps.numord.l, function (x) 
    get_taxa_unique(x, taxonomic.rank = rank))
  names(taxonNa) <- names(ps.numord.l)
  unique.taxonNa.cumsum <- sapply(0:length(taxonNa),
                                  function (i) length(unique(unlist(taxonNa[0:i]))))
  return(unique.taxonNa.cumsum)
}

###Data
if(!exists("primerInput")){
  source("~/GitProjects/AA_Primer_Evaluation/R/General_Primer_information.R") ##   
}

if(!exists("PS.l")){
  PS.l<- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Primer_Evaluation/MergedPhyloSeqList.Rds") ##   
}

SSU_18S <- F ## Analysis just with 18S primers
SSU_LSU <- F ## Analysis with 18S and 28S primers
SSU_LSU_COI <- F ## Analysis with 18S, 28S and COI primers

####Cummulative curves by taxa with all information
num.taxa <-sapply(c("species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PS.l, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads <- unlist(lapply(PS.l, function (x) 
  sum(otu_table(x))))

PrimTax <- as.data.frame(cbind(num.taxa, num.reads))

PrimTax <- as.data.frame(apply(PrimTax, 2, unlist))
PrimTax[,8] <- row.names(PrimTax)  
colnames(PrimTax)[8] <- "Primer_comb_ID"

###First merge PrimTax with primerInput
PrimTax$Primer_comb_ID <- gsub(pattern = " ", replacement = "", x = PrimTax$Primer_comb_ID) ### check if the primer names have extra spaces
intersect(PrimTax$Primer_comb_ID, primerInput$Primer_comb_ID) 
#setdiff(PrimTax$Primer_name, primerInput$Primer_name)
PrimTax<- join(PrimTax, primerInput, by = "Primer_comb_ID")

ps.numord.l <- PS.l[order(PrimTax$num.reads, decreasing=TRUE)]

cum.tax <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l, rank)
})

colnames(cum.tax) <- paste0("cum.", colnames(cum.tax))

cum.tax <- as.data.frame(cbind(1:nrow(cum.tax), cum.tax))
#prnames <- names(ps.numord.l) ## get names from the primers
#cum.tax[2:40,7] <- prnames ## add primer names 
cum.tax[2:40,7] <- PrimTax$Gen[order(PrimTax$num.reads, decreasing=TRUE)]
colnames(cum.tax)[1] <- "Primer_comb_ID"
colnames(cum.tax)[7] <- "Gen"

a <- ggplot(cum.tax) +
  geom_step(aes(Primer_comb_ID, cum.species), color="red") +
  geom_text(aes(20, 5000, label="Species"), color="red")+ 
  geom_step(aes(Primer_comb_ID, cum.genus), color="#0D0887FF") +
  geom_text(aes(20, 2600, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_comb_ID, cum.family), color="#7E03A8FF") +
  geom_text(aes(20, 1200, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_comb_ID, cum.order), color="#CC4678FF") +
  geom_text(aes(20, 500, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F89441FF") +
  #geom_text(aes(20, 200, label="Classes"), color="#F89441FF") +
  geom_step(aes(Primer_comb_ID, cum.phylum), color="#F89441FF") +
  geom_text(aes(20, 60, label="Phyla"), color="#F89441FF") +
  geom_point(aes(Primer_comb_ID, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"))+
  scale_y_log10("Cummulative count of taxa") +
  #scale_x_continuous("Number of primers considered \n(starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "none", 
        axis.title.x = element_blank(), text = element_text(size=20))+
  labs(tag = "A)")+
  coord_cartesian(ylim = c(20, 5000))

### Parasites 

Primtax.comb2 <- subset(PrimTax, Gen=="28S" | Gen== "18S" | Gen== "COI")

ps.numord.l.comb2 <- PS.l[c(Primtax.comb2$Primer_comb_ID)][order(Primtax.comb2$num.reads, decreasing=TRUE)]

ps.numord.l.comb2[["Euk_COI_04"]]<- NULL ### No parasites or fungi 

PM.para.comb2 <- lapply(ps.numord.l.comb2, function (x) {
  subset_taxa(x, phylum%in%c("Nematoda",
                             "Apicomplexa",
                             "Platyhelminthes"))
})

num.taxa.para <-sapply(c( "species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.para.comb2, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.para <- unlist(lapply(PM.para.comb2, function (x) 
  sum(otu_table(x))))

PrimTax.para <-as.data.frame(cbind(num.taxa.para,as.data.frame(num.reads.para)))
PrimTax.para$Primer_comb_ID<- rownames(PrimTax.para)
markers <- PrimTax[,c(8,17)]
PrimTax.para <- join(PrimTax.para, markers)
PrimTax.para <- PrimTax.para[order(-PrimTax.para$num.reads.para),] 

ps.numord.l.para <- PM.para.comb2[c(PrimTax.para$Primer_comb_ID)][order(as.vector(unlist(PrimTax.para$num.reads.para)), decreasing=TRUE)]

cum.tax.para <- sapply(c("species","genus", "family","order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.para, rank)
})

colnames(cum.tax.para) <- paste0("cum.", colnames(cum.tax.para))

cum.tax.para <- as.data.frame(cbind(1:nrow(cum.tax.para), cum.tax.para))
cum.tax.para[2:28,7] <- PrimTax.para$Gen ## add marker amplified 
colnames(cum.tax.para)[1] <- "Primer_comb_ID"
colnames(cum.tax.para)[7] <- "Gen"

b <- ggplot(cum.tax.para) +
  geom_step(aes(Primer_comb_ID, cum.species), color="red") +
  geom_text(aes(26, 380, label="Species"), color="red") +
  geom_step(aes(Primer_comb_ID, cum.genus), color="#0D0887FF") +
  geom_text(aes(26, 180, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_comb_ID, cum.family), color="#7E03A8FF") +
  geom_text(aes(26, 110, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_comb_ID, cum.order), color="#CC4678FF") +
  geom_text(aes(26, 30, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 10, label="Classes"), color="#F0F921FF") +
  geom_step(aes(Primer_comb_ID, cum.phylum), color="#F89441FF") +
  geom_text(aes(26, 4, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa (Parasites)") +
  #scale_x_continuous("Number of primers considerd \n(starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  geom_point(aes(Primer_comb_ID, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))+
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = 'none',
        axis.title.x = element_blank(), text = element_text(size=20))+
  labs(tag = "B)")

###Fungi
PM.fungi.comb2 <- lapply(ps.numord.l.comb2, function (x) {
  subset_taxa(x, phylum%in%c("Ascomycota",
                             "Basidiomycota",
                             "Zygomycota", "Mucoromycota", "Zoopagomycota", "Blastocladiomycota", "Cryptomycota", "Microsporidia"))
})

####Cummulative curve for fungi
num.taxa.fungi <-sapply(c("species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.fungi.comb2, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.fungi <- unlist(lapply(PM.fungi.comb2, function (x) 
  sum(otu_table(x))))

PrimTax.fungi <-as.data.frame(cbind(num.taxa.fungi,as.data.frame(num.reads.fungi)))
PrimTax.fungi$Primer_comb_ID<- rownames(PrimTax.fungi)
markers <- PrimTax[,c(8,17)]
PrimTax.fungi <- join(PrimTax.fungi, markers)
PrimTax.fungi <- PrimTax.fungi[order(-PrimTax.fungi$num.reads.fungi),] 

ps.numord.l.fungi <- PM.fungi.comb2[c(PrimTax.fungi$Primer_comb_ID)][order(as.vector(unlist(PrimTax.fungi$num.reads.fungi)), decreasing=TRUE)]

cum.tax.fungi <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.fungi, rank)
})

colnames(cum.tax.fungi) <- paste0("cum.", colnames(cum.tax.fungi))

cum.tax.fungi <- as.data.frame(cbind(1:nrow(cum.tax.fungi), cum.tax.fungi))
cum.tax.fungi[2:28,7] <- PrimTax.fungi$Gen ## add marker amplified 
colnames(cum.tax.fungi)[1] <- "Primer_comb_ID"
colnames(cum.tax.fungi)[7] <- "Gen"

c <- ggplot(cum.tax.fungi) +
  geom_step(aes(Primer_comb_ID, cum.species), color="red") +
  geom_text(aes(26, 2000, label="Species"), color="red") +
  geom_point(aes(Primer_comb_ID, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("#440154FF", "#21908CFF", "#FDE725FF"))+
  geom_step(aes(Primer_comb_ID, cum.genus), color="#0D0887FF") +
  geom_text(aes(26, 1000, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_comb_ID, cum.family), color="#7E03A8FF") +
  geom_text(aes(26, 440, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_comb_ID, cum.order), color="#CC4678FF") +
  geom_text(aes(26, 160, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 50, label="Classes"), color="") +
  geom_step(aes(Primer_comb_ID, cum.phylum), color="#F89441FF") +
  geom_text(aes(26, 8, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa (Fungi)") +
  #scale_x_continuous("Number of primers considerd \n(starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = 'none', 
        axis.title.x = element_blank(), text = element_text(size=20))+
  labs(tag = "C)")

###Other eukaryotic elements ("Diet" and passing material)
##Working approach 
PM.diet.comb3 <- lapply(ps.numord.l, function (x) {
  subset_taxa(x, phylum%in%c("Arthropoda", "Chlorophyta", "Streptophyta", "Chordata", 
                             "Annelida", "Mollusca", "Bryozoa", "Picozoa", "Porifera",
                             "Rotifera", "Tardigrada"))
})

num.taxa.diet <-sapply(c("species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.diet.comb3, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.diet <- unlist(lapply(PM.diet.comb3, function (x) 
  sum(otu_table(x))))

PrimTax.diet <-as.data.frame(cbind(num.taxa.diet,as.data.frame(num.reads.diet)))
PrimTax.diet$Primer_comb_ID<- rownames(PrimTax.diet)
markers <- PrimTax[,c(8,17)]
PrimTax.diet <- join(PrimTax.diet, markers)
PrimTax.diet <- PrimTax.diet[order(-PrimTax.diet$num.reads.diet),] 

ps.numord.l.diet <- PM.diet.comb3[c(PrimTax.diet$Primer_comb_ID)][order(as.vector(unlist(PrimTax.diet$num.reads.diet)), decreasing=TRUE)]

cum.tax.diet <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.diet, rank)
})

colnames(cum.tax.diet) <- paste0("cum.", colnames(cum.tax.diet))

cum.tax.diet <- as.data.frame(cbind(1:nrow(cum.tax.diet), cum.tax.diet))
cum.tax.diet[2:40,7] <- PrimTax.diet$Gen ## add marker amplified 
colnames(cum.tax.diet)[1] <- "Primer_comb_ID"
colnames(cum.tax.diet)[7] <- "Gen"

d<- ggplot(cum.tax.diet) +
  geom_step(aes(Primer_comb_ID, cum.species), color="red") +
  geom_text(aes(35, 1200, label="Species"), color="red") +
  geom_step(aes(Primer_comb_ID, cum.genus), color="#0D0887FF") +
  geom_text(aes(35, 700, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_comb_ID, cum.family), color="#7E03A8FF") +
  geom_text(aes(35, 400, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_comb_ID, cum.order), color="#CC4678FF") +
  geom_text(aes(35, 160, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 50, label="Classes"), color="") +
  geom_step(aes(Primer_comb_ID, cum.phylum), color="#F89441FF") +
  geom_text(aes(35, 15, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa \n (non fungi or parasite)") +
  geom_point(aes(Primer_comb_ID, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("#E3DAC9","pink","#440154FF", "#21908CFF", "#FDE725FF", "#C46210", "#D0FF14"))+
  #scale_x_continuous("Number of primers considerd \n(starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = 'none', 
        axis.title.x = element_blank(), text = element_text(size=20))+
  labs(tag = "D)")

###Compile and save all the figure

ab <- grid.arrange(a,b, nrow= 1, ncol= 2,  
             bottom= textGrob("Number of primers considerd (starting with the one with highest read count)", 
                              gp= gpar(fontsize= 20)))

cd <- grid.arrange(c,d, nrow= 1, ncol= 2,  
                   bottom= textGrob("Number of primers considerd (starting with the one with highest read count)", 
                                    gp= gpar(fontsize= 20)))

#pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_1.pdf", width = 12, height = 15)
#grid.arrange(ab, cd, nrow= 2, ncol= 1)
#dev.off()

####Just 18S
if(SSU_18S){
  Primtax.comb18 <- subset(PrimTax, Gen== "18S")
  
  ps.numord.l.comb18 <- PS.l[c(Primtax.comb18$Primer_name)][order(Primtax.comb18$num.reads, decreasing=TRUE)]
  
  
  cum.tax.comb18 <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
    unique.taxa.cumsum(ps.numord.l.comb18, rank)
  })
  
  colnames(cum.tax.comb18) <- paste0("cum.", colnames(cum.tax.comb18))
  
  cum.tax.comb18 <- as.data.frame(cbind(1:nrow(cum.tax.comb18), cum.tax.comb18))
  #prnames28s <- names(ps.numord.l.28S) ## get names from the primers
  #cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
  colnames(cum.tax.comb18)[1] <- "Primer_name"
  
  cum.plot.comb18 <- ggplot(cum.tax.comb18) +
    geom_step(aes(Primer_name, cum.species), color="red") +
    geom_text(aes(20, 2200, label="Species"), color="red")+ 
    geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
    geom_text(aes(20, 1600, label="Genera"), color="#0D0887FF") +
    geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
    geom_text(aes(20, 840, label="Families"), color="#7E03A8FF") +
    geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
    geom_text(aes(20, 380, label="Orders"), color="#CC4678FF") +
    #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
    #geom_text(aes(20, 135, label="Classes"), color="#F0F921FF") +
    geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
    geom_text(aes(20, 40, label="Phyla"), color="#F89441FF") +
    scale_y_log10("Cummulative count of taxa") +
    #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
    annotation_logticks(sides="l") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), text = element_text(size=20))+
    labs(tag = "B)")+
    coord_cartesian(ylim = c(20, 4000))
}
####18S and 28S
if(SSU_LSU){
  Primtax.comb28 <- subset(PrimTax, Gen=="28S" | Gen== "18S")
  
  ps.numord.l.comb28 <- PS.l[c(Primtax.comb28$Primer_name)][order(Primtax.comb28$num.reads, decreasing=TRUE)]
  
  cum.tax.comb28 <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
    unique.taxa.cumsum(ps.numord.l.comb28, rank)
  })
  
  colnames(cum.tax.comb28) <- paste0("cum.", colnames(cum.tax.comb28))
  
  cum.tax.comb28 <- as.data.frame(cbind(1:nrow(cum.tax.comb28), cum.tax.comb28))
  #prnames28s <- names(ps.numord.l.28S) ## get names from the primers
  #cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
  colnames(cum.tax.comb28)[1] <- "Primer_name"
  
  cum.plot.comb28 <- ggplot(cum.tax.comb28) +
    geom_step(aes(Primer_name, cum.species), color="red") +
    geom_text(aes(22, 3300, label="Species"), color="red")+ 
    geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
    geom_text(aes(22, 2000, label="Genera"), color="#0D0887FF") +
    geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
    geom_text(aes(22, 1000, label="Families"), color="#7E03A8FF") +
    geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
    geom_text(aes(22, 420, label="Orders"), color="#CC4678FF") +
    #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
    #geom_text(aes(20, 135, label="Classes"), color="#F0F921FF") +
    geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
    geom_text(aes(22, 42, label="Phyla"), color="#F89441FF") +
    scale_y_log10("Cummulative count of taxa") +
    #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
    annotation_logticks(sides="l") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), text = element_text(size=20))+
    labs(tag = "C)")+
    coord_cartesian(ylim = c(20, 4000))
  
}
####28S, 18S, COI####
if(SSU_LSU_COI){
  Primtax.comb2 <- subset(PrimTax, Gen=="28S" | Gen== "18S" | Gen== "COI")
  
  ps.numord.l.comb2 <- PS.l[c(Primtax.comb2$Primer_name)][order(Primtax.comb2$num.reads, decreasing=TRUE)]
  
  
  cum.tax.comb2 <- sapply(c("species","genus", "family", "order", "phylum"), function (rank){
    unique.taxa.cumsum(ps.numord.l.comb2, rank)
  })
  
  colnames(cum.tax.comb2) <- paste0("cum.", colnames(cum.tax.comb2))
  
  
  cum.tax.comb2 <- as.data.frame(cbind(1:nrow(cum.tax.comb2), cum.tax.comb2))
  #prnames28s <- names(ps.numord.l.28S) ## get names from the primers
  #cum.tax.18S[2:22,1] <- prnames18s ## add primer names 
  colnames(cum.tax.comb2)[1] <- "Primer_name"
  
  cum.plot.comb2 <- ggplot(cum.tax.comb2) +
    geom_step(aes(Primer_name, cum.species), color="red") +
    geom_text(aes(26, 3800, label="Species"), color="red")+ 
    geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
    geom_text(aes(26, 2200, label="Genera"), color="#0D0887FF") +
    geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
    geom_text(aes(26, 1100, label="Families"), color="#7E03A8FF") +
    geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
    geom_text(aes(26, 450, label="Orders"), color="#CC4678FF") +
    #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
    #geom_text(aes(20, 135, label="Classes"), color="#F0F921FF") +
    geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
    geom_text(aes(26, 43, label="Phyla"), color="#F89441FF") +
    scale_y_log10("Cummulative count of taxa") +
    #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
    annotation_logticks(sides="l") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), text = element_text(size=20))+
    labs(tag = "D)")+
    coord_cartesian(ylim = c(20, 4000))
  
  pdf(file = "~/AA_Primer_evaluation/Figures/Figure_1.pdf", width = 10, height = 8)
  grid.arrange(cum.plot.all, cum.plot.comb18, cum.plot.comb28, cum.plot.comb2, nrow= 2, ncol= 2,
               bottom= textGrob("Number of primers considerd (starting with the one with highest read count)", 
                                gp= gpar(fontsize= 20)))
  dev.off()
  
}

rm(a,b,ab,c,cd,d, 
   cum.tax, cum.tax.diet, cum.tax.fungi, cum.tax.para,
   num.taxa, num.taxa.diet, num.taxa.fungi, num.taxa.para,
   PrimTax, Primtax.comb2, PrimTax.diet, PrimTax.fungi, PrimTax.para,
   num.reads, num.reads.diet, num.reads.fungi, num.reads.para,
   PM.diet.comb3, PM.fungi.comb2, PM.para.comb2,
   ps.numord.l, ps.numord.l.comb2, ps.numord.l.diet, ps.numord.l.fungi, ps.numord.l.para, markers)