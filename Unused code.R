##Unused code 

####Emanuel's function (don't work)
#sumSeqByTax2 <- function (Phy, tax) {
#  counts <- data.frame(cbind(asvCount=colSums(otu_table(Phy)), tax_table(Phy)))
#  counts$asvCount <- as.numeric(as.character(counts$asvCount))
#  tapply(counts$asvCount, counts[, tax], sum)
#}
#function for calculating relative amp efficiency for each taxon (Wrong! it can't be calculated for our results, at least not in this way)
##Mod from Kelly et al. 2019
EFFICIENCY <- function(Relative_abundance, Total_ASV, numberCycles, ...){
  temp <- unlist(
    (Total_ASV/Relative_abundance)^(1/numberCycles) - 1)
  temp[which(is.infinite(temp))] <- NA
  temp <- temp/(max(temp, na.rm=T))  #normalize to maximum, to create relative index of amp efficiency
  temp[temp<=0] <- 0
  as.numeric(temp)
}
##old pipeline
#MAF <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Fox/MAF_complete.RDS") ##Without taxa data
###No chimeras data 
#STNCF <- getSequenceTableNoChime(MAF)
###Annotation list 
#annot.list.f <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Fox/Fox_blast_tax_complete.Rds") 
##Discard primers without annotation 
#keepf <- unlist(lapply(annot.list.f, nrow))>0
#names(STNCF)[keepf]<-gsub(pattern = "-", replacement = "_", x= names(STNCF[keepf]))
#names(STNCF)[keepf]<-gsub(pattern = " ", replacement = "", x= names(STNCF[keepf]))
#names(PS.l.f) <- names(STNCF)[keepf]

##Old pipeline
#MAR <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_Raccoon/MAR_complete.RDS") 
###No chimeras data 
#STNCR <- getSequenceTableNoChime(MAR)
###Annotation list 
#annot.list.r <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_Raccoon/Racc_blast_tax_all.Rds")
##Discard primers without annotation 
#keepr <- unlist(lapply(annot.list.r, nrow))>0
#names(STNCR)[keepr]<-gsub(pattern = "-", replacement = "_", x= names(STNCR[keepr]))
#names(STNCR)[keepr]<-gsub(pattern = " ", replacement = "", x= names(STNCR[keepr]))

###Merge both phyloseq lists ... Use as reference Raccoon list 
#(in theory the same primer pairs worked but "515F_99.R1073_103" and "16SA_L_110.16SB_H_110" didn't work for raccoons)

##Let's eliminate them to merge them correclty
##if both lists have the same primer should work properly 
#PS.l.f[["16SA_L_110_F.16SB_H_110_R"]]<- NULL
#PS.l.f[["515F_99_F.R1073_103_R"]]<- NULL

#readNumByPhylumF<- lapply(getTaxonTable(MAF3), function (x) table(as.vector(x[, "phylum"])))
#readNumByPhylumR<- lapply(getTaxonTable(MAR2), function (x) table(as.vector(x[, "phylum"])))
# Then it is necessary to eliminate all empty primer pairs to make the rest of the code work (annoying)

#Function is not working anymore :( find out why
#readNumByGenus <- lapply(PS.l, sumSeqByTax, "genus") 
#readNumByFam <- lapply(PS.l, sumSeqByTax, "family")
#readNumbySpecies <-lapply(PS.l, sumSeqByTax, "species")

###richeness analysis 
###Calculate alpha diversity indexes for every primer pair 
#richness <- data.frame()

#for (i in 1:length(PS.l)){

#  a <- data.frame()

#  b <- data.frame(estimate_richness(PS.l[[i]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson")))

#  a<-colMeans(b)

#  richness<- rbind(richness, a)
#}

#richness[,7] <- names(PS.l)

#colnames(richness) <- c("Observed","Chao1","se_Chao1", "Shannon", "Simpson", "InvSimpson", "Primer_name")


### Obtain the number of samples 
#nsamp <- data.frame()

#nsamp[1:41,1] <-names(PS.l)

#for (i in 1:length(PS.l))
#  nsamp[i,2]<- nrow(otu_table(PS.l.f[[i]])) + nrow(otu_table(PS.l.r[[i]]))

#colnames(nsamp) <- c("Primer_name", "Number_samples")


###Check both elements and merge! 
#richness$Primer_name <- gsub(pattern = " ", replacement = "", x = richness$Primer_name) ### check if the primer names have extra spaces
#nsamp$Primer_name <- gsub(pattern = " ", replacement = "", x= nsamp$Primer_name)
#rawcounts <- as.data.frame(num.reads)
#rawcounts[,2] <- rownames(rawcounts)
#colnames(rawcounts)[2] <- "Primer_name"


#primers$Primer_name <- gsub(pattern = " ", replacement = "", x= primers$Primer_name)

#setdiff(richness$Primer_name, nsamp$Primer_name)

#richness <- merge(richness, nsamp, by= "Primer_name")

#richness <- richness %>% distinct(Primer_name, .keep_all = TRUE)

#richness$Primer_name <- gsub(pattern = "-", replacement = "_", x = richness$Primer_name)

#richness <- merge(richness, rawcounts, "Primer_name", all= T) %>% distinct(Primer_name, .keep_all = TRUE)

#richness <- merge(richness, primerInput, "Primer_name")

#richness <- richness[which(richness$num.reads>8),] 

#Div <- richness$Shannon
#Rich <- log(richness$Observed + 1)
#evenness = Div/Rich
#richness$Evenness = evenness

###Plots

#ob <- ggplot(data = richness, aes(x = Gen, y = na.omit(Observed), color= Target)) +
#  geom_point() +
#facet_grid(.~Gen) +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Observed richness")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#ggplot(data = richness, aes(x = Gen, y = na.omit(Observed), fill= Target)) +
#  geom_boxplot() +
#  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
#facet_grid(.~Gen) +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Observed richness")+
#scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#ch1<- ggplot(data = richness, aes(x = Gen, y = na.omit(Chao1), color= Target)) +
#  geom_point() +
#facet_grid(.~Gen) +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Chao-1 richness")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#sh <- ggplot(data = richness, aes(x = Gen, y = na.omit(Shannon), color= Target)) +
#  geom_point() +
#facet_grid(.~Gen) +
#  theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Diversity (Shannon Index)")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#sim<- ggplot(data = richness, aes(x = Gen, y = na.omit(Simpson), color= Target)) +
#  geom_point() +
#facet_grid(.~Gen) +
# theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Diversity (Simpson Index)")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#isim<-ggplot(data = richness, aes(x = Gen, y = na.omit(InvSimpson), color= Target)) +
#  geom_point() +
#facet_grid(.~Gen) +
# theme_classic() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#  scale_y_continuous(name = "Diversity (Inverse Simpson Index)")+
#  scale_colour_brewer(palette = "Set1") +
#  theme(legend.position="right") 

#eve<- ggplot(data = richness, aes(x = Gen, y = na.omit(Evenness), color= Target)) +
#  geom_point() +
#facet_grid(.~Gen) +
# theme_classic() +
# theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
# scale_y_continuous(name = "Evenness (Shannon index/ log(Observed richness +1))")+
# scale_colour_brewer(palette = "Set1") +
# theme(legend.position="right") 


#grid.arrange(ob, ch1, sh, sim, isim, eve, nrow= 2, ncol= 3)

#linob <- ggplot(data = richness, aes(x = log10(num.reads), y = Observed, color= Target)) +
#       geom_point() +
#facet_grid(.~Gen) +
#       theme_classic() +
#       geom_smooth(method = "lm", col = "black") +
#theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#        scale_x_continuous(name = "log10 (Sequencing reads)")+
#       scale_colour_brewer(palette = "Set1") +
#       theme(legend.position="right") 

#linshan <- ggplot(data = richness, aes(x = log10(num.reads), y = Shannon, color= Target)) +
#         geom_point() +
#facet_grid(.~Gen) +
#         theme_classic() +
#         geom_smooth(method = "lm", col = "black") +
#theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#         scale_x_continuous(name = "log10 (Sequencing reads)")+
#         scale_colour_brewer(palette = "Set1") +
#         theme(legend.position="right") 


#lineve <-  ggplot(data = richness, aes(x = log10(num.reads), y = Evenness, color= Target)) +
#          geom_point() +
#facet_grid(.~Gen) +
#          theme_classic() +
#          geom_smooth(method = "lm", col = "black") +
#theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#         scale_x_continuous(name = "log10 (Sequencing reads)")+
#          scale_colour_brewer(palette = "Set1") +
#          theme(legend.position="right") 

#summary(lm(log10(richness$Raw_counts)~ richness$Observed))$adj.r.squared

###Plot raw counts by primer pair 
#ggplot(richness, aes(x= reorder(Primer_name, -num.reads), y=num.reads, color= Gen, fill= Gen)) + 
#  geom_bar(stat = "identity")+
#  scale_y_log10(name = "log10 Total sequencing reads")+
#  theme_minimal() +
#  labs(x= "Primer name")+
#  geom_hline(yintercept=100, linetype="dashed", color = "black") +
#  theme(axis.text.x = element_text(angle = 90, vjust = 1))

####Estimate efficiency (function from Kelly et al. 2019) ###This is not possible with our data :(

#AbPhy$Efficieny <- EFFICIENCY(Relative_abundance = AbPhy$Relative_abundance, Total_ASV = AbPhy$Total_ASV, numberCycles = AbPhy$numberCycles)
#AbPhy$New_efficieny <- EFFICIENCY(Relative_abundance = AbPhy$Relative_abundance, Total_ASV = AbPhy$Total_reads, numberCycles = AbPhy$numberCycles)
#AbPhyF$Efficieny <- EFFICIENCY(Relative_abundance = AbPhyF$Relative_abundance, Total_ASV = AbPhyF$Total_ASV, numberCycles = AbPhyF$numberCycles)
#AbPhyR$Efficieny <- EFFICIENCY(Relative_abundance = AbPhyR$Relative_abundance, Total_ASV = AbPhyR$Total_ASV, numberCycles = AbPhyR$numberCycles)
#AbGen$Efficiency <- EFFICIENCY(Relative_abundance = AbGen$Relative_abundance, Total_ASV = AbGen$Total_ASV, numberCycles = AbGen$numberCycles)
#AbGen$numberCycles <- rep(35,nrow(AbGen))

###Plot distribution of efficiency by primer (By Phylum)

#require("ggridges")
#pdf(file = "~/AA_Primer_evaluation/Figures/Amplification_efficiency_FoxvsRaccoon.pdf", width = 10, height = 20)
#a <- ggplot(AbPhy, aes(x = Efficieny, y = reorder(Primer_name, Total_reads))) +
#      geom_density_ridges(aes(fill = Gen, alpha=0.5))+
#      geom_jitter(aes(alpha= 0.5))+
#      scale_x_continuous(limits = c(0,1))+
#      geom_vline(aes(xintercept=0.5),
#             linetype="dashed", size=1, colour="red")+
#      xlab("Amplification efficiency")+
#      ylab("Primer name")+
#      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#      labs(tag= "A)")

#b <- ggplot(AbPhyF, aes(x = Efficieny, y = reorder(Primer_name, Raw_counts))) +
#      geom_density_ridges(aes(fill = Gen, alpha=0.5))+
#      geom_jitter(aes(alpha= 0.5))+
#      scale_x_continuous(limits = c(0,1))+
#      geom_vline(aes(xintercept=0.5),
#             linetype="dashed", size=1, colour="red")+
#      xlab("Amplification efficiency")+
#      ylab("Primer name")+
#      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#      labs(tag = "B)")

#c <- ggplot(AbPhyR, aes(x = Efficieny, y = reorder(Primer_name, Raw_counts))) +
#      geom_density_ridges(aes(fill = Gen, alpha=0.5))+
#      geom_jitter(aes(alpha= 0.5))+
#      scale_x_continuous(limits = c(0,1))+
#      geom_vline(aes(xintercept=0.5),
#             linetype="dashed", size=1, colour="red")+
#      xlab("Amplification efficiency")+
#      ylab("Primer name")+
#      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#      labs(tag = "C)")

#grid.arrange(a, b, c, ncol= 1, nrow= 3)
#dev.off()

#ggplot(AbPhy, aes(x = New_efficieny, y = reorder(Primer_name, Total_reads))) +
#  geom_density_ridges(aes(fill = Gen, alpha=0.5))+
#geom_jitter(aes(alpha= 0.5))+
#  scale_x_continuous(limits = c(0,1))+
#  geom_vline(aes(xintercept=.5),
#             linetype="dashed", size=1, colour="red")+
#  xlab("Amplification efficiency")+
#  ylab("Primer name")+
#  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#  labs(tag= "A)")

#d <- ggplot(AbGen, aes(x = Efficiency, y = reorder(Primer_name, Total_reads))) +
#      geom_density_ridges(aes(fill = Gen, alpha=0.5))+
#      geom_jitter(aes(alpha= 0.005))+
#      scale_x_continuous(limits = c(0,1))+
#      geom_vline(aes(xintercept=0.5),
#             linetype="dashed", size=1, colour="red")+
#      xlab("Amplification efficiency")+
#      ylab("Primer name")+
#      theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
#      labs(tag= "B)") #+ annotate("density", x = 1, y = "L2513_64_F.H2714_64_R", label = "OK")

###The lower efficiency, the higher amplification probability 

#to fit to a Beta, make strictly non-zero, non-one
#unZeroOne <- function(x){
#  x[x<1e-5] <- 1e-5  #set zero values to near zero
#  x[x> (1 - 1e-5)] <- (1 - 1e-5) #set one values to near one
#  x
#}

#######cool function to analyse diversity by primer pair########

#names(PS.l)<- names(STNC[keep])


#setdiff(names(PS.l), as.vector(primerInput$Primer_name))
#setdiff(names(PS.l.f), as.vector(primerInput$Primer_name))
#setdiff(names(PS.l.r), as.vector(primerInput$Primer_name))

#Index<- "Shannon" ###Change according to the index is required 
#vector_of_primers<- as.vector(primerInput$Primer_name)

#Pdiversity <- function(vector_of_primers, PS.l, Index)
#{
#  finalPrime <- list(list())
#  for(i in vector_of_primers)
#  {
#    if(i %in% names(PS.l))
#    {
#      id <- which(i == names(PS.l))
#      finalPrime[[i]][["richness_index"]] <- data.frame(estimate_richness(PS.l[[id]], measures = c("Observed","Chao1", "Shannon", "Simpson", "InvSimpson")))
#      Div <- finalPrime[[i]][["richness_index"]]$Shannon
#      Rich <- log(finalPrime[[i]][["richness_index"]]$Observed + 1)
#      evenness = Div/Rich
#      finalPrime[[i]][["richness_index"]]$Evenness = evenness
#      finalPrime[[i]][["number_of_samples"]] <- nrow(otu_table(PS.l[[id]])) ##Not correct 
#      mean_ind <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = mean)
#      std_dev <- apply(finalPrime[[i]][["richness_index"]], 2, FUN = sd)
#      ci <- qnorm(0.95)*std_dev/sqrt(finalPrime[[i]][["number_of_samples"]])
#      upper <- mean_ind+ci
#      lower <- mean_ind-ci
#      finalPrime[[i]][["central_tendency"]] <- data.frame(Index_name = c("Observed","Chao1", "se.chao1", "Shannon", "Simpson", "InvSimpson", "Evenness"), 
#                                                          Mean = mean_ind, 
#                                                          SD = std_dev,
#                                                          Upper = upper, 
#                                                          Lower = lower)
#     finalPrime[[i]][["primer_info"]] <- data.frame(primerInput[ primerInput$Primer_name %in% i, c(6:ncol(primerInput))])
#    }
#  }
#  finalPrime[[1]] <- NULL
# plots
#  plot_df <- data.frame()
#  x <- data.frame()
#  for(j in names(finalPrime))
#  {
#    x <- finalPrime[[j]][["richness_index"]]
#    x[1:nrow(x),8] <- j
#    x[1:nrow(x),9] <- finalPrime[[j]][["primer_info"]]$Gen
#    x[1:nrow(x),10] <- finalPrime[[j]][["primer_info"]]$Target
#    plot_df <- rbind(plot_df, x)
#  }
#  colnames(plot_df)[8] <- c("primer_name")
#  colnames(plot_df)[9] <- "gene"
#  colnames(plot_df)[10] <- "target"

#  ggplot(plot_df, aes(x= primer_name, y = plot_df[,Index], fill = gene))+
#    geom_boxplot(outlier.colour="black", alpha = 0.5)+
#    geom_jitter(shape=1)+
#    xlab("Primer name")+
#    ylab(paste0(Index, " Diversity Index", collapse = ''))+ ###Change according to the index that is ploted 
#    theme_classic()+
#    theme(axis.text.x = element_text(angle = 90, vjust = 1))


# upsetr
#  require(UpSetR)
# We need the list of sample names for a given vector of primers
#  upset_list <- list()
#  for(k in vector_of_primers)
#  {
#    upset_list[[k]] <- rownames(finalPrime[[k]][["richness_index"]])
#  }
#  upset(fromList(upset_list), order.by = "freq", mainbar.y.label = "Number of samples", sets.x.label = "Sample size", nintersects = NA, nsets = 41)

#}

#####Parasites All primers (Dodgy solution, but working!!!!)
Primtax.comb <- subset(PrimTax, Gen=="28S" | Gen== "18S" | Gen== "COI" | Gen=="16S")
ps.numord.l.comb <- PS.l[c(Primtax.comb$Primer_name)][order(Primtax.comb$num.reads, decreasing=TRUE)]

ps.numord.l.comb[["L2513_64_F.H2714_64_R"]] <- NULL
ps.numord.l.comb[["Ins16S_1shortF_70_F.Ins16S_1shortR_70_R"]] <- NULL
ps.numord.l.comb[["Coleop_16Sc_68_F.Coleop_16Sd_68_R"]] <- NULL
ps.numord.l.comb[["F1073_101_F.1492R_100_R"]] <- NULL
ps.numord.l.comb[["ADM330_18_F.Klin0785_CR_19_R"]] <- NULL
ps.numord.l.comb[["16S.1100.F16_100_F.1492R_100_R"]] <- NULL
ps.numord.l.comb[["ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R"]]<- NULL

PM.para.comb <- lapply(ps.numord.l.comb, function (x) {
  subset_taxa(x, phylum%in%c("Nematoda",
                             "Apicomplexa",
                             "Platyhelminthes"))
})


PG.para.comb <- lapply(PM.para.comb, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})


mat.para.comb <- lapply(PG.para.comb, function (x) {
  as.matrix(tax_table(x))
}) 

P.comb <- lapply(PG.para.comb, function (y) {
  make.names(tax_table(y)[, 5])
}) 

P.comb$L2513_64_F.H2714_64_R<- list()
P.comb$Ins16S_1shortF_70_F.Ins16S_1shortR_70_R<- list()
P.comb$Coleop_16Sc_68_F.Coleop_16Sd_68_R<- list()
P.comb$Ave12F_91_F.Ave12R_92_R<- list()
P.comb$F1073_101_F.1492R_100_R<- list()
P.comb$ADM330_18_F.Klin0785_CR_19_R<- list()
P.comb$`16S.1100.F16_100_F.1492R_100_R`<- list()
P.comb$S2F_80_F.S3R_81_R<- list()
P.comb$F52_69_F.R193_69_R<- list()
P.comb$ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R<- list()

####compair all the primers for multiple intersections

P.intersection <- sapply(seq_len(length(P.comb)), function(x)
  sapply(seq_len(length(P.comb)), function(y) length(intersect(unlist(P.comb[x]), unlist(P.comb[y]))))
)

rownames(P.intersection)<- names(P.comb)
colnames(P.intersection)<- names(P.comb)

P.clust <- hclust(dist(P.intersection), method = "complete") ##Dendogram
#library(dendextend)
as.dendrogram(P.clust) %>%
  plot(horiz = TRUE)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))
P.col$Primer_name <- rownames(P.col)

P.col <- merge(P.col, PrimTax, by="Primer_name", sort= F)

col_groups <- P.col %>%
  select("Primer_name", "Gen", "Region") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Primer_name

col_groups$Primer_name<- NULL

colour_groups <- list( Gen= c("18S"= "#440154FF", 
                              "28S"= "#21908CFF", 
                              "COI"= "#FDE725FF", 
                              "16S"= "pink",
                              "12S"= "#E3DAC9",
                              "ITS"= "#C46210",
                              "rbcL" = "#D0FF14"))

paraheatmap <- pheatmap(P.intersection, 
                        color = plasma(100),
                        border_color = NA,
                        annotation_col = col_groups, 
                        #annotation_row = col_groups,
                        annotation_colors = colour_groups,
                        #cutree_rows = 2,
                        #cutree_cols = 2,
                        show_rownames = F,
                        show_colnames = F,
                        main= "Redundant genus of parasites amplified")

ggsave("~/AA_Primer_evaluation/Figures/Parasites_redundancy_allprimers.pdf", plot = paraheatmap, dpi = 450, width = 14, height = 10)

####Cummulative curve for parasites (All primers)
num.taxa.para <-sapply(c( "species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.para.comb, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.para <- unlist(lapply(PM.para.comb, function (x) 
  sum(otu_table(x))))

PrimTax.para <-as.data.frame(cbind(num.taxa.para,as.data.frame(num.reads.para)))
PrimTax.para$Primer_name<- rownames(PrimTax.para)
markers <- PrimTax[,c(8,16)]
PrimTax.para <- join(PrimTax.para, markers)
PrimTax.para <- PrimTax.para[order(-PrimTax.para$num.reads.para),] 

ps.numord.l.para <- PM.para.comb[c(PrimTax.para$Primer_name)][order(as.vector(unlist(PrimTax.para$num.reads.para)), decreasing=TRUE)]

cum.tax.para <- sapply(c("species","genus", "family","order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.para, rank)
})

colnames(cum.tax.para) <- paste0("cum.", colnames(cum.tax.para))

cum.tax.para <- as.data.frame(cbind(1:nrow(cum.tax.para), cum.tax.para))
cum.tax.para[2:30,7] <- PrimTax.para$Gen ## add marker amplified 
colnames(cum.tax.para)[1] <- "Primer_name"
colnames(cum.tax.para)[7] <- "Gen"
cum.tax.para<- add_row(cum.tax.para, Primer_name= 31:40)
cum.tax.para[31:40,2:6]<-cum.tax.para[30,2:6] 
cum.tax.para[31:40,7]<-col_groups[30:39,1] 

cum.plot.para <- ggplot(cum.tax.para) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(26, 380, label="Species"), color="red") +
  geom_point(aes(Primer_name, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("#E3DAC9","pink","#440154FF","#21908CFF","#FDE725FF","#C46210","#D0FF14"))+
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(26, 180, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(26, 110, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(26, 30, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 10, label="Classes"), color="#F0F921FF") +
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(26, 4, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa (Parasites)") +
  #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = 'none', axis.title.x = element_blank())+
  labs(tag = "B)")#+
#coord_cartesian(ylim = c(20, 3800))
#ggsave("~/AA_Primer_evaluation/Figures/Parasites_cummulative.pdf", plot = cum.plot.para, dpi = 450, width = 14, height = 10)

#####Fungi with all primers 

ps.numord.l.comb[["Klin0341_19_F.Klin0785_CR_19_R"]]<- NULL
ps.numord.l.comb[["27M_F_98_F.CCR_98_R"]]<- NULL

PM.fungi.comb <- lapply(ps.numord.l.comb, function (x) {
  subset_taxa(x, phylum%in%c("Ascomycota",
                             "Basidiomycota",
                             "Zygomycota", "Mucoromycota", "Zoopagomycota", "Blastocladiomycota", "Cryptomycota", "Microsporidia"))
})

PG.fungi.comb <- lapply(PM.fungi.comb, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})

mat.fungi.comb <- lapply(PG.fungi.comb, function (x) {
  as.matrix(tax_table(x))
}) 


F.comb <- lapply(PG.fungi.comb, function (y) {
  make.names(tax_table(y)[,5])
})

###Add empty primers
F.comb$L2513_64_F.H2714_64_R<- list()
F.comb$Ins16S_1shortF_70_F.Ins16S_1shortR_70_R<- list()
F.comb$Coleop_16Sc_68_F.Coleop_16Sd_68_R<- list()
F.comb$Ave12F_91_F.Ave12R_92_R<- list()
F.comb$F1073_101_F.1492R_100_R<- list()
F.comb$ADM330_18_F.Klin0785_CR_19_R<- list()
F.comb$`16S.1100.F16_100_F.1492R_100_R`<- list()
F.comb$S2F_80_F.S3R_81_R<- list()
F.comb$F52_69_F.R193_69_R<- list()
F.comb$ZBJ_ArtF1c_67_F.ZBJ_ArtR2c_67_R<- list()
F.comb$Klin0341_19_F.Klin0785_CR_19_R<- list()
F.comb$`27M_F_98_F.CCR_98_R`<- list()

####compair all the primers for multiple intersections

F.intersection <- sapply(seq_len(length(F.comb)), function(x)
  sapply(seq_len(length(F.comb)), function(y) length(intersect(unlist(F.comb[x]), unlist(F.comb[y]))))
)

rownames(F.intersection)<- names(F.comb)
colnames(F.intersection)<- names(F.comb)

F.clust <- hclust(dist(F.intersection), method = "complete") ##Dendogram
#library(dendextend)
as.dendrogram(F.clust) %>%
  plot(horiz = TRUE)

F.col <- cutree(tree = F.clust, k = 2)
F.col  <- data.frame(cluster = ifelse(test = F.col  == 1, yes = "cluster 1", no = "cluster 2"))
F.col$Primer_name <- rownames(F.col)

F.col <- merge(F.col, PrimTax, by="Primer_name", sort= F)

col_groups <- F.col %>%
  select("Primer_name", "Gen", "Region") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Primer_name

col_groups$Primer_name<- NULL

colour_groups <- list( Gen= c("18S"= "#440154FF", 
                              "28S"= "#21908CFF", 
                              "COI"= "#FDE725FF", 
                              "16S"= "pink",
                              "12S"= "#E3DAC9",
                              "ITS"= "#C46210",
                              "rbcL" = "#D0FF14"))

fungheatmap <- pheatmap(F.intersection, 
                        color = plasma(100),
                        border_color = NA,
                        annotation_col = col_groups, 
                        #annotation_row = col_groups,
                        annotation_colors = colour_groups,
                        #cutree_rows = 2,
                        #cutree_cols = 2,
                        show_rownames = F,
                        show_colnames = F,
                        main= "Redundant genus of fungi amplified")

ggsave("~/AA_Primer_evaluation/Figures/Fungi_redundancy_allprimers.pdf", plot = fungheatmap, dpi = 450, width = 14, height = 10)

####Cummulative curve for fungi (All primers)
num.taxa.fungi <-sapply(c( "species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PM.fungi.comb, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads.fungi <- unlist(lapply(PM.fungi.comb, function (x) 
  sum(otu_table(x))))

PrimTax.fungi <-as.data.frame(cbind(num.taxa.fungi,as.data.frame(num.reads.fungi)))
PrimTax.fungi$Primer_name<- rownames(PrimTax.fungi)
markers <- PrimTax[,c(8,16)]
PrimTax.fungi <- join(PrimTax.fungi, markers)
PrimTax.fungi <- PrimTax.fungi[order(-PrimTax.fungi$num.reads.fungi),] 

ps.numord.l.fungi <- PM.fungi.comb[c(PrimTax.fungi$Primer_name)][order(as.vector(unlist(PrimTax.fungi$num.reads.fungi)), decreasing=TRUE)]

cum.tax.fungi <- sapply(c("species","genus", "family","order", "phylum"), function (rank){
  unique.taxa.cumsum(ps.numord.l.fungi, rank)
})

colnames(cum.tax.fungi) <- paste0("cum.", colnames(cum.tax.fungi))

cum.tax.fungi <- as.data.frame(cbind(1:nrow(cum.tax.fungi), cum.tax.fungi))
cum.tax.fungi[2:28,7] <- PrimTax.fungi$Gen ## add marker amplified 
colnames(cum.tax.fungi)[1] <- "Primer_name"
colnames(cum.tax.fungi)[7] <- "Gen"
cum.tax.fungi<- add_row(cum.tax.fungi, Primer_name= 29:40)
cum.tax.fungi[29:40,2:6]<-cum.tax.fungi[28,2:6] 
cum.tax.fungi[29:40,7]<-col_groups[28:39,1] 

cum.plot.fungi <- ggplot(cum.tax.fungi) +
  geom_step(aes(Primer_name, cum.species), color="red") +
  geom_text(aes(26, 2000, label="Species"), color="red") +
  geom_point(aes(Primer_name, cum.species, colour = Gen, size= 2))+
  scale_color_manual(values = c("#E3DAC9","pink","#440154FF","#21908CFF","#FDE725FF","#C46210","#D0FF14"))+
  geom_step(aes(Primer_name, cum.genus), color="#0D0887FF") +
  geom_text(aes(26, 1000, label="Genera"), color="#0D0887FF") +
  geom_step(aes(Primer_name, cum.family), color="#7E03A8FF") +
  geom_text(aes(26, 440, label="Families"), color="#7E03A8FF") +
  geom_step(aes(Primer_name, cum.order), color="#CC4678FF") +
  geom_text(aes(26, 160, label="Orders"), color="#CC4678FF") +
  #geom_step(aes(Primer_name, cum.class), color="#F0F921FF") +
  #geom_text(aes(20, 50, label="Classes"), color="") +
  geom_step(aes(Primer_name, cum.phylum), color="#F89441FF") +
  geom_text(aes(26, 8, label="Phyla"), color="#F89441FF") +
  scale_y_log10("Cummulative count of taxa (Fungi)") +
  #scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = 'none', axis.title.x = element_blank())+
  labs(tag = "C)")#+

CummParaFung <- grid.arrange(cum.plot.all, cum.plot.para, cum.plot.fungi, nrow= 1, ncol= 3,  
                             bottom= textGrob("Number of primers considerd (starting with the one with highest read count)"))
###Put parasites and fungi together 
ggsave("~/AA_Primer_evaluation/Figures/ParaFung_cummulative_allprimers.pdf", plot = CummParaFung, dpi = 450, width = 14, height = 10)


#lapply(AbPhy, unZeroOne)
#AbPhy <- AbPhy[-c(373),]
#AbPhy$numberCycles <- rep(35,nrow(AbPhy))