#library(BioPhysConnectoR) It doesn't exist anymore :S 
library(ggplot2)
library(viridis)
library(parallel)
library(DECIPHER)
library(seqinr)
library(Bios2cor)


## Why is this needed here for the primers agains Silva?
## if(!exists("primerF")){
##   source("R/1_generalAA.R")
## }

## had to change this path
aln <- read.fasta("/SAN/db/RDP_Silva/Silva_123/silva.nr_v123_EUK.align")

test <- aln[1:10]
e <- entropy(test)
o <- omes(test)
o <- o$normalized

entropy_graph(entropy = e, o, filter = NULL, elite = 25, high = 275, csv = TRUE, 
              name = "Figures/test_entropy_graph.pdf")

keep <- !apply(aln, 2, function (x) all(x %in% c("-", ".")) )

aln <- aln$ali[,keep]
ent <- get.entropy(aln, gapchar = "-")
ent <- data.frame(entropy=ent, plain.length=1:length(ent))

## get the length of non gap characters for each sequence in the
## alingment
non.gap.lenght <- t(apply(aln, 1, function (x) {
  cumsum(x != "-")
}))
mean.no.gap.length <- colMeans(non.gap.lenght)

ent$no.gap.l <- mean.no.gap.length


ent$perc.gap <- apply(aln, 2,
                      function(x) sum(x=="-")/length(x) *100)


ent$trimmed.ent <- ifelse(ent$perc.gap < 50, ent$entropy, NA)

ent$smooth.ent <- ent$trimmed.ent

ent$smooth.ent[!is.na(ent$trimmed.ent)] <-
  runmed(ent$trimmed.ent[!is.na(ent$trimmed.ent)], 21)


## mapping the primers against sequence

matches <- mclapply(seq_along(primerF), function (i){
  pmatches <- apply(aln, 1, function (x) {
    nogap <- x[x != "-"]
    seq <- DNAString(paste(nogap, collapse=""))
    hitF <- matchPattern(primerF[[i]],
                         seq, fixed=FALSE)
    hitR <- matchPattern(reverseComplement(DNAString(primerR[[i]])),
                         seq, fixed=FALSE)
    list(start(hitF), end(hitR))
  })
  return(pmatches)
}, mc.cores = 20)

names(matches) <- paste(names(primerF), names(primerR), sep=".")

mat.list <- lapply(matches, function (x) do.call(rbind, x))

## something like this for getting the tax scope of forward and
## reverse primer binding

## foo <- lapply(mat.list, function(x) !isEmpty(x[, 1]) & !isEmpty(x[,
## 2]))


## for now only hte start and end points
meanStart <- unlist(lapply(mat.list,
                           function(x) mean(unlist(x[, 1]), na.rm=TRUE)))

meanEnd <- unlist(lapply(mat.list,
                         function(x) mean(unlist(x[, 2]), na.rm=TRUE)))

Pranges <- as.data.frame(cbind(meanStart, meanEnd))

## DODGY FIX for EukB... watch out for this primer!
Pranges[grepl(".EukB", rownames(Pranges)), ]$meanEnd <-
  max(ent$no.gap.l)
## DODGY FIX for Medin... watch out for this primer!
Pranges[grepl("Medin.", rownames(Pranges)), ]$meanStart <-
  min(ent$no.gap.l)

Pranges <- Pranges[order(Pranges$meanStart, Pranges$meanEnd), ]

Pranges$y.pos <- seq(2.5, 100, by = 0.1)[1:nrow(Pranges)]

Pranges <- merge(Pranges, PrimTax, by=0)

pdf("figures/entropy_primers_norm.pdf", width=10, height=6)
ggplot(ent, aes(no.gap.l, trimmed.ent)) +
  geom_hex(binwidth=c(31, 0.1)) +
  scale_fill_viridis(option = "viridis") +
  geom_segment(mapping=aes(x = meanStart, y = y.pos, xend = meanEnd,
                           yend = y.pos, color=log10(num.reads)),
               size=2,
               data=Pranges)+
  scale_color_viridis(option = "plasma")+
  geom_text(aes(x = meanStart, y = y.pos+0.04, label = Row.names), Pranges, size=2) +
  scale_x_continuous("mean non-gapped alignement length")+
  theme_bw()
dev.off()
## EukB Primers are at wrong location!


devtools::source_gist("524eade46135f6348140",
                      filename = "ggplot_smooth_func.R")

pdf("figures/size_vs_num.pdf", width=8, height=6)
ggplot(Pranges, aes(meanEnd-meanStart, num.reads)) +
  geom_point(aes(size=Genus, color=Phylum))+
  scale_color_viridis(option = "plasma")+
  scale_x_continuous("lenght of amplicon") +
  scale_y_log10("number of sequencing reads")+
  stat_smooth_func(geom="text", method="lm", hjust = -1.5,  parse=TRUE) +
  stat_smooth(method="lm", se=FALSE) +
  annotation_logticks(sides="l") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
dev.off()
