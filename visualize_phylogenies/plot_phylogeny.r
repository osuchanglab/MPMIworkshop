#If the required R packages are not found, install them
#If they are installed you can delete these three lines
if (!require(phangorn)) install.packages('phangorn')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(ggtree)) install.packages('ggtree')

#load the required packages
library(phangorn)
library(ggtree)
library(ggplot2)

#load color definitions from a separate file so you only have to change them in one place
source("./figure_colors.r")

#load a table of strain information/phenotypes to plot on the tree
datasettable <- read.table("./strainphenotypes.txt", sep = "\t", header=T, stringsAsFactors = F, quote = "")


#load a tree file in newick or nexus tree format
mytree <- read.tree("./myphylogeny.tre")

#midpoint root the tree
mytreemid<-phangorn::midpoint(mytree,node.labels='label')


#######################################################################
#
#  Plot a phylogeny with nothing on it
#
#######################################################################

#save visual plot output as a pdf (png also works)
pdf("Figure_mlsa_blank.pdf",width=12,height=12)
#plot the phylogeny using ggtree + ggplot2
ggtree(mytreemid, ladderize = FALSE, layout="circular")  
#finish the pdf figure
dev.off()


#######################################################################
#
#   Plot a phylogeny with tip labels
#
#######################################################################

pdf("Figure_mlsa_tiplabels.pdf",width=12,height=12)

ggtree(mytreemid, ladderize = FALSE, layout="circular") %<+% datasettable +
    geom_text2(aes(subset=!isTip, label=label),size=2) +
    geom_tiplab2(aes(label=label,angle=angle),offset=0.01,size=4) 

dev.off()


#######################################################################
#
#  Plot a phylogeny with colored points at the tips
#
#######################################################################

pdf("Figure_mlsa_points.pdf",width=12,height=12)

ggtree(mytreemid, ladderize = FALSE, layout="circular") %<+% datasettable +
    geom_text2(aes(subset=!isTip, label=label),size=2) +
    geom_tippoint(aes(color=Genomospecies), size=3) +
    scale_fill_manual("Taxonomy", values=c("Unknown"="gray", "Unknown"="gray")) +
    scale_color_manual("Taxonomy",values = taxonomycolors) +
    theme(legend.position=c(0.9, 0.2), plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

dev.off()


#######################################################################
#
#  Plot a phylogeny with a heatmap or gene presence/absence information
#
#######################################################################

#load a gene presence/absence table
genetable <- read.table("./genetable.txt", sep = "\t", header=T, stringsAsFactors = F, quote = "")
row.names(genetable) <- genetable$strain
genetable$strain <- NULL
genetable[]<-lapply(genetable,function(x) {factor(x,levels=c(1,0),labels=c("Present","Absent"))})

genecolors<-c("dodgerblue3", "gray90")
names(genecolors)<-c("Present","Absent")

allcolors<-c(genecolors,taxonomycolors)

pdf("Figure_mlsa_heatmap.pdf",width=12,height=12)

p<-ggtree(mytreemid, ladderize = FALSE, layout="circular") %<+% datasettable +
    geom_text2(aes(subset=!isTip, label=label),size=2) +
    geom_tippoint(aes(color=Genomospecies), size=3) +
    scale_color_manual("Taxonomy",values = allcolors) +
    theme(legend.position=c(0.9, 0.2), plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))

gheatmap(p, genetable, colnames=TRUE, colnames_offset_y=0, colnames_position="top", colnames_angle=90, font.size=2,width=0.3,color="black") + scale_fill_manual("Taxonomy",values = allcolors) 

dev.off()
