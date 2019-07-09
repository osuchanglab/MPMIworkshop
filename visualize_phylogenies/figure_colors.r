if (!require(RColorBrewer)) install.packages('RColorBrewer')
library(RColorBrewer)

#a function to generate color palettes in the ggplot2 default colors
#example, gg_color_hue(30) returns a 30 color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#Species colors
specieslist<-c("G1","G4","G7","G7_Other","G8","G9","larrymoorei","rhizogenes","rubi","skierniewicense","vitis")
speciescolors<-gg_color_hue(11)
names(speciescolors)<-specieslist

#Biovar colors
biovarlist<-c("BiovarI","BiovarII","Brucella","BiovarIII","Rhizobium","Rhizobium_Other","Neorhizobium","Sinorhizobium-Ensifer","undicola","Mesorhizobium", "arsenijevicii", "G2", "blue", "black", "Ochrobactrum", "Shinella", "Aureimonas", "Martelella","rhizogenes-like")
biovarcolors<-c("#E41A1C","#4DAF4A","darkgreen","#377EB8","#04048C","cadetblue3","darkcyan","#984EA3","orchid","gray", "salmon", "purple", "blue", "black", "darkolivegreen1", "orange", "lightgoldenrod1", "pink", "dodgerblue3")
names(biovarcolors)<-biovarlist

#Combined biovar/species colors
taxonomycolors<-c(speciescolors,biovarcolors)
