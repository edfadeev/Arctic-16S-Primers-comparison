#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries 
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(cowplot); packageVersion("cowplot")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(gridExtra); packageVersion("gridExtra")

#load scripts
source("./scripts/color_palettes.R")

#load datasets 
V3V4_data.BAC <- readRDS("./Data/V3V4_data_BAC.rds")
V4V5_data.BAC <- readRDS("./Data/V4V5_data_BAC.rds")

#####################################
#Preprocess bacterial OTU by prevalence of each taxa
#####################################
#V3V4
#in how many samples did each taxa appear at least once
prev0 = apply(X = otu_table(V3V4_data.BAC),
              MARGIN = ifelse(taxa_are_rows(V3V4_data.BAC), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(V3V4_data.BAC)

# Execute prevalence filter, using `prune_taxa()` function
V3V4_data.BAC.pruned <-  prune_taxa((prev0 > prevalenceThreshold), V3V4_data.BAC)

#V4V5
#in how many samples did each taxa appear at least once
prev0 = apply(X = otu_table(V4V5_data.BAC),
              MARGIN = ifelse(taxa_are_rows(V4V5_data.BAC), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(V4V5_data.BAC)

# Execute prevalence filter, using `prune_taxa()` function
V4V5_data.BAC.pruned <-  prune_taxa((prev0 > prevalenceThreshold), V4V5_data.BAC)

#####################################
#Bar plots
#####################################
selected_taxa <- c("Alphaproteobacteria",
                   "Acidimicrobiia",
                   "Bacteroidia",
                   "Deltaproteobacteria",
                   "Dehalococcoidia",
                   "Gammaproteobacteria",
                   "Marinimicrobia (SAR406 clade)_unclassified", 
                   "Nitrospinia", 
                   "Planctomycetacia",
                   "Verrucomicrobiae")

#V3V4
#calculate relative abundance
V3V4_data.BAC.ra <- transform_sample_counts(V3V4_data.BAC.pruned, function(x) x / sum(x) )
#merge the extracted OTU
V3V4_data.BAC.melt <- psmelt(V3V4_data.BAC.ra)
names(V3V4_data.BAC.melt)[names(V3V4_data.BAC.melt) == "class"] <- "Taxonomy"
V3V4_data.BAC.melt$Taxonomy <- as.character(V3V4_data.BAC.melt$Taxonomy)
V3V4_data.BAC.melt$Taxonomy[!V3V4_data.BAC.melt$Taxonomy %in% selected_taxa] <- "Other"
V3V4_data.BAC.melt <- arrange(V3V4_data.BAC.melt, Sample, Taxonomy)
V3V4_data.BAC.melt$SampleID <- paste(V3V4_data.BAC.melt$Expedition,V3V4_data.BAC.melt$Station,V3V4_data.BAC.melt$Depth,V3V4_data.BAC.melt$Type, sep="+")

#aggregate by taxonomy
V3V4_data.BAC.melt.agg <- aggregate(Abundance~SampleID+Taxonomy+merging, V3V4_data.BAC.melt, FUN= sum)
V3V4_data.BAC.melt.agg$Abundance <- V3V4_data.BAC.melt.agg$Abundance*100
V3V4_data.BAC.melt.agg$Taxonomy <- as.factor(V3V4_data.BAC.melt.agg$Taxonomy)
V3V4_data.BAC.melt.agg$Primer <- "V3V4"


#V4V5
#calculate relative abundance
V4V5_data.BAC.ra <- transform_sample_counts(V4V5_data.BAC.pruned, function(x) x / sum(x) )
#merge the extracted OTU
V4V5_data.BAC.melt <- psmelt(V4V5_data.BAC.ra)
names(V4V5_data.BAC.melt)[names(V4V5_data.BAC.melt) == "class"] <- "Taxonomy"
V4V5_data.BAC.melt$Taxonomy <- as.character(V4V5_data.BAC.melt$Taxonomy)
V4V5_data.BAC.melt$Taxonomy[!V4V5_data.BAC.melt$Taxonomy %in% selected_taxa] <- "Other"
V4V5_data.BAC.melt <- arrange(V4V5_data.BAC.melt, Sample, Taxonomy)
V4V5_data.BAC.melt$SampleID <- paste(V4V5_data.BAC.melt$Expedition,V4V5_data.BAC.melt$Station,V4V5_data.BAC.melt$Depth, V4V5_data.BAC.melt$Type, sep="+")

#aggregate by taxonomy
V4V5_data.BAC.melt.agg <- aggregate(Abundance~SampleID+Taxonomy+merging, V4V5_data.BAC.melt, FUN= sum)
V4V5_data.BAC.melt.agg$Abundance <- V4V5_data.BAC.melt.agg$Abundance*100
V4V5_data.BAC.melt.agg$Taxonomy <- as.factor(V4V5_data.BAC.melt.agg$Taxonomy)
V4V5_data.BAC.melt.agg$Primer <- "V4V5"

BAC.melt.agg <- rbind(V3V4_data.BAC.melt.agg,V4V5_data.BAC.melt.agg)
#BAC.melt.agg %>% separate(SampleID, c("Expedition", "Station","Depth","Type"),"\\+") ->BAC.melt.agg

BAC.melt.combined.agg <- aggregate(Abundance~SampleID+merging+Primer+Taxonomy, BAC.melt.agg, mean)
BAC.melt.combined.agg$Taxonomy <- gsub("_unclassified| Incertae Sedis", "",BAC.melt.combined.agg$Taxonomy) 

BAC.melt.combined.agg$Taxonomy <- factor(BAC.melt.combined.agg$Taxonomy, levels= c("Alphaproteobacteria","Acidimicrobiia","Bacteroidia","Deltaproteobacteria","Dehalococcoidia","Gammaproteobacteria","Marinimicrobia (SAR406 clade)", 
                                                                           "Nitrospinia", "Planctomycetacia","Verrucomicrobiae","Other"))
BAC.melt.combined.agg$merging <- factor(BAC.melt.combined.agg$merging, levels= c("Sea-ice", "Surface-water","Deep-water","Sediment"))

p <- list()
for (i in levels(BAC.melt.combined.agg$merging)){
  sub <- BAC.melt.combined.agg[BAC.melt.combined.agg$merging == i,]
  p[[i]] <-ggplot(sub, aes(x = SampleID, y = Abundance, fill = Taxonomy)) + 
    facet_grid(~Primer) + labs(title = i)+
    geom_bar(stat = "identity", position="fill", colour = "black") +
    scale_y_continuous(breaks=NULL)+
    scale_fill_manual(values = taxa_col) +
    #coord_polar("y", start=0)+
    theme_classic(base_size = 10)+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          strip.background = element_blank())
}

plot_grid(p[[1]],p[[2]],p[[3]], p[[4]], ncol = 2, align = "hv")

ggsave("./figures/bar_plots.pdf",
       plot = ggplot2::last_plot(),
       scale = 1,
       units = "cm",
       width = 17.8,
       #height = 17.4,
       dpi = 300)

#####################################
#Statistical tests for representation of the different groups between the primers
#####################################+
library("ggpubr")
#test normal distribution
shapiro.test(BAC.melt.combined.agg[BAC.melt.combined.agg$Primer=="V3V4",]$Abundance) #if significant not normal
#density plot
ggdensity(V3V4_data.BAC.melt.agg$Abundance, 
          main = "",
          xlab = "")

ggdensity(V4V5_data.BAC.melt.agg$Abundance, 
          main = "",
          xlab = "")

#because the abundance distribution is not normal, using wilcox
wilcox.test(V3V4_data.BAC.melt.agg$Abundance, 
            V4V5_data.BAC.melt.agg$Abundance)


#test for each environment
wilcox.test(V3V4_data.BAC.melt.agg[V3V4_data.BAC.melt.agg$merging =="Sea-ice",]$Abundance, 
            V4V5_data.BAC.melt.agg[V4V5_data.BAC.melt.agg$merging =="Sea-ice",]$Abundance)

wilcox.test(V3V4_data.BAC.melt.agg[V3V4_data.BAC.melt.agg$merging =="Surface-water",]$Abundance, 
            V4V5_data.BAC.melt.agg[V4V5_data.BAC.melt.agg$merging =="Surface-water",]$Abundance)

wilcox.test(V3V4_data.BAC.melt.agg[V3V4_data.BAC.melt.agg$merging =="Deep-water",]$Abundance, 
            V4V5_data.BAC.melt.agg[V4V5_data.BAC.melt.agg$merging =="Deep-water",]$Abundance)

wilcox.test(V3V4_data.BAC.melt.agg[V3V4_data.BAC.melt.agg$merging =="Sediment",]$Abundance, 
            V4V5_data.BAC.melt.agg[V4V5_data.BAC.melt.agg$merging =="Sediment",]$Abundance)

#test for each taxonomic group
sapply(selected_taxa, function(i){
  wilcox.test(V3V4_data.BAC.melt.agg[V3V4_data.BAC.melt.agg$Taxonomy == i,]$Abundance, 
              V4V5_data.BAC.melt.agg[V4V5_data.BAC.melt.agg$Taxonomy == i,]$Abundance, exact=FALSE )$p.value
}) 

ggdensity(V3V4_data.BAC.melt.agg[V3V4_data.BAC.melt.agg$Taxonomy == "Deltaproteobacteria",]$Abundance, 
          main = "",
          xlab = "")




