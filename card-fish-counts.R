#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

library("ggplot2")
library("tidyr")
library("reshape2")

#load scripts
source("./scripts/color_palettes.R")

#import raw data
raw.counts <- read.csv("../CARD-FISH/card-fish-counts.csv", dec = ".", stringsAsFactors = FALSE)
raw.counts$conc.FISH <- as.numeric(raw.counts$conc.FISH)


#summarize counts
counts.agg <- as.data.frame(as.list(aggregate(conc.FISH~StationName+Depth+Domain,raw.counts,FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))))

#overview plot
counts.agg.abs <- counts.agg[counts.agg$Domain!="EUB",]

FISH.plot <- ggplot() + 
  geom_bar(aes(y = conc.FISH.mean, x = Domain, fill =Domain ), colour = "black", data = counts.agg.abs, stat="identity")+
  facet_grid(Depth~StationName)+
  #scale_y_log10()+
  scale_fill_manual(values = cbbPalette)+
  theme_classic(base_size = 12)

levels(counts.agg.abs$Depth) <- c("ice","DCM","MESO","SED")
counts.agg.abs$Domain <-  factor(counts.agg.abs$Domain,
                                 levels= c("SAR11","BACT","POL","GAM","ALT"))


p <- list()
for (i in levels(counts.agg.abs$Depth)){
  p[[i]] <- ggplot()+
    geom_bar(aes(y = conc.FISH.mean, x = Domain, fill =Domain), colour = "black", data = counts.agg.abs[counts.agg.abs$Depth==i,], stat="identity")+
    facet_grid(.~StationName)+
    #scale_y_log10()+
    ggtitle(label = i)+
    geom_errorbar(aes(x = Domain, ymin = conc.FISH.mean-conc.FISH.se, ymax = conc.FISH.mean+conc.FISH.se), data = counts.agg.abs[counts.agg.abs$Depth==i,], width=.2,position=position_dodge(.9)) +
    scale_fill_manual(values = cbbPalette)+
    theme_classic(base_size = 12)+
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(),
          legend.position="none",
          panel.grid.major=element_blank(),
          strip.background = element_blank())+
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    
}
ml <- marrangeGrob(p, nrow=4, ncol=1)
ml
ggsave("./figures/comparison_FISH.pdf",
       plot = ml,
       scale = 1,
       units = "cm",
       #width = 17.8,
       #height = 17.4,
       dpi = 300)



#calculate count proportions 
subset(counts.agg, select = -conc.FISH.se)%>% spread(Domain, conc.FISH.mean)-> counts.agg.wide

counts.agg.wide.ra <- counts.agg.wide
counts.agg.wide.ra$SAR11 <- counts.agg.wide.ra$SAR11/counts.agg.wide.ra$EUB
counts.agg.wide.ra$GAM <- (counts.agg.wide.ra$GAM-counts.agg.wide.ra$ALT)/counts.agg.wide.ra$EUB
counts.agg.wide.ra$ALT <- counts.agg.wide.ra$ALT/counts.agg.wide.ra$EUB
counts.agg.wide.ra$BACT <- (counts.agg.wide.ra$BACT-counts.agg.wide.ra$POL)/counts.agg.wide.ra$EUB
counts.agg.wide.ra$POL <- counts.agg.wide.ra$POL/counts.agg.wide.ra$EUB

counts.agg.ra <- melt(subset(counts.agg.wide.ra, select = -EUB))


FISH.ra.plot <- ggplot() + 
  geom_bar(aes(y = value, x = variable, fill =variable ), colour = "black", data = counts.agg.ra, stat="identity")+
  facet_grid(Depth~StationName)+
  scale_fill_manual(values = cbbPalette)+
  theme_classic(base_size = 12)


#get the relevant 16S samples
#V3V4
V3V4_data.order.ra <- transform_sample_counts(V3V4_data.BAC.pruned, function(x) x / sum(x) )
V3V4_data.order.melt <- psmelt(V3V4_data.order.ra)

V3V4_data.melt.fish <- na.omit(V3V4_data.order.melt[V3V4_data.order.melt$FISH==1,])

V3V4_class.agg <- aggregate(Abundance~Station+merging+class,V3V4_data.melt.fish, sum)
names(V3V4_class.agg) <- c("Station","Env","Taxa","Abundance")
V3V4_order.agg <- aggregate(Abundance~Station+merging+order,V3V4_data.melt.fish, sum)
names(V3V4_order.agg) <- c("Station","Env","Taxa","Abundance")
V3V4_family.agg <- aggregate(Abundance~Station+merging+family,V3V4_data.melt.fish, sum)
names(V3V4_family.agg) <- c("Station","Env","Taxa","Abundance")
V3V4_genus.agg <- aggregate(Abundance~Station+merging+genus,V3V4_data.melt.fish, sum)
names(V3V4_genus.agg) <- c("Station","Env","Taxa","Abundance")
V3V4_genus.agg$Abundance <- gsub("Polaribacter.*|\\[Polaribacter.*","Polaribacter",V3V4_genus.agg$Abundance)
                                 
V3V4_summary.agg <-rbind(V3V4_class.agg,V3V4_order.agg,V3V4_family.agg,V3V4_genus.agg)
V3V4_summary.agg$Method <- "V3V4"


#V4V5
V4V5_data.order.ra <- transform_sample_counts(V4V5_data.BAC.pruned, function(x) x / sum(x) )
V4V5_data.order.melt <- psmelt(V4V5_data.order.ra)

V4V5_data.melt.fish <- na.omit(V4V5_data.order.melt[V4V5_data.order.melt$FISH==1,])

V4V5_class.agg <- aggregate(Abundance~Station+merging+class,V4V5_data.melt.fish, sum)
names(V4V5_class.agg) <- c("Station","Env","Taxa","Abundance")
V4V5_order.agg <- aggregate(Abundance~Station+merging+order,V4V5_data.melt.fish, sum)
names(V4V5_order.agg) <- c("Station","Env","Taxa","Abundance")
V4V5_family.agg <- aggregate(Abundance~Station+merging+family,V4V5_data.melt.fish, sum)
names(V4V5_family.agg) <- c("Station","Env","Taxa","Abundance")

V4V5_data.melt.fish$genus <- gsub("Polaribacter.*|\\[Polaribacter.*","Polaribacter",V4V5_data.melt.fish$genus)
V4V5_genus.agg <- aggregate(Abundance~Station+merging+genus,V4V5_data.melt.fish, sum)
names(V4V5_genus.agg) <- c("Station","Env","Taxa","Abundance")

V4V5_summary.agg <-rbind(V4V5_class.agg,V4V5_order.agg,V4V5_family.agg,V4V5_genus.agg)
V4V5_summary.agg$Method <- "V4V5"

#extract the taxonomic groups with counts
FISH_taxa <- c("SAR11 clade", "Alteromonadales",
               "Bacteroidia","Flavobacteriaceae","Gammaproteobacteria")

#FISH_taxa <- c("SAR11 clade", "Alteromonadales","Polaribacter")
V3V4_summary.agg.sub <- V3V4_summary.agg[V3V4_summary.agg$Taxa %in% FISH_taxa,]
V4V5_summary.agg.sub <- V4V5_summary.agg[V4V5_summary.agg$Taxa %in% FISH_taxa,]


summary.agg.sub <- rbind(V3V4_summary.agg.sub,V4V5_summary.agg.sub)
summary.agg.sub$Abundance <- as.numeric(summary.agg.sub$Abundance)
summary.agg.sub$Taxa <- factor(summary.agg.sub$Taxa,
                               levels = c("SAR11 clade", "Bacteroidia","Flavobacteriaceae",
                                          "Gammaproteobacteria", "Alteromonadales"))
q <- list()
for (i in levels(summary.agg.sub$Env)){
  q[[i]] <- ggplot() + 
  geom_bar(aes(y = Abundance, x = Taxa, fill =Taxa ), colour = "black", 
           data = summary.agg.sub[summary.agg.sub$Env==i,], stat="identity")+
  facet_grid(~Method+Station)+
  #ylim(0,0.8)+
  scale_fill_manual(values = cbbPalette)+
    ggtitle(label = i)+
    theme_classic(base_size = 12)+
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(),
          legend.position="bottom",
          panel.grid.major=element_blank(),
          strip.background = element_blank())
}

ml <- marrangeGrob(q, nrow=4, ncol=1)
ml

ggsave("./figures/comparison_V3V4.pdf",
       plot = ml,
       scale = 1,
       units = "cm",
       #width = 17.8,
       #height = 17.4,
       dpi = 300)
