#load libraries
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(DESeq2); packageVersion("DESeq2")

#####################################
#PCA plot
#####################################
#V3V4
#variance stabilization
V3V4_data.BAC.prev.dds <- phyloseq_to_deseq2(V3V4_data.BAC.prev, ~1)
varianceStabilizingTransformation(V3V4_data.BAC.prev.dds, blind = TRUE, fitType = "parametric")
V3V4_data.BAC.prev.dds <- estimateSizeFactors(V3V4_data.BAC.prev.dds)
V3V4_data.BAC.prev.dds <- estimateDispersions(V3V4_data.BAC.prev.dds)
otu.vst <- getVarianceStabilizedData(V3V4_data.BAC.prev.dds)

V3V4_data.BAC.prev.vst<-V3V4_data.BAC.prev
otu_table(V3V4_data.BAC.prev.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

#ordinaion
V3V4.ord <- ordinate(V3V4_data.BAC.prev.vst, method = "RDA", distance = "euclidean")
V3V4.ord.df <- plot_ordination(V3V4_data.BAC.prev.vst, V3V4.ord, axes = c(1,2,3),justDF = TRUE)

##extract explained variance
V3V4.ord.evals <- 100 * (V3V4.ord$CA$eig/ sum(V3V4.ord$CA$eig))
V3V4.ord.df$ID <- rownames(V3V4.ord.df)

V3V4.ord.p <- ggplot(data = V3V4.ord.df, aes(x =PC1, y=PC2, shape = merging))+
  geom_point(colour="black",size = 4)+
  geom_point(size = 3)+
  #geom_polygon(data=PS107.ord.df,aes(x=NMDS1,y=NMDS2,fill=Type,group=Type),alpha=0.30) +
  #geom_text(aes(label = StationName), colour = "black", nudge_y= -0.5,  size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(V3V4.ord.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(V3V4.ord.evals[2], 2)), shape = "Depth", color = "Community")+
  #annotate(geom="text", size = 5, x=-1.2, y=0.7, label= paste("Stress =", round(PS107.ord$stress, 3), sep = " "),color="black")+
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")

#V4V5
#variance stabilization
V4V5_data.BAC.prev.dds <- phyloseq_to_deseq2(V4V5_data.BAC.prev, ~1)
varianceStabilizingTransformation(V4V5_data.BAC.prev.dds, blind = TRUE, fitType = "parametric")
V4V5_data.BAC.prev.dds <- estimateSizeFactors(V4V5_data.BAC.prev.dds)
V4V5_data.BAC.prev.dds <- estimateDispersions(V4V5_data.BAC.prev.dds)
otu.vst <- getVarianceStabilizedData(V4V5_data.BAC.prev.dds)

V4V5_data.BAC.prev.vst<-V4V5_data.BAC.prev
otu_table(V4V5_data.BAC.prev.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)

#ordinaion
V4V5.ord <- ordinate(V4V5_data.BAC.prev.vst, method = "RDA", distance = "euclidean")
V4V5.ord.df <- plot_ordination(V4V5_data.BAC.prev.vst, V4V5.ord, axes = c(1,2,3),justDF = TRUE)

##extract explained variance
V4V5.ord.evals <- 100 * (V4V5.ord$CA$eig/ sum(V4V5.ord$CA$eig))
V4V5.ord.df$ID <- rownames(V4V5.ord.df)

V4V5.ord.p <- ggplot(data = V4V5.ord.df, aes(x =PC1, y=PC2, shape = merging))+
  geom_point(colour="black",size = 4)+
  geom_point(size = 3)+
  labs(x = sprintf("PC1 [%s%%]", round(V4V5.ord.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(V4V5.ord.evals[2], 2)), shape = "Depth", color = "Community")+
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")

plot_grid(V3V4.ord.p,V4V5.ord.p, ncol = 2)

ggsave("./figures/PCA_prev.pdf", 
       last_plot(), 
       dpi = 300, 
       #width = 11.4, height = 23, 
       units = "cm")


#check whether the environments definition are significant
V3V4_metadata <- as(sample_data(V3V4_data.BAC.prev.vst), "data.frame")
V3V4.dist <- phyloseq::distance(V3V4_data.BAC.prev.vst, "euclidean")
V3V4_comm.adonis <- adonis(V3V4.dist ~  merging, V3V4_metadata)
V3V4_comm.adonis

V4V5_metadata <- as(sample_data(V4V5_data.BAC.prev.vst), "data.frame")
V4V5.dist <- phyloseq::distance(V4V5_data.BAC.prev.vst, "euclidean")
V4V5_comm.adonis <- adonis(V4V5.dist ~  merging, V4V5_metadata)
V4V5_comm.adonis


##############################
# Mantel test between the datasets
##############################
#create sample names to merge between the datasets
sample_names(V3V4_data.BAC.prev.vst) <- paste(sample_data(V3V4_data.BAC.prev.vst)$Expedition,sample_data(V3V4_data.BAC.prev.vst)$Station,sample_data(V3V4_data.BAC.prev.vst)$Depth,sample_data(V3V4_data.BAC.prev.vst)$Type, sep="+")
sample_names(V4V5_data.BAC.prev.vst) <- paste(sample_data(V4V5_data.BAC.prev.vst)$Expedition,sample_data(V4V5_data.BAC.prev.vst)$Station,sample_data(V4V5_data.BAC.prev.vst)$Depth,sample_data(V4V5_data.BAC.prev.vst)$Type, sep="+")

#export otu tables
V3V4_otu.vst <- t(otu_table(V3V4_data.BAC.prev.vst))
V3V4_otu.vst <- V3V4_otu.vst[sort(rownames(V3V4_otu.vst)),]

V4V5_otu.vst <- t(otu_table(V4V5_data.BAC.prev.vst))
V4V5_otu.vst <- V4V5_otu.vst[sort(rownames(V4V5_otu.vst)),]

#calculate distances
V3V4_data.BAC.prev.vst.dist <- vegdist(V3V4_otu.vst, 
                                       method="euclidean")

V4V5_data.BAC.prev.vst.dist <- vegdist(V4V5_otu.vst, 
                                       method="euclidean")


#run mantel test
mantel(V3V4_data.BAC.prev.vst.dist, V4V5_data.BAC.prev.vst.dist, method="pearson", permutations=999, strata = NULL,
       na.rm = FALSE, parallel = getOption("mc.cores"))


