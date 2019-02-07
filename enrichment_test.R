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
library(DESeq2); packageVersion("DESeq2")

#load scripts
source("./scripts/color_palettes.R")


#load datasets 
V3V4_data.BAC <- readRDS("./Data/V3V4_data_BAC.rds")
V4V5_data.BAC <- readRDS("./Data/V4V5_data_BAC.rds")

#calculate relative abundance and agglomerate to family level 
#V3V4
V3V4_data.BAC.ra.glom <- tax_glom(V3V4_data.BAC, taxrank = "family")
V3V4_data.BAC.ra.glom <-prune_taxa(taxa_sums(V3V4_data.BAC.ra.glom)>100,V3V4_data.BAC.ra.glom)
V3V4_data.BAC.melt <- psmelt(V3V4_data.BAC.ra.glom)

#V4V5
V4V5_data.BAC.ra.glom <- tax_glom(V4V5_data.BAC, taxrank = "family")
V4V5_data.BAC.ra.glom <-prune_taxa(taxa_sums(V4V5_data.BAC.ra.glom)>100,V4V5_data.BAC.ra.glom)
V4V5_data.BAC.melt <- psmelt(V4V5_data.BAC.ra.glom)

#conduct enrichment test for each environment separately
environment <- c("Sea-ice","Surface-water","Deep-water","Sediment")
merged_enrichment_tests <- data.frame()
#run DEseq2
for (i in 1:4){
#subset the selected environment
V3V4_data.BAC.melt.sub <- V3V4_data.BAC.melt[V3V4_data.BAC.melt$merging == environment[i],]
V4V5_data.BAC.melt.sub <- V4V5_data.BAC.melt[V4V5_data.BAC.melt$merging == environment[i],]

#generate abundance tables
#V3V4
V3V4_data.BAC.melt.sub$Sample <- paste("V3V4",V3V4_data.BAC.melt.sub$Station,
                                       V3V4_data.BAC.melt.sub$Type,
                                       V3V4_data.BAC.melt.sub$Depth, sep = ".")
#generate abundance table for Family*Sample
V3V4_data.BAC.agg <- aggregate(Abundance~Sample+family, V3V4_data.BAC.melt.sub,sum)
V3V4_data.BAC.agg[,c("Sample","family","Abundance")] %>% spread(Sample, Abundance)-> V3V4_data.BAC.wide

#V4V5
V4V5_data.BAC.melt.sub$Sample <- paste("V4V5",V4V5_data.BAC.melt.sub$Station,
                                       V4V5_data.BAC.melt.sub$Type,
                                       V4V5_data.BAC.melt.sub$Depth, sep = ".")

#generate abundance table for Family*Sample
V4V5_data.BAC.agg <- aggregate(Abundance~Sample+family, V4V5_data.BAC.melt.sub,sum)
V4V5_data.BAC.agg[,c("Sample","family","Abundance")] %>% spread(Sample, Abundance)-> V4V5_data.BAC.wide

#merge abundance tables
deseq.test <- merge(V3V4_data.BAC.wide,V4V5_data.BAC.wide, by= "family")
rownames(deseq.test) <- deseq.test$family
deseq.test <- subset(deseq.test, select= -c(family))

#run DEseq on sed fraction
primer <- factor(c(rep("V3V4",dim(deseq.test)[2]/2),rep("V4V5",dim(deseq.test)[2]/2)))
dds <- DESeqDataSetFromMatrix(deseq.test, DataFrame(primer),~primer)
varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
BAC_sed.ddsMat <- estimateSizeFactors(dds)
BAC_sed.ddsMat <- estimateDispersions(BAC_sed.ddsMat)
BAC_sed.DEseq <- DESeq(BAC_sed.ddsMat, fitType="parametric")
BAC_sed.DEseq.res <- results(BAC_sed.DEseq)

#extract only significant OTU
BAC_sed.DEseq.res.sig <- BAC_sed.DEseq.res[which(BAC_sed.DEseq.res$padj < 0.1), ]
BAC_sed.DEseq.res.sig.sub <- BAC_sed.DEseq.res.sig[abs(BAC_sed.DEseq.res.sig$log2FoldChange)>1,]
BAC_sed.DEseq.res.sig.sub$Env <- environment[i]

#add taxonomy
BAC_sed.DEseq.res.sig.sub$family <-rownames(BAC_sed.DEseq.res.sig.sub) 
BAC_sed.DEseq.res.sig.sub <- as(BAC_sed.DEseq.res.sig.sub,"data.frame")
BAC_sed.DEseq.res.sig.sub <- merge(BAC_sed.DEseq.res.sig.sub, unique(V3V4_data.BAC.melt[,c("family","class")]), by = "family")
rownames(BAC_sed.DEseq.res.sig.sub) <- NULL


merged_enrichment_tests <- rbind(merged_enrichment_tests,BAC_sed.DEseq.res.sig.sub)
}

merged_enrichment_tests$class <- as.character(merged_enrichment_tests$class)
merged_enrichment_tests$class[!merged_enrichment_tests$class %in% selected_taxa] <- "Other"

merged_enrichment_tests$class <- factor(merged_enrichment_tests$class, ordered = TRUE,
                                          levels = unique(merged_enrichment_tests$class[order(merged_enrichment_tests$class)]))

merged_enrichment_tests$family <- factor(merged_enrichment_tests$family, ordered = TRUE,
                                           levels = unique(merged_enrichment_tests$family[order(merged_enrichment_tests$class)]))

#plot
BAC_daOTU.p <- ggplot(data= merged_enrichment_tests[merged_enrichment_tests$class=="Other",], aes(x = log2FoldChange, y = family))+ 
  #geom_point(data= BAC_sed.DEseq.res.sig.sub[BAC_sed.DEseq.res.sig.sub$class %in% selected_taxa,], aes(x = family, y = log2FoldChange), colour = "black",size = 6)+
  #geom_point(data= BAC_sed.DEseq.res.sig.sub[BAC_sed.DEseq.res.sig.sub$class %in% selected_taxa,], aes(x = family, y = log2FoldChange, colour = class), size = 4)+
  #geom_point(data= BAC_sed.DEseq.res.sig.sub[!BAC_sed.DEseq.res.sig.sub$class %in% selected_taxa,], aes(x = family, y = log2FoldChange), colour = "black", alpha = 0.5, size = 3)+
  geom_point(data= merged_enrichment_tests[merged_enrichment_tests$class=="Other",], aes(x = log2FoldChange, y = family), colour = "black",size = 5)+
  geom_point(data= merged_enrichment_tests[merged_enrichment_tests$class=="Other",], aes(x = log2FoldChange, y = family, colour = class), size = 4)+
  ylab("log2foldchange")+
  xlim(-10,10)+
  geom_vline(aes(xintercept=0), linetype="dashed")+
  #geom_hline(aes(yintercept=11), linetype="solid")+
  #geom_hline(aes(yintercept=-11), linetype="solid")+
  scale_colour_manual(values = taxa_col)+
  #facet_grid(Env~.)+
  theme_classic(base_size = 10)+
  #scale_x_discrete(labels=unique(deseq_res_all[order(deseq_res_all$class),"family"]))+
  theme(axis.text.x =element_text(angle=90) , legend.position = "bottom")

ggsave("./figures/enrichments_SI.pdf",
       plot = BAC_daOTU.p,
       scale = 1,
       units = "cm",
       width = 30,
       height = 30,
       dpi = 300)

as.data.frame(as.list(aggregate(pvalue~class,merged_enrichment_tests,FUN = function(x) count=length(x))))
