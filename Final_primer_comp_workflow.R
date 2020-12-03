#load libraries
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(iNEXT); packageVersion("iNEXT")
library(rstatix); packageVersion("rstatix")
library(ggpubr); packageVersion("ggpubr")
library(vegan); packageVersion("vegan")
library(VennDiagram); packageVersion("VennDiagram")

set.seed(1586)
theme_set(theme_bw())

#function to add total number of observations to a ggboxplot
give.n <- function(x){
  return(c(y = median(x)*1.15, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

#load colour pallette
source("./scripts/colours.R")

#####################################
#Import dada2 output into phyloseq
#####################################
#load V3V4 dataset
OTU<- read.csv("./data/V3V4/dada2_seqtab_nochim2.txt", h=T, sep="\t")
#correct sample names in OTU table
colnames(OTU)<- gsub("_F_filt.fastq.gz","",colnames(OTU))

TAX<- as.matrix(read.csv("./data/V3V4/dada2_taxonomy_table.txt", h=T,sep = "\t"))

#metadata
ENV <- read.csv("./data/V3V4/V3V4_samples_list.csv", header = TRUE)
rownames(ENV)<- paste("X",ENV$SampleID, sep ="")

# Check order of samples
all.equal(rownames(OTU), rownames(TAX))

#creating Phyloseq dataset
OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
meta <- sample_data(ENV)
V3V4_ps <- phyloseq(OTU, TAX, meta)

#add reference sequence and replace variants with ASVs
dna_V3V4 <- Biostrings::DNAStringSet(taxa_names(V3V4_ps))
names(dna_V3V4) <- taxa_names(V3V4_ps)
V3V4_ps <- merge_phyloseq(V3V4_ps, dna_V3V4)
taxa_names(V3V4_ps) <- paste0("V3V4_ASV", seq(ntaxa(V3V4_ps)))

#load V4V5 dataset
OTU<- read.csv("./data/V4V5/dada2_seqtab_nochim2.txt", h=T, sep="\t")
#correct sample names in OTU table
colnames(OTU)<- gsub("_F_filt.fastq.gz","",colnames(OTU))

TAX<- as.matrix(read.csv("./data/V4V5/dada2_taxonomy_table.txt", h=T,sep = "\t"))

#metadata
ENV <- read.csv("./data/V4V5/V4V5_samples_list.csv", header = TRUE)
rownames(ENV)<- paste("X",ENV$SampleID, sep ="")

# Check order of samples
all.equal(rownames(OTU), rownames(TAX))

#creating Phyloseq dataset
OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
meta <- sample_data(ENV)
V4V5_ps <- phyloseq(OTU, TAX, meta)

#add reference sequence and replace variants with ASVs
dna_V4V5 <- Biostrings::DNAStringSet(taxa_names(V4V5_ps))
names(dna_V4V5) <- taxa_names(V4V5_ps)
V4V5_ps <- merge_phyloseq(V4V5_ps, dna_V4V5)
taxa_names(V4V5_ps) <- paste0("V4V5_ASV", seq(ntaxa(V4V5_ps)))

#####################################
#Generate an overview table
#####################################
summary_table <- data.frame()

for (p in c("V3V4","V4V5")){
i <- paste(p, "_ps", sep = "")
  
phy_obj <- get(i)
#metadata
meta <- as(sample_data(phy_obj),"data.frame")
meta$Row.names<-sample_names(sample_data(phy_obj))

#dada2 workflow
reads.tab <- read.csv2(file.path("data",p,"libs_summary_table.txt"),header = TRUE, sep = "\t")
reads.tab$Row.names <-paste("X",row.names(reads.tab), sep ="")

#check how many ASVs were unclassified on phylum level, or assigned to chloroplast and Mitochondria
ps_euk <- as(sample_sums(subset_taxa(phy_obj, Kingdom %in% c("Eukaryota"))),"vector")
ps_uncl <- as(sample_sums(subset_taxa(phy_obj, Phylum %in% c("Bacteria_uncl","Archaea_uncl","NA_uncl"))),"vector")
ps_chl <- as(sample_sums(subset_taxa(phy_obj, Order %in% c("Chloroplast"))),"vector")
ps_mit <- as(sample_sums(subset_taxa(phy_obj, Family %in% c("Mitochondria"))),"vector")

pruned_seq_sums <- data.frame(Eukaryota = ps_euk, Phyl_uncl = ps_uncl,
                              Chloroplast = ps_chl,Mitochondria = ps_mit)

pruned_seq_sums$Row.names <- rownames(pruned_seq_sums)

#remove unclassified on phylum level, chloroplast and Mitochondrial sequence variants
phy_obj0 <- subset_taxa(phy_obj, !Kingdom %in% c("Eukaryota") &!Phylum %in% c("Bacteria_uncl","Archaea_uncl","NA_uncl") & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") )

#alpha diversity indexes
ps0_alpha <- estimate_richness(phy_obj0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
ps0_alpha$Row.names<-rownames(ps0_alpha)

#merge all together
ps0_summary_table <- merge(meta,reads.tab,by ="Row.names") %>%
  merge(pruned_seq_sums,by ="Row.names")%>%
  merge(ps0_alpha,by ="Row.names")%>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(Seq.prop = round(tabled/input,2),
         Euk.Seq.prop = round(Eukaryota/tabled,2),
         Phyl_uncl.Seq.prop = round(Phyl_uncl/tabled,2),
         Chloroplast.Seq.prop = round(Chloroplast/tabled,2),
         Mitochondria.Seq.prop = round(Mitochondria/tabled,2)) %>%
  select("sample_title","Primer_set","Type","Depth", #metadata
         "input","merged", "tabled","Seq.prop", #dada2 
         "Euk.Seq.prop", "Phyl_uncl.Seq.prop","Chloroplast.Seq.prop","Mitochondria.Seq.prop", #taxa
         "Observed","Chao1","Shannon","InvSimpson")   #alpha div


summary_table <- rbind(summary_table,ps0_summary_table)

#save the clean phyloseq object to a proper variable
assign(paste(i,"0", sep =""),phy_obj0)
}

write.table(summary_table, "./Tables/Micro_overview_table.txt" , sep = "\t", quote = F)

#####################################
#Test patterns in Alpha diversity
#####################################
#not normally distributed, therefore applying Kruskal Wallis test
kruskal.test(Chao1 ~ Primer_set, data = summary_table)

#posthoc wilcoxon test to compared the different depths
Chao1_Wilcox<- summary_table  %>%
  group_by(Type) %>%
  wilcox_test(Chao1 ~ Primer_set, p.adjust.method = "bonferroni") %>%
  add_significance()

mean_chao1<- summary_table %>%
              filter(Type %in% c("Deep water", "Sediment")) %>%
  group_by(Type, Primer_set)%>%
  summarise(Chao1_mean = mean(Chao1))

#####################################
#compare richness of Bacteria and Archaea separately in Deep water
#####################################
#V3V4
V3V4_ps_bac <- subset_samples(V3V4_ps, Type %in% c("Deep water"))
V3V4_ps_bac <- subset_taxa(V3V4_ps_bac, Kingdom %in% c("Bacteria"))
V3V4_ps_bac<- prune_taxa(taxa_sums(V3V4_ps_bac)>0,V3V4_ps_bac)

V3V4_ps_Sed <- subset_samples(V3V4_ps, Type %in% c("Sediment"))
V3V4_ps_Sed <- subset_taxa(V3V4_ps_Sed, Kingdom %in% c("Bacteria"))
V3V4_ps_Sed<- prune_taxa(taxa_sums(V3V4_ps_Sed)>0,V3V4_ps_Sed)
#V4V5
V4V5_ps_bac <- subset_samples(V4V5_ps, Type %in% c("Deep water"))
V4V5_ps_bac <- subset_taxa(V4V5_ps_bac, Kingdom %in% c("Bacteria"))
V4V5_ps_bac<- prune_taxa(taxa_sums(V4V5_ps_bac)>0,V4V5_ps_bac)

V4V5_ps_Sed <- subset_samples(V4V5_ps, Type %in% c("Sediment"))
V4V5_ps_Sed <- subset_taxa(V4V5_ps_Sed, Kingdom %in% c("Bacteria"))
V4V5_ps_Sed<- prune_taxa(taxa_sums(V4V5_ps_Sed)>0,V4V5_ps_Sed)


#####################################
#Plot rarefaction
#####################################
iNEXT.rare.line <- data.frame()
iNEXT.rare.point<- data.frame()

for (p in c("V3V4","V4V5")){
  i <- paste(p, "_ps0", sep = "")
  phy_obj <- get(i)
  
  ps0.iNEXT <- iNEXT(as.data.frame(otu_table(phy_obj)), q=0, datatype="abundance", knots = 40)
  
  ps0.iNEXT.rare <-fortify(ps0.iNEXT, type=1)
  ps0.iNEXT.meta <- as(sample_data(phy_obj), "data.frame")
  ps0.iNEXT.meta$site <- rownames(ps0.iNEXT.meta)
  
  ps0.iNEXT.rare$Type <- ps0.iNEXT.meta$Type[match(ps0.iNEXT.rare$site, ps0.iNEXT.meta$site)] 
  ps0.iNEXT.rare$SampleName <- ps0.iNEXT.meta$sample_title[match(ps0.iNEXT.rare$site, ps0.iNEXT.meta$site)] 
  
  ps0.iNEXT.rare.point <- ps0.iNEXT.rare[which(ps0.iNEXT.rare$method == "observed"),]
  ps0.iNEXT.rare.line <- ps0.iNEXT.rare[which(ps0.iNEXT.rare$method != "observed"),]
  ps0.iNEXT.rare.line$method <- factor (ps0.iNEXT.rare.line$method,
                                        c("interpolated", "extrapolated"),
                                        c("interpolation", "extrapolation"))
  
  ps0.iNEXT.rare.line$Primer <- p
  ps0.iNEXT.rare.point$Primer <- p
  
  iNEXT.rare.line <- rbind(iNEXT.rare.line,ps0.iNEXT.rare.line)
  iNEXT.rare.point <- rbind(iNEXT.rare.point,ps0.iNEXT.rare.point)
}


iNEXT.rare.line$Type <- factor(iNEXT.rare.line$Type, levels = c("Sea ice", "Surface water", "Deep water", "Sediment trap","Sediment"))
iNEXT.rare.point$Type <- factor(iNEXT.rare.point$Type, levels = c("Sea ice", "Surface water", "Deep water", "Sediment trap","Sediment"))

rare.p <- ggplot(iNEXT.rare.line, aes(x=x, y=y, shape = site))+
  geom_line(aes(linetype = method, colour = Primer), lwd = 0.5, data= iNEXT.rare.line)+
  geom_point(shape = 21, size =3, colour = "black", data= iNEXT.rare.point)+
  labs(x = "Sample size", y = "Species richness")+
  xlim(0,3e5)+
  theme_classic(base_size = 12)+theme(legend.position="none")+
  facet_grid(Type~Primer, scales = "free",as.table = TRUE)+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) + 
  coord_cartesian(clip="off")


ggsave("Figures/rarefactions.pdf", rare.p)




#####################################
# Plot community composition
#####################################
#remove singletons
#V3V4
#in how many samples did each taxa appear at least once
prev0_V3V4 <- apply(X = otu_table(V3V4_ps0),
              MARGIN = ifelse(taxa_are_rows(V3V4_ps0), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf_V3V4 <- data.frame(Prevalence = prev0_V3V4,
                    TotalAbundance = taxa_sums(V3V4_ps0),
                    tax_table(V3V4_ps0),
                    Primer = "V3V4")

prevdf_V3V4<- prevdf_V3V4 %>% filter (Prevalence>1)
V3V4_ps0.prev <-  prune_taxa((prev0_V3V4 > 1), V3V4_ps0)

#V4V5
#in how many samples did each taxa appear at least once
prev0_V4V5 <- apply(X = otu_table(V4V5_ps0),
                    MARGIN = ifelse(taxa_are_rows(V4V5_ps0), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf_V4V5 <- data.frame(Prevalence = prev0_V4V5,
                          TotalAbundance = taxa_sums(V4V5_ps0),
                          tax_table(V4V5_ps0),
                          Primer = "V4V5")

prevdf_V4V5<- prevdf_V4V5 %>% filter (Prevalence>1)
V4V5_ps0.prev <-  prune_taxa((prev0_V4V5 > 1), V4V5_ps0)


bar_plot <- list()

for (p in c("V3V4","V4V5")){
  i <- paste(p, "_ps0.prev", sep = "")
  phy_obj <- get(i)

  #transform data
  phy_obj.ra <- transform_sample_counts(phy_obj, function(x) x / sum(x))
  
  #melt phyloseq object
  phy_obj.ra.long <- psmelt(phy_obj.ra)
  phy_obj.ra.long$Abundance <- phy_obj.ra.long$Abundance*100
  
  #fix unclassified lineages 
  phy_obj.ra.long$Class <- as.character(phy_obj.ra.long$Class)
  phy_obj.ra.long$Class[is.na(phy_obj.ra.long$Class)] <- paste(phy_obj.ra.long$Phylum[is.na(phy_obj.ra.long$Class)],"uc", sep = "_")

  #calculate abundance for each Class
  phy_obj.ra.long.agg <- phy_obj.ra.long %>% 
    select(Overlap_ID,Primer_set,Type,Class,Abundance)%>%
    group_by(Overlap_ID,Primer_set,Type,Class) %>%
    dplyr::summarise(Abund.total= sum(Abundance)) 
  
  #remove below 1% ra
  taxa_classes <- unique(phy_obj.ra.long.agg$Class)
  phy_obj.ra.long.agg$Class[phy_obj.ra.long.agg$Abund.total<2] <- "Other taxa <2%"
  phy_obj.ra.long.agg$Class <- factor(phy_obj.ra.long.agg$Class,
                                     levels=c(taxa_classes,"Other taxa <2%"))

  phy_obj.ra.long.agg$Type<- factor(phy_obj.ra.long.agg$Type,
                                               levels= c("Sea ice", "Surface water","Sediment trap", "Deep water", "Sediment"))
  
  bar_plot[[p]] <-ggplot(phy_obj.ra.long.agg, aes(x = Primer_set, y = Abund.total, fill = Class))+
    facet_grid(Type~Primer_set, space= "fixed")+
    geom_bar(stat = "identity", position="fill") +
    scale_y_continuous(breaks=NULL)+
    scale_fill_manual(values = phyla.col)+ 
    coord_polar("y")
}

ggarrange(bar_plot$V3V4, bar_plot$V4V5, #heights = c(2,1.2),
          ncol = 2, nrow = 1, align = "hv", legend = "left",
          legend.grob = do.call(rbind, c(list(get_legend(bar_plot["V3V4"]),get_legend(bar_plot["V4V5"])), size="first")))


ggsave("./figures/pie_charts.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
# Explore ASVs prevalence on phylum level
#####################################
#compare on Phylum level
prevdf_V3V4_sum_phylum <- plyr::ddply(prevdf_V3V4, "Phylum", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

#define top phyla
top_Phylum_V3V4 <- prevdf_V3V4_sum_phylum %>% 
  top_n(10, Proportion)


#explore datasets on phylum levels
prevdf_V4V5_sum_phylum <- plyr::ddply(prevdf_V4V5, "Phylum", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

#define top phyla
top_Phylum_V4V5 <- prevdf_V4V5_sum_phylum %>% 
  top_n(10, Proportion)


#merge the top phyla from both datasets
top_Phylum <- union(top_Phylum_V3V4$Phylum, top_Phylum_V4V5$Phylum)

#remove singletons
prevdf_V3V4_top <- prevdf_V3V4 %>%
  filter(Phylum %in% top_Phylum)


#remove singletons 
prevdf_V4V5_top <- prevdf_V4V5 %>%
  filter(Phylum %in% top_Phylum)

prevdf_both_top <- rbind(prevdf_V3V4_top,prevdf_V4V5_top)

# #plot
prev_plot_phyl <- ggplot(prevdf_both_top, aes(TotalAbundance, Prevalence,color=Primer)) +
  # # Include a guess for parameter
  geom_point(size = 2, alpha = 0.3) +
  geom_smooth(method  = "loess")+
  scale_x_log10() +  xlab("Total Abundance") + ylab("Number of samples") +
  facet_wrap(~Phylum) + theme(legend.position="bottom")


ggsave("./figures/Phyla_prevalance.pdf", 
       plot = prev_plot_phyl,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
# Explore ASVs distribution on class level
#####################################
#define top classes
prevdf_V3V4_sum<- prevdf_V3V4 %>%
  plyr::ddply( "Class", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

prevdf_V4V5_sum<- prevdf_V4V5 %>%
  plyr::ddply( "Class", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

#define top classes
top_class_V3V4 <- prevdf_V3V4_sum %>% 
  filter(Proportion > 1)

top_class_V4V5 <- prevdf_V4V5_sum %>% 
  filter(Proportion > 1)

top_class <- union(top_class_V3V4$Class,top_class_V4V5$Class)

# merge the datasets and plot ASV abundance
prevdf_both <- rbind(prevdf_V3V4,prevdf_V4V5)

class_dist.p<- ggplot(prevdf_both[prevdf_both$Class %in% top_class,], aes(Primer, TotalAbundance, fill=Primer))+
  geom_jitter(color = "gray50", alpha =.1)+
  geom_boxplot()+
  geom_signif(comparisons = list(c("V3V4", "V4V5")),
             map_signif_level=TRUE, test = "wilcox.test", color = "black")+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, 
               position = position_dodge(width = 1))+
  scale_y_log10()+
  facet_wrap(~Class, ncol = 3)+
  theme(legend.position="none",
    strip.background.x = element_rect(
      color=NULL, fill=NULL, size=0, linetype=NULL))

ggsave("Figures/class_ASVS_ dist.pdf", 
       plot = class_dist.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#explore the differences inside each Class
V3V4_total_seqs<- sum(sample_sums(V3V4_ps0.prev))
V4V5_total_seqs<- sum(sample_sums(V4V5_ps0.prev))

prevdf_both_diff_Class <-  prevdf_both %>% 
  filter(Class %in% top_class) %>% 
  mutate(TotalAbundance = case_when(Primer =="V3V4" ~ TotalAbundance/V3V4_total_seqs,
                                    Primer =="V4V5" ~ TotalAbundance/V4V5_total_seqs)) %>% 
  group_by(Primer,Class) %>%
  dplyr::summarise(Seq.prop= round(sum(TotalAbundance),2)*100, ASVs.num = length(TotalAbundance)) %>% 
  tidyr::pivot_wider(id_cols = Class, 
                     names_from = Primer, 
                     values_from = c(Seq.prop,ASVs.num)) %>% 
  mutate(ASVs_prop = ASVs.num_V3V4/ASVs.num_V4V5)

prevdf_both_diff_order <-  prevdf_both %>% 
  filter(Class %in% top_class) %>% 
  mutate(TotalAbundance = case_when(Primer =="V3V4" ~ TotalAbundance/V3V4_total_seqs,
                                    Primer =="V4V5" ~ TotalAbundance/V4V5_total_seqs)) %>% 
  group_by(Primer,Class,Order) %>%
  dplyr::summarise(Seq.prop= round(sum(TotalAbundance),2)*100, ASVs.num = length(TotalAbundance)) %>% 
  tidyr::pivot_wider(id_cols = c(Class,Order), 
              names_from = Primer, 
              values_from = c(Seq.prop,ASVs.num)) %>% 
  mutate(ASVs_prop = ASVs.num_V3V4/ASVs.num_V4V5)


prevdf_both_diff_Family <-  prevdf_both %>% 
  filter(Class %in% top_class) %>% 
  mutate(TotalAbundance = case_when(Primer =="V3V4" ~ TotalAbundance/V3V4_total_seqs,
                                    Primer =="V4V5" ~ TotalAbundance/V4V5_total_seqs)) %>% 
  group_by(Primer,Class,Family) %>%
  dplyr::summarise(Seq.prop= round(sum(TotalAbundance),2)*100, ASVs.num = length(TotalAbundance)) %>% 
  tidyr::pivot_wider(id_cols = c(Class,Family), 
                     names_from = Primer, 
                     values_from = c(Seq.prop,ASVs.num))  %>% 
  mutate(ASVs_prop = ASVs.num_V3V4/ASVs.num_V4V5)


prevdf_both_diff_Genus <-  prevdf_both %>% 
  filter(Class %in% top_class) %>% 
  mutate(TotalAbundance = case_when(Primer =="V3V4" ~ TotalAbundance/V3V4_total_seqs,
                                    Primer =="V4V5" ~ TotalAbundance/V4V5_total_seqs)) %>% 
  group_by(Primer,Class,Genus) %>%
  dplyr::summarise(Seq.prop= round(sum(TotalAbundance),2)*100, ASVs.num = length(TotalAbundance)) %>% 
  tidyr::pivot_wider(id_cols = c(Class,Genus), 
                     names_from = Primer, 
                     values_from = c(Seq.prop,ASVs.num))  %>% 
  mutate(ASVs_prop = ASVs.num_V3V4/ASVs.num_V4V5)

write.csv(prevdf_both_diff_Family, "Tables/Family_ASVs.csv")

#####################################
# Merge datasets on a Genus level
#####################################
#summarize genera in V3V4 dataset
V3V4_tax_table.agg <- as.data.frame(as.list(aggregate(Species~Class+Order+Family+Genus, 
                                                      as.data.frame(tax_table(V3V4_ps0.prev)),
                                                  FUN = function(x) c(Number = length(x)))))

V3V4_ps0_glom <- tax_glom(V3V4_ps0.prev, taxrank="Genus")

#correct duplicated name of Unknown family in both Alpha- and Gamma- proteobacteria
tax_table(V3V4_ps0_glom)[,"Genus"][tax_table(V3V4_ps0_glom)[,"Order"] =="Alphaproteobacteria_Incertae_Sedis"]<- "Alpha_unk_family"
taxa_names(V3V4_ps0_glom) <- tax_table(V3V4_ps0_glom)[,"Genus"]
sample_names(V3V4_ps0_glom)<- paste(sample_names(V3V4_ps0_glom),"V3V4", sep ="_")


#summarize genera in V4V5 dataset
V4V5_tax_table.agg <- as.data.frame(as.list(aggregate(Species~Class+Order+Family+Genus, 
                                                      as.data.frame(tax_table(V4V5_ps0.prev)),
                                                      FUN = function(x) c(Number = length(x)))))

V4V5_ps0_glom <- tax_glom(V4V5_ps0.prev, taxrank="Genus")
#correct duplicated name of Unknown family in both Alpha- and Gamma- proteobacteria
tax_table(V4V5_ps0_glom)[,"Genus"][tax_table(V4V5_ps0_glom)[,"Order"] =="Alphaproteobacteria_Incertae_Sedis"]<- "Alpha_unk_family"
taxa_names(V4V5_ps0_glom) <- tax_table(V4V5_ps0_glom)[,"Genus"]
sample_names(V4V5_ps0_glom)<- paste(sample_names(V4V5_ps0_glom),"V4V5", sep ="_")

#merge both datasets
ps0_glom_merged <- merge_phyloseq(V3V4_ps0_glom, V4V5_ps0_glom)

#calculate sequence proportion and melt the dataframe
ps0_glom_merged.prop <- transform_sample_counts(ps0_glom_merged, function(otu) otu/sum(otu))

ps0_glom_merged.melt <- psmelt(ps0_glom_merged.prop)

#Valvulate overlap of lineages in the sequence abundance of the unique ones
genus_over <- list()
genus_over$V3V4 <- as.vector(tax_table(V3V4_ps0_glom)[,"Genus"])
genus_over$V4V5 <- as.vector(tax_table(V4V5_ps0_glom)[,"Genus"])

genus_over.venn <- venn.diagram(genus_over,
                              lwd = 1, fill = c("red", "yellow"),  alpha = c(0.5, 0.5), cex = 2,
                              #cat.fontface = 4,
                              lty =2, filename = NULL, scaled = TRUE,
                              inverted = FALSE, print.mode = c("raw","percent"))

grid.newpage()
grid.draw(genus_over.venn)

#identify unique genera in each dataset
V3V4_unique_genera <- taxa_names(V3V4_ps0_glom)[!taxa_names(V3V4_ps0_glom) %in% taxa_names(V4V5_ps0_glom)]
V4V5_unique_genera <- taxa_names(V4V5_ps0_glom)[!taxa_names(V4V5_ps0_glom) %in% taxa_names(V3V4_ps0_glom)]

#calculate abundance for each taxa
ps0_glom_Genus <- as.data.frame(as.list(aggregate(Abundance~Overlap_ID+Primer_set+Type+Class+Order+Family+Genus, ps0_glom_merged.melt,
                                                  FUN = function(x) c(sum = sum(x), count=length(x)))))
ps0_glom_Genus$Abundance.sum <- ps0_glom_Genus$Abundance.sum*100
ps0_glom_Genus<- ps0_glom_Genus[ps0_glom_Genus$Abundance.sum>0,]

ps0_glom_Genus_by_primer <- as.data.frame(as.list(aggregate(Abundance.sum~Primer_set+Type+Class+Order+Family+Genus, ps0_glom_Genus,
                                                            FUN = function(x) c(mean = mean(x), sd=sd(x)))))


V3V4_unique_glom_Genus <-ps0_glom_Genus_by_primer %>% filter (Genus %in% V3V4_unique_genera)


V4V5_unique_glom_Genus <-ps0_glom_Genus_by_primer %>% filter (Genus %in% V4V5_unique_genera)


#############################
# Dissimilarity
#############################
#statistical significance of the groups
df <- as(sample_data(ps0_glom_merged.prop), "data.frame")
d <- phyloseq::distance(ps0_glom_merged.prop, "jsd")
adonis_all <- adonis2(d ~ Type+Primer_set , data= df, perm = 999)
adonis_all


# plot
ps0.prop.jsd_slv <- ordinate(ps0_glom_merged.prop, method = "NMDS", 
                             distance = "jsd")
stressplot(ps0.prop.jsd_slv)

plot_ordination(ps0_glom_merged.prop, ps0.prop.jsd_slv, color = "Type")+
  stat_ellipse(aes( group = Type), type="norm", color = "black", alpha = 0.5)+
  geom_point(aes(shape = Primer_set),color ="black", size = 6)+
  geom_point(aes(shape = Primer_set), size = 4)+
  scale_color_manual(values = c("Sea ice" = "white", "Surface water" = "green",
                                "Deep water" = "darkblue", "Sediment trap" = "grey",
                                "Sediment" = "brown"))+
  #geom_text(aes( label = Overlap_ID), color ="black", size = 3)+
  coord_fixed()


#############################
# Abundance comparison to CARD-FISH
#############################
#import raw data
raw.counts <- read.csv("./data/CARD-FISH/card-fish-counts-SH.csv", dec = ".", stringsAsFactors = FALSE)
raw.counts<- raw.counts[raw.counts$conc.FISH>0,]
raw.counts$conc.FISH <- as.numeric(raw.counts$conc.FISH)

#microscope calculation factor
calc.factor <- 99515.5458411807

#Split sample name and remove SV2 station and ABYSS depth and calculate counts per FOV
cell.counts<- raw.counts %>% 
  mutate(DAPI_conc= (DAPI_Nr_Set*calc.factor)/Volume,
         FISH_conc = (DAPI_Nr_SubSet*calc.factor)/Volume)


#summarize counts
counts.agg.FISH <- as.data.frame(as.list(aggregate(conc.FISH~Overlap_ID+Type+Taxa,raw.counts,FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))))

counts.agg.DAPI <- as.data.frame(as.list(aggregate(conc.DAPI~Overlap_ID+Type+Taxa,raw.counts,FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))))

counts.agg.Station <- merge(counts.agg.FISH, counts.agg.DAPI) %>%
  group_by(Type,Overlap_ID)%>%
  filter(Type %in% c("Sea ice","Surface water"),
         !Taxa %in% c("EUB"))%>%
  mutate(DAPI_mean = mean(conc.DAPI.mean),
         Abundance.sum = 100*conc.FISH.mean/DAPI_mean,
         Primer_set ="FISH")



counts_prop<- raw.counts %>%
  mutate(Proportion = conc.FISH/conc.DAPI)%>%
  group_by(Type, Taxa)%>%
  summarise(Proportion.mean= mean(Proportion),
            Proportion.sd= sd(Proportion))


#subset only samples that were used for CARD-FISH
ps0_glom_merged.FISH<- ps0_glom_merged.melt[ps0_glom_merged.melt$FISH=="yes",]
ps0_glom_merged.FISH<- ps0_glom_merged.FISH[ps0_glom_merged.FISH$Type %in% c("Sea ice","Surface water"), ]

#calculate abundance for each Genus
ps0_glom_Genus <- as.data.frame(as.list(aggregate(Abundance~Overlap_ID+Primer_set+Type+Class+Order+Family+Genus, ps0_glom_merged.FISH,
                                                   FUN = function(x) c(sum = sum(x), count=length(x)))))
ps0_glom_Genus$Abundance.sum <- ps0_glom_Genus$Abundance.sum*100
ps0_glom_Genus<- ps0_glom_Genus[ps0_glom_Genus$Abundance.sum>0,]
ps0_glom_Genus$Taxa<- ps0_glom_Genus$Genus

#calculate abundance for each Order
ps0_glom_Order <- as.data.frame(as.list(aggregate(Abundance~Overlap_ID+Primer_set+Type+Class+Order, ps0_glom_merged.FISH,
                                                   FUN = function(x) c(sum = sum(x), count=length(x)))))
ps0_glom_Order$Abundance.sum <- ps0_glom_Order$Abundance.sum*100
ps0_glom_Order<- ps0_glom_Order[ps0_glom_Order$Abundance.sum>0,]
ps0_glom_Order$Taxa<- ps0_glom_Order$Order

#calculate abundance for each Class
ps0_glom_Class <- as.data.frame(as.list(aggregate(Abundance~Overlap_ID+Primer_set+Type+Class, ps0_glom_merged.FISH,
                                                  FUN = function(x) c(sum = sum(x), count=length(x)))))
ps0_glom_Class$Abundance.sum <- ps0_glom_Class$Abundance.sum*100
ps0_glom_Class<- ps0_glom_Class[ps0_glom_Class$Abundance.sum>0,]
ps0_glom_Class$Taxa<- ps0_glom_Class$Class



seq_prop_selected_taxa <- rbind(ps0_glom_Class[ps0_glom_Class$Taxa %in% c("Gammaproteobacteria", "Bacteroidia"),c("Overlap_ID", "Primer_set", "Type","Taxa","Abundance.sum")],
              ps0_glom_Order[ps0_glom_Order$Taxa %in% c("Alteromonadales", "SAR11_clade"),c("Overlap_ID", "Primer_set", "Type","Taxa","Abundance.sum")],
              ps0_glom_Genus[ps0_glom_Genus$Taxa %in% c("Polaribacter"),c("Overlap_ID", "Primer_set", "Type","Taxa","Abundance.sum")],
              counts.agg.Station[,c("Overlap_ID", "Primer_set","Type", "Taxa", "Abundance.sum") ])



seq_prop_selected_taxa_Wilcox <- seq_prop_selected_taxa   %>%
  rstatix::wilcox_test(Abundance.sum ~ Primer_set, p.adjust.method = "BH",paired = TRUE) %>%
  add_significance()


seq_prop_selected_taxa_wide <- seq_prop_selected_taxa %>% 
                    tidyr::pivot_wider(id_cols = Overlap_ID, 
                   names_from = c(Taxa,Primer_set), 
                   values_from = Abundance.sum)

