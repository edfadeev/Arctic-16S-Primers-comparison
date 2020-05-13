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
set.seed(1586)
theme_set(theme_bw())

#function to add total number of observations to a ggboxplot
give.n <- function(x){
  return(c(y = median(x)*1.15, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

#define colour palletes
phyla.col <- c("Acidimicrobiia"="#AA4488",
               "Actinobacteria" = "#DD1232",
               "Alphaproteobacteria"= "#771155",
               "Bacilli"="#117744", 
               "Bacteroidia"= "#77AADD",
               "Campylobacteria" = "#CC1234",
               "Chlamydiae"= "#DDAA77" ,
               "Clostridia"= "#77CCCC", 
               "Dehalococcoidia"= "#AAAA44",
               "Desulfobacteria" ="#FF2222", 
               "Desulfobulbia"= "#117744", 
               "Gammaproteobacteria"="#117777",
               "Kiritimatiellae"= "#44AA77",
               "Lentisphaeria" = "#DD7788",
               "Marinimicrobia_(SAR406_clade)_uncl"= "#34ABAA", 
               "Methanosarcinia"= "#11BBCC", 
               "NB1-j_uncl" = "#774411",
               "Nitrososphaeria" = "#E69F00",
               "Phycisphaerae"= "#88CCAA", 
               "Planctomycetes"= "#777711",
               "OM190"= "#009E73",
               "SAR324_clade(Marine_group_B)_uncl"="#CC99BB",
               "Spirochaetia" = "#AACC45",
               "Thermoplasmata" = "#0072B2",
               "Verrucomicrobiae" = "#AA7744",
               "Vicinamibacteria" ="#DDDD77",
               "Other taxa"= "#114477")

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
taxa_names(V3V4_ps) <- paste0("ASV", seq(ntaxa(V3V4_ps)))

#remove the unsuccesful sample 9
V3V4_ps<- subset_samples(V3V4_ps, !sample_title %in% c("FRAM_2963_Fevi_32_oben_2",
                                                       "FRAM_1078_Fevi_3_oben_1",
                                                       "FRAM_1098_Fevi_3_oben_6",
                                                       "FRAM_1122_Fevi_3_oben_12",
                                                       "FRAM_1150_Fevi_3_oben_19"))

#remove unobserved ASVs
V3V4_ps<- prune_taxa(taxa_sums(V3V4_ps)>0, V3V4_ps)

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

V4V5_ps<- subset_samples(V4V5_ps, !sample_title %in% c("FRAM_2963_Fevi_32_oben_2",
                                                       "FRAM_1078_Fevi_3_oben_1",
                                                       "FRAM_1098_Fevi_3_oben_6",
                                                       "FRAM_1122_Fevi_3_oben_12",
                                                       "FRAM_1150_Fevi_3_oben_19"))

#remove unobserved ASVs
V4V5_ps<- prune_taxa(taxa_sums(V4V5_ps)>0, V4V5_ps)


#add reference sequence and replace variants with ASVs
dna_V4V5 <- Biostrings::DNAStringSet(taxa_names(V4V5_ps))
names(dna_V4V5) <- taxa_names(V4V5_ps)
V4V5_ps <- merge_phyloseq(V4V5_ps, dna_V4V5)
taxa_names(V4V5_ps) <- paste0("ASV", seq(ntaxa(V4V5_ps)))


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

#alpha diversity indeces
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
         "Observed","Chao1","Shannon","InvSimpson") #alpha div

summary_table <- rbind(summary_table,ps0_summary_table)

#save the clean phyloseq object to a proper variable
assign(paste(i,"0", sep =""),phy_obj0)
}

write.table(summary_table, "./Tables/Micro_overview_table.txt" , sep = "\t", quote = F)


#Statistical comparison by type of samples
Chao1_Wilcox <- summary_table   %>%
  group_by(Type) %>%
  rstatix::wilcox_test(Chao1 ~ Primer_set, p.adjust.method = "BH",paired = TRUE) %>%
  add_significance()

Shannon_Wilcox <- summary_table   %>%
  group_by(Type) %>%
  rstatix::wilcox_test(Shannon ~ Primer_set, p.adjust.method = "BH",paired = TRUE) %>%
  add_significance()

InvSimpson_Wilcox <- summary_table   %>%
  group_by(Type) %>%
  rstatix::wilcox_test(InvSimpson ~ Primer_set, p.adjust.method = "BH",paired = TRUE) %>%
  add_significance()

#####################################
#Plot rarefaction
####################################
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
# ASVs distribution
#####################################
#in how many samples did each taxa appear at least once
prev0_V3V4 <- apply(X = otu_table(V3V4_ps0),
              MARGIN = ifelse(taxa_are_rows(V3V4_ps0), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf_V3V4 <- data.frame(Prevalence = prev0_V3V4,
                    TotalAbundance = taxa_sums(V3V4_ps0),
                    tax_table(V3V4_ps0),
                    Primer = "V3V4")

#remove singletons
prevdf_V3V4<- prevdf_V3V4 %>% filter (Prevalence>1)

#compare on Phylum level
prevdf_V3V4_sum_phylum <- plyr::ddply(prevdf_V3V4, "Phylum", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))


#define top phyla
top_Phylum_V3V4 <- prevdf_V3V4_sum_phylum %>% 
  top_n(10, Proportion)


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

#remove singletons
prevdf_V4V5<- prevdf_V4V5 %>% filter (Prevalence>1)

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
  # geom_smooth(method  = "lm")+
  scale_x_log10() +  xlab("Total Abundance") + ylab("Number of samples") +
  facet_wrap(~Phylum) + theme(legend.position="none")

prev_plot_phyl

#explore results on class level
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
  filter(Proportion > 2)

top_class_V4V5 <- prevdf_V4V5_sum %>% 
  filter(Proportion > 2)

top_class <- union(top_class_V3V4$Class,top_class_V4V5$Class)

# merge and extract ASVs that were observed in at least two samples
prevdf_both <- rbind(prevdf_V3V4,prevdf_V4V5)

#box plot of ASVs abundance
ggplot(prevdf_both[prevdf_both$Class %in% top_class,], aes(Primer, TotalAbundance, fill=Primer))+
  geom_jitter(color = "gray50", alpha =.1)+
  geom_boxplot()+
  geom_signif(comparisons = list(c("V3V4", "V4V5")),
             map_signif_level=TRUE, test = "wilcox.test", color = "black")+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, 
               position = position_dodge(width = 1))+
  scale_y_log10()+
  facet_wrap(~Class, ncol = 3)+
  theme_classic(base_size = 12)+theme(legend.position="none")

# Remove singletons
V3V4_ps0.prev <-  prune_taxa((prev0_V3V4 > 1), V3V4_ps0)
V4V5_ps0.prev <-  prune_taxa((prev0_V4V5 > 1), V4V5_ps0)

#####################################
# Dissimilarities between the communities on Genus level
#####################################

#summarize genera in V3V4 dataset
V3V4_tax_table <- as.data.frame(as(tax_table(V3V4_ps0.prev),"matrix"))

V3V4_tax_table.agg <- as.data.frame(as.list(aggregate(Species~Class+Order+Family+Genus, 
                                                      V3V4_tax_table,
                                                  FUN = function(x) c(Number = length(x)))))

V3V4_ps0_glom <- tax_glom(V3V4_ps0.prev, taxrank="Genus")
#correct duplicated name of Unknown family in both Alpha- and Gamma- proteobacteria
tax_table(V3V4_ps0_glom)[,"Genus"][tax_table(V3V4_ps0_glom)[,"Order"] =="Alphaproteobacteria_Incertae_Sedis"]<- "Alpha_unk_family"
taxa_names(V3V4_ps0_glom) <- tax_table(V3V4_ps0_glom)[,"Genus"]
sample_names(V3V4_ps0_glom)<- paste(sample_names(V3V4_ps0_glom),"V3V4", sep ="_")


#summarize genera in V4V5 dataset
V4V5_tax_table <- as.data.frame(as(tax_table(V4V5_ps0.prev),"matrix"))

V4V5_tax_table.agg <- as.data.frame(as.list(aggregate(Species~Class+Order+Family+Genus, 
                                                      V4V5_tax_table,
                                                      FUN = function(x) c(Number = length(x)))))

V4V5_ps0_glom <- tax_glom(V4V5_ps0.prev, taxrank="Genus")
#correct duplicated name of Unknown family in both Alpha- and Gamma- proteobacteria
tax_table(V4V5_ps0_glom)[,"Genus"][tax_table(V4V5_ps0_glom)[,"Order"] =="Alphaproteobacteria_Incertae_Sedis"]<- "Alpha_unk_family"
taxa_names(V4V5_ps0_glom) <- tax_table(V4V5_ps0_glom)[,"Genus"]
sample_names(V4V5_ps0_glom)<- paste(sample_names(V4V5_ps0_glom),"V4V5", sep ="_")


#identify unique genera in each dataset
V3V4_unique_genera <- taxa_names(V3V4_ps0_glom)[!taxa_names(V3V4_ps0_glom) %in% taxa_names(V4V5_ps0_glom)]
V4V5_unique_genera <- taxa_names(V4V5_ps0_glom)[!taxa_names(V4V5_ps0_glom) %in% taxa_names(V3V4_ps0_glom)]

#merge both datasets
ps0_glom_merged <- merge_phyloseq(V3V4_ps0_glom, V4V5_ps0_glom)

#disimilarity
ps0_glom_merged.prop <- transform_sample_counts(ps0_glom_merged, function(otu) otu/sum(otu))

#melt data
ps0_glom_merged.melt <- psmelt(ps0_glom_merged.prop)

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
# Taxonomic compositions on class level
#############################
#calculate abundance for each taxa
ps0_glom_merged.melt.agg <- aggregate(Abundance~Overlap_ID+Primer_set+Type+Class, ps0_glom_merged.melt,
                                      FUN = "sum")

ps0_glom_merged.melt.agg$Abundance <- ps0_glom_merged.melt.agg$Abundance*100
ps0_glom_merged.melt.agg<- ps0_glom_merged.melt.agg[ps0_glom_merged.melt.agg$Abundance>0,]

#remove below 3% ra
threshold<- 3
ps0_glom_merged.melt.agg$Class <- as.character(ps0_glom_merged.melt.agg$Class)

taxa_classes <- unique(ps0_glom_merged.melt.agg$Class[!ps0_glom_merged.melt.agg$Abundance<threshold])

ps0_glom_merged.melt.agg$Class[ps0_glom_merged.melt.agg$Abundance<threshold] <- "Other taxa"

ps0_glom_merged.melt.agg$Class <- factor(ps0_glom_merged.melt.agg$Class,
                                         levels=c(taxa_classes,"Other taxa"))


ps0_glom_merged.melt.agg$Type <- factor(ps0_glom_merged.melt.agg$Type,
                                        levels= c("Sea ice","Surface water","Sediment trap", "Deep water","Sediment"))

p <- list()
for (i in c("V3V4","V4V5")){
  sub <- ps0_glom_merged.melt.agg[ps0_glom_merged.melt.agg$Primer_set == i,]
  p[[i]] <-ggplot(sub, aes(x = Primer_set, y = Abundance, fill = Class))+
    facet_grid(Type~., space= "fixed")+
    geom_bar(stat = "identity", position="fill") +
    scale_y_continuous(breaks=NULL)+
    scale_fill_manual(values = phyla.col)+ 
    coord_polar("y")
}
ml <- marrangeGrob(p, nrow=1, ncol=2)
#do.call(grid.arrange,p)
ml

#calculate abundance for each Order
ps0_glom_Order <- as.data.frame(as.list(aggregate(Abundance~Overlap_ID+Primer_set+Type+Class+Order, ps0_glom_merged.melt,
                                                  FUN = function(x) c(sum = sum(x), count=length(x)))))
ps0_glom_Order$Abundance.sum <- ps0_glom_Order$Abundance.sum*100
ps0_glom_Order<- ps0_glom_Order[ps0_glom_Order$Abundance.sum>0,]

ps0_glom_Order.by_primer <- as.data.frame(as.list(aggregate(Abundance.sum~Primer_set+Type+Class+Order, ps0_glom_Order,
                                                            FUN = function(x) c(mean = mean(x), sd = sd(x), count=length(x)))))


#calculate abundance for each Family
ps0_glom_Family <- as.data.frame(as.list(aggregate(Abundance~Overlap_ID+Primer_set+Type+Class+Order+Family, ps0_glom_merged.melt,
                                                   FUN = function(x) c(sum = sum(x), count=length(x)))))
ps0_glom_Family$Abundance.sum <- ps0_glom_Family$Abundance.sum*100
ps0_glom_Family<- ps0_glom_Family[ps0_glom_Family$Abundance.sum>0,]

ps0_glom_Family.by_primer <- as.data.frame(as.list(aggregate(Abundance.sum~Primer_set+Type+Class+Order+Family, ps0_glom_Family,
                                                             FUN = function(x) c(mean = mean(x), sd = sd(x), count=length(x)))))


#statistical significance of the groups
df <- as(sample_data(ps0_glom_merged.prop), "data.frame")
d <- phyloseq::distance(ps0_glom_merged.prop, "jsd")
adonis_all <- adonis2(d ~ Type+Primer_set , data= df, perm = 999)
adonis_all


#alternative plot
ps0.prop.jsd_slv <- ordinate(ps0_glom_merged.prop, method = "NMDS", 
                                  distance = "jsd")
stressplot(ps0.prop.jsd_slv)

plot_ordination(ps0_glom_merged.prop, ps0.prop.jsd_slv, color = "Type")+
  stat_ellipse(type="norm")+
  geom_point(aes(shape = Primer_set),color ="black", size = 5)+
  geom_point(aes(shape = Primer_set), size = 4)+
    coord_fixed()




#############################
# CARD-FISH
#############################
#import raw data
raw.counts <- read.csv("./data/CARD-FISH/card-fish-counts-SH.csv", dec = ".", stringsAsFactors = FALSE)
raw.counts$conc.FISH <- as.numeric(raw.counts$conc.FISH)


#summarize counts
counts.agg.FISH <- as.data.frame(as.list(aggregate(conc.FISH~StationName+Depth+Domain,raw.counts,FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))))


counts.agg.DAPI <- as.data.frame(as.list(aggregate(conc.DAPI~StationName+Depth+Domain,raw.counts,FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))))



counts.agg <- merge(counts.agg.FISH, counts.agg.DAPI) %>%
  mutate(Proportion = conc.FISH.mean/conc.DAPI.mean)

















#overview plot
counts.agg.abs <- counts.agg[counts.agg$Domain!="EUB",]

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
    scale_y_log10()+
    ylim(0,10e8)+
    theme_classic(base_size = 12)+
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(),
          legend.position="none",
          panel.grid.major=element_blank(),
          strip.background = element_blank())+
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
  
}
ml <- marrangeGrob(p, nrow=4, ncol=1)
ml

