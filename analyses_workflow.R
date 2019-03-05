#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")
library(reshape2); packageVersion("reshape2")
library(iNEXT); packageVersion("iNEXT")
library(tidyr); packageVersion("tidyr")
library(dplyr); packageVersion("dplyr")
library(plyr); packageVersion("plyr")
library(ggsignif); packageVersion("ggsignif")

#load scripts
source("./scripts/color_palettes.R")
source("./scripts/summarize_taxa.R")
source("./scripts/miseqR.R")
source("./scripts/lm_eqn.R")
source('./scripts/ReadAmplicon_with_cl_EF.R')
source("./scripts/phyloseq_to_vegan.R")

#####################################
#Parse the OTU contigency tables for Phyloseq
#####################################
#V3V4
X <- ReadAmplicon_with_cl(otu = "./Data/Swarm-output/V3V4/OTU_contingency_table.csv", tax = "./Data/Swarm-output/V3V4/amplicons_seeds_taxonomy.txt", 
                          silva = "./silva_tax/silva132_tax_ssu_curated.tsv", 
                          domain = "Bacteria", silva.sep = "\t", singletons = F, unclassified = T, write.files = F)

OTU<- X$OTU
TAX<- X$TAX

# Check order of samples
all.equal(rownames(OTU), rownames(TAX))


# Contextual data
ENV <- read.table("./Data/Swarm-output/V3V4/sample_list_V3V4.txt",header = TRUE, sep = "\t")
rownames(ENV) <- ENV$X
ENV <- ENV[,2:8]

#creating Phyloseq dataset
OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
meta <- sample_data(ENV)
V3V4_data.raw <- phyloseq(OTU, TAX, meta)


#V4V5
X <- ReadAmplicon_with_cl(otu = "./Data/Swarm-output/V4V5/OTU_contingency_table.csv", tax = "./Data/Swarm-output/V4V5/amplicons_seeds_taxonomy.txt", 
                          silva = "./silva_tax/silva132_tax_ssu_curated.tsv", 
                          domain = "Bacteria", silva.sep = "\t", singletons = F, unclassified = T, write.files = F)

OTU<- X$OTU
TAX<- X$TAX

# Check order of samples
all.equal(rownames(OTU), rownames(TAX))


# Contextual data
ENV <- read.table("./Data/Swarm-output/V4V5/sample_list_V4V5.txt",header = TRUE, sep = "\t")
rownames(ENV) <- ENV$X
ENV <- ENV[,2:8]

#creating Phyloseq dataset
OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
meta <- sample_data(ENV)
V4V5_data.raw <- phyloseq(OTU, TAX, meta)

#remove mitochondria and Chloroplast sequences
V4V5_data.BAC <- subset_taxa(V4V5_data.raw, family !="Mitochondria")
V3V4_data.BAC <- subset_taxa(V3V4_data.raw, family !="Mitochondria")

V4V5_data.BAC <- subset_taxa(V4V5_data.raw, order !="Chloroplast")
V3V4_data.BAC <- subset_taxa(V3V4_data.raw, order !="Chloroplast")

#####################################
#Qunatify reads retainment in the bionformatic pipeline
#####################################
#V3V4
#Read counts from the QC
V3V4_swarm_output <- read.table("./Data/Swarm-output/V3V4/nSeqs_V3V4_all.txt", sep="\t", header = FALSE)

#remove paths of files
V3V4_swarm_output$SampleID <- V3V4_swarm_output$V1
V3V4_swarm_output$SampleID <- gsub("./Renamed/","",V3V4_swarm_output$SampleID)
V3V4_swarm_output$SampleID <- gsub("_R1.*","",V3V4_swarm_output$SampleID)
V3V4_swarm_output$SampleID <- paste("X",V3V4_swarm_output$SampleID, sep = "")

names(V3V4_swarm_output) <- c("Raw","Clipped","Trimmed","Merged","SampleID")
for (n in c("Raw","Clipped","Trimmed","Merged")){
  V3V4_swarm_output[,c(n)] <- as.numeric(gsub(".*fastq:","",V3V4_swarm_output[,c(n)]))
}

#add final number of sequences
V3V4_swarm_output <- V3V4_swarm_output[V3V4_swarm_output$SampleID %in% sample_names(V3V4_data.BAC),]
V3V4_swarm_output$Final <- sample_sums(V3V4_data.BAC)
V3V4_swarm_output$Primer <- "V3V4"

#calculate proportions for each sample
V3V4_swarm_output.prop <- sweep(V3V4_swarm_output[, -c(5:6)], 1, V3V4_swarm_output[, "Raw"], "/")
V3V4_swarm_output.prop <- cbind(V3V4_swarm_output[, c(5:6)], V3V4_swarm_output.prop)

#V4V5
#Read counts from the QC
V4V5_swarm_output <- read.table("./Data/Swarm-output/V4V5/nSeqs_V4V5_all.txt", sep="\t", header = FALSE)

#remove paths of files
V4V5_swarm_output$SampleID <- V4V5_swarm_output$V1
V4V5_swarm_output$SampleID <- gsub("./Renamed/","",V4V5_swarm_output$SampleID)
V4V5_swarm_output$SampleID <- gsub("_R1.*","",V4V5_swarm_output$SampleID)
V4V5_swarm_output$SampleID <- paste("X",V4V5_swarm_output$SampleID, sep = "")

names(V4V5_swarm_output) <- c("Raw","Clipped","Trimmed","Merged","SampleID")
for (n in c("Raw","Clipped","Trimmed","Merged")){
  V4V5_swarm_output[,c(n)] <- as.numeric(gsub(".*fastq:","",V4V5_swarm_output[,c(n)]))
}

#add final number of sequences
V4V5_swarm_output$Final <- sample_sums(V4V5_data.BAC)
V4V5_swarm_output <- V4V5_swarm_output[V4V5_swarm_output$SampleID %in% sample_names(V4V5_data.BAC),]
V4V5_swarm_output$Primer <- "V4V5"

#calculate proportions for each sample
V4V5_swarm_output.prop <- sweep(V4V5_swarm_output[, -c(5:6)], 1, V4V5_swarm_output[, "Raw"], "/")
V4V5_swarm_output.prop <- cbind(V4V5_swarm_output[, c(5:6)], V4V5_swarm_output.prop)

#merge tables from both primers
swarm_output <- rbind(V3V4_swarm_output.prop,V4V5_swarm_output.prop)
swarm_output <- melt(subset(swarm_output, select = -c(SampleID)),"Primer")

#summarize
swarm_output.agg <- as.data.frame(as.list(aggregate(value~Primer+variable,swarm_output,FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))))
swarm_output.agg <- cbind(swarm_output.agg[,c(1:2)],sweep(swarm_output.agg[, -c(1:2)], 1, 100, "*"))

#plot
overview.p <- ggplot(swarm_output.agg, aes(x=variable, y = value.mean,group = Primer, fill = Primer))+
  geom_bar(position=position_dodge(),stat="identity")+ 
  geom_errorbar(aes(ymin = value.mean-value.se, ymax = value.mean+value.se,group = Primer), width=.2,position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0,100))+
  ylab("Proportion of retained reads (%)")+
  xlab("Workflow steps")+
  theme_classic(base_size = 10)

#####################################
#Counts of mitochondrial and chloroplast-related sequences
V3V4_mit_chl_counts <- data.frame(Primer = "V3V4",Station = sample_data(V3V4_data.raw)$Station,Type = sample_data(V3V4_data.raw)$Type,Depth = sample_data(V3V4_data.raw)$Depth, 
                                  Mitochondria = 1-sample_sums(subset_taxa(V3V4_data.raw, family !="Mitochondria"))/sample_sums(V3V4_data.raw), 
                                  Chloroplast = 1-sample_sums(subset_taxa(V3V4_data.raw, order !="Chloroplast"))/sample_sums(V3V4_data.raw))

V4V5_mit_chl_counts <- data.frame(Primer = "V4V5", Station = sample_data(V4V5_data.raw)$Station,Type = sample_data(V4V5_data.raw)$Type,Depth = sample_data(V4V5_data.raw)$Depth, 
                                  Mitochondria = 1-sample_sums(subset_taxa(V4V5_data.raw, family !="Mitochondria"))/sample_sums(V4V5_data.raw), 
                                  Chloroplast = 1-sample_sums(subset_taxa(V4V5_data.raw, order !="Chloroplast"))/sample_sums(V4V5_data.raw))

#summarize for both primers
mit.agg <- as.data.frame(as.list(aggregate(Mitochondria~Primer,rbind(V3V4_mit_chl_counts[V3V4_mit_chl_counts$Depth<50,],V4V5_mit_chl_counts[V4V5_mit_chl_counts$Depth<50,]),FUN = function(x) c(mean = mean(x),median = median(x), se = sd(x)/sqrt(length(x))))))
chl.agg <- as.data.frame(as.list(aggregate(Chloroplast~Primer,rbind(V3V4_mit_chl_counts[V3V4_mit_chl_counts$Depth<50,],V4V5_mit_chl_counts[V4V5_mit_chl_counts$Depth<50,]),FUN = function(x) c(mean = mean(x),median = median(x), se = sd(x)/sqrt(length(x))))))


#####################################
#Alpha diversity calculations
#####################################
#V3V4
#Rarefy the dataset by the smallest sample
#rare_size <- min(sample_sums(V3V4_data.BAC))
rare_size <- 12791
# Initialize matrices to store alpha diversity indeces
#Observed number of OTUs
OTUs <- matrix(nrow = nsamp, ncol = trials)
row.names(OTUs) <- sample_names(V3V4_data.BAC)

#Chao1 richness
richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(V3V4_data.BAC)

#Shannon div. index
Shannon <- matrix(nrow = nsamp, ncol = trials)
row.names(Shannon) <- sample_names(V3V4_data.BAC)

#Simpson div.index
Simpson <- matrix(nrow = nsamp, ncol = trials)
row.names(Simpson) <- sample_names(V3V4_data.BAC)


#Resample the diversity indeces 100 times
set.seed(12345)
nsamp <-  nsamples(V3V4_data.BAC)
trials <- 100

for (i in 1:100) {
  V3V4_data.BAC.rare <- rarefy_even_depth(V3V4_data.BAC, sample.size = rare_size,
                                          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  
  
  V3V4_alpha.div <- estimate_richness(V3V4_data.BAC.rare, split = TRUE, measures = NULL)
  
  # OTU num.
  OTUs[ ,i] <- V3V4_alpha.div$Observed
  
  # Calculate richness
  richness[ ,i] <- V3V4_alpha.div$Chao1
  
  # Calculate richness
  Shannon[ ,i] <- V3V4_alpha.div$Shannon
  
  # Calculate richness
  Simpson[ ,i] <- V3V4_alpha.div$Simpson
  
  # Calculate evenness
  evenness[ ,i] <- V3V4_alpha.div$Shannon/log(V3V4_alpha.div$Observed)
}

# Create a new dataframe to hold the means and standard deviations of observed OTUs
SampleID <- row.names(OTUs)
mean <- apply(OTUs, 1, mean)
sd <- apply(OTUs, 1, sd)
measure <- rep("OTUs", nsamp)
OTUs_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of Shannon div. index estimates
SampleID <- row.names(Shannon)
mean <- apply(Shannon, 1, mean)
sd <- apply(Shannon, 1, sd)
measure <- rep("Shannon", nsamp)
Shannon_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of Simpson div. index estimates
SampleID <- row.names(Simpson)
mean <- apply(Simpson, 1, mean)
sd <- apply(Simpson, 1, sd)
measure <- rep("Simpson", nsamp)
Simpson_stats <- data.frame(SampleID, mean, sd, measure)

#summarize in a table
V3V4_comm.char<- data.frame(PANGAEA_sampleID = paste(sample_data(V3V4_data.BAC)$Expedition,sample_data(V3V4_data.BAC)$Station, sep = "/"),
                            #StationName = sample_data(V3V4_data.BAC)$Station,
                            Type = sample_data(V3V4_data.BAC)$Type,
                            Environment = sample_data(V3V4_data.BAC)$merging,
                            Depth = sample_data(V3V4_data.BAC)$Depth,
                            Sequences= sample_sums(V3V4_data.BAC),
                            Observed = round(OTUs_stats$mean,digits=0),
                            Chao1 = rich_stats$mean,
                            Completness = round(100*(OTUs_stats$mean/rich_stats$mean),digits=1),
                            Shannon = Shannon_stats$mean,
                            Simpson = Simpson_stats$mean,
                            Primer="V3-V4")

#####################################
#V4V5

#Rarefy the dataset by the smallest sample
#rare_size <- min(sample_sums(V4V5_data.BAC))
rare_size <- 12791
# Initialize matrices to store alpha diversity indeces
#Observed number of OTUs
OTUs <- matrix(nrow = nsamp, ncol = trials)
row.names(OTUs) <- sample_names(V4V5_data.BAC)

#Chao1 richness
richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(V4V5_data.BAC)

#Shannon div. index
Shannon <- matrix(nrow = nsamp, ncol = trials)
row.names(Shannon) <- sample_names(V4V5_data.BAC)

#Simpson div.index
Simpson <- matrix(nrow = nsamp, ncol = trials)
row.names(Simpson) <- sample_names(V4V5_data.BAC)


#Resample the diversity indeces 100 times
set.seed(12345)
nsamp <-  nsamples(V4V5_data.BAC)
trials <- 100

for (i in 1:100) {
  V4V5_data.BAC.rare <- rarefy_even_depth(V4V5_data.BAC, sample.size = rare_size,
                                          rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
  
  
  V4V5_alpha.div <- estimate_richness(V4V5_data.BAC.rare, split = TRUE, measures = NULL)
  
  # OTU num.
  OTUs[ ,i] <- V4V5_alpha.div$Observed
  
  # Calculate richness
  richness[ ,i] <- V4V5_alpha.div$Chao1
  
  # Calculate richness
  Shannon[ ,i] <- V4V5_alpha.div$Shannon
  
  # Calculate richness
  Simpson[ ,i] <- V4V5_alpha.div$Simpson
  
  # Calculate evenness
  evenness[ ,i] <- V4V5_alpha.div$Shannon/log(V4V5_alpha.div$Observed)
}

# Create a new dataframe to hold the means and standard deviations of observed OTUs
SampleID <- row.names(OTUs)
mean <- apply(OTUs, 1, mean)
sd <- apply(OTUs, 1, sd)
measure <- rep("OTUs", nsamp)
OTUs_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of Shannon div. index estimates
SampleID <- row.names(Shannon)
mean <- apply(Shannon, 1, mean)
sd <- apply(Shannon, 1, sd)
measure <- rep("Shannon", nsamp)
Shannon_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of Simpson div. index estimates
SampleID <- row.names(Simpson)
mean <- apply(Simpson, 1, mean)
sd <- apply(Simpson, 1, sd)
measure <- rep("Simpson", nsamp)
Simpson_stats <- data.frame(SampleID, mean, sd, measure)

#summarize in a table
V4V5_comm.char<- data.frame(PANGAEA_sampleID = paste(sample_data(V4V5_data.BAC)$Expedition,sample_data(V4V5_data.BAC)$Station, sep = "/"),
                            #StationName = sample_data(V4V5_data.BAC)$Station,
                            Type = sample_data(V4V5_data.BAC)$Type,
                            Environment = sample_data(V4V5_data.BAC)$merging,
                            Depth = sample_data(V4V5_data.BAC)$Depth,
                            Sequences= sample_sums(V4V5_data.BAC),
                            Observed = round(OTUs_stats$mean,digits=0),
                            Chao1 = rich_stats$mean,
                            Completness = round(100*(OTUs_stats$mean/rich_stats$mean),digits=1),
                            Shannon = Shannon_stats$mean,
                            Simpson = Simpson_stats$mean,
                            Primer="V4-V5")

#####################################
#export alpha diversity merged table
write.csv(rbind(V3V4_comm.char,V4V5_comm.char), "./Data/alpha_table_12K.csv")

#Alpha diversity significance tests
alpha_together <- rbind(V3V4_comm.char,V4V5_comm.char)

alpha_together$Environment <- factor(alpha_together$Environment,
                                     levels=c("Sea-ice","Surface-water","Deep-water","Sediment"))
alpha_together$Primer <- factor(alpha_together$Primer,
                                     levels=c("V4-V5","V3-V4"))

#Plot alpha diversity 
chao1.p <- ggplot(alpha_together, aes(x = Primer, y = Chao1)) +
  theme_classic(base_size = 12) +
  labs(x = "Primer")+
  geom_boxplot(outlier.color = "black", notch = FALSE)+
  geom_jitter(position=position_jitter(0), alpha =0.2, colour = "gray50")+
  facet_grid(Environment~.)+
  coord_flip()+
  geom_signif(comparisons = list(c("V3-V4", "V4-V5")), 
              map_signif_level=TRUE, test = "wilcox.test")

shannon.p <- ggplot(alpha_together, aes(x = Primer, y = Shannon)) +
  theme_classic(base_size = 12) +
  labs(x = "Primer")+
  geom_boxplot(outlier.color = "black", notch = FALSE)+
  geom_jitter(position=position_jitter(0), alpha =0.2, colour = "gray50")+
  facet_grid(Environment~.)+
  coord_flip()+
  geom_signif(comparisons = list(c("V3-V4", "V4-V5")), 
              map_signif_level=TRUE, test = "wilcox.test")

simpson.p <- ggplot(alpha_together, aes(x = Primer, y = Simpson)) +
  theme_classic(base_size = 12) +
  labs(x = "Primer")+
  geom_boxplot(outlier.color = "black", notch = FALSE)+
  geom_jitter(position=position_jitter(0), alpha =0.2, colour = "gray50")+
  facet_grid(Environment~.)+
  coord_flip()+
  geom_signif(comparisons = list(c("V3-V4", "V4-V5")), 
                           map_signif_level=TRUE, test = "wilcox.test")

plot_grid(chao1.p, shannon.p,simpson.p, ncol =3)


#####################################
#Community composition barplots
#####################################

#Prevalence filter to remove rare taxa
#V3V4
#in how many samples did each taxa appear at least once
prev0 = apply(X = otu_table(V3V4_data.BAC),
              MARGIN = ifelse(taxa_are_rows(V3V4_data.BAC), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

# Define prevalence threshold as 10% of total samples
prevalenceThreshold = 0.1 * nsamples(V3V4_data.BAC)

# Execute prevalence filter, using `prune_taxa()` function
V3V4_data.BAC.pruned <-  prune_taxa((prev0 > prevalenceThreshold), V3V4_data.BAC)

#V4V5
#in how many samples did each taxa appear at least once
prev0 = apply(X = otu_table(V4V5_data.BAC),
              MARGIN = ifelse(taxa_are_rows(V4V5_data.BAC), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

# Define prevalence threshold as 10% of total samples
prevalenceThreshold = 0.1 * nsamples(V4V5_data.BAC)

# Execute prevalence filter, using `prune_taxa()` function
V4V5_data.BAC.pruned <-  prune_taxa((prev0 > prevalenceThreshold), V4V5_data.BAC)

#####################################
#list of the major taxonomic groups
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
#calculate relative abundance and melt to a data frame
V3V4_data.BAC.ra <- transform_sample_counts(V3V4_data.BAC.pruned, function(x) x / sum(x) )
V3V4_data.BAC.melt <- psmelt(V3V4_data.BAC.ra)

#create sample names to merge between the datasets
V3V4_data.BAC.melt$SampleID <- paste(V3V4_data.BAC.melt$Expedition,V3V4_data.BAC.melt$Station,V3V4_data.BAC.melt$Depth,V3V4_data.BAC.melt$Type, sep="+")

#aggregate by taxonomy
names(V3V4_data.BAC.melt)[names(V3V4_data.BAC.melt) == "class"] <- "Taxonomy"
V3V4_data.BAC.melt$Taxonomy <- as.character(V3V4_data.BAC.melt$Taxonomy)
V3V4_data.BAC.melt$Taxonomy[!V3V4_data.BAC.melt$Taxonomy %in% selected_taxa] <- "Other"
V3V4_data.BAC.melt <- arrange(V3V4_data.BAC.melt, Sample, Taxonomy)
V3V4_data.BAC.melt.agg <- aggregate(Abundance~SampleID+Taxonomy+merging, V3V4_data.BAC.melt, FUN= sum)
V3V4_data.BAC.melt.agg$Abundance <- V3V4_data.BAC.melt.agg$Abundance*100
V3V4_data.BAC.melt.agg$Taxonomy <- as.factor(V3V4_data.BAC.melt.agg$Taxonomy)
V3V4_data.BAC.melt.agg$Primer <- "V3V4"


#V4V5
#calculate relative abundance and melt to a data frame
V4V5_data.BAC.ra <- transform_sample_counts(V4V5_data.BAC.pruned, function(x) x / sum(x) )
V4V5_data.BAC.melt <- psmelt(V4V5_data.BAC.ra)

#create sample names to merge between the datasets
V4V5_data.BAC.melt$SampleID <- paste(V4V5_data.BAC.melt$Expedition,V4V5_data.BAC.melt$Station,V4V5_data.BAC.melt$Depth, V4V5_data.BAC.melt$Type, sep="+")

#aggregate by taxonomy
names(V4V5_data.BAC.melt)[names(V4V5_data.BAC.melt) == "class"] <- "Taxonomy"
V4V5_data.BAC.melt$Taxonomy <- as.character(V4V5_data.BAC.melt$Taxonomy)
V4V5_data.BAC.melt$Taxonomy[!V4V5_data.BAC.melt$Taxonomy %in% selected_taxa] <- "Other"
V4V5_data.BAC.melt <- arrange(V4V5_data.BAC.melt, Sample, Taxonomy)
V4V5_data.BAC.melt.agg <- aggregate(Abundance~SampleID+Taxonomy+merging, V4V5_data.BAC.melt, FUN= sum)
V4V5_data.BAC.melt.agg$Abundance <- V4V5_data.BAC.melt.agg$Abundance*100
V4V5_data.BAC.melt.agg$Taxonomy <- as.factor(V4V5_data.BAC.melt.agg$Taxonomy)
V4V5_data.BAC.melt.agg$Primer <- "V4V5"

#merge both datasets togather and redefine factors order (for plotting)
BAC.melt.agg <- rbind(V3V4_data.BAC.melt.agg,V4V5_data.BAC.melt.agg)
#BAC.melt.agg %>% separate(SampleID, c("Expedition", "Station","Depth","Type"),"\\+") ->BAC.melt.agg

BAC.melt.combined.agg <- aggregate(Abundance~SampleID+merging+Primer+Taxonomy, BAC.melt.agg, mean)
BAC.melt.combined.agg$Taxonomy <- gsub("_unclassified| Incertae Sedis", "",BAC.melt.combined.agg$Taxonomy) 
BAC.melt.combined.agg$Taxonomy <- factor(BAC.melt.combined.agg$Taxonomy, levels= c("Alphaproteobacteria","Acidimicrobiia","Bacteroidia","Deltaproteobacteria","Dehalococcoidia","Gammaproteobacteria","Marinimicrobia (SAR406 clade)","Nitrospinia", "Planctomycetacia","Verrucomicrobiae","Other"))
BAC.melt.combined.agg$merging <- factor(BAC.melt.combined.agg$merging, levels= c("Sea-ice", "Surface-water","Deep-water","Sediment"))

#generate plots for each environment separately 
p <- list()
for (i in levels(BAC.melt.combined.agg$merging)){
  sub <- BAC.melt.combined.agg[BAC.melt.combined.agg$merging == i,]
  p[[i]] <-ggplot(sub, aes(x = SampleID, y = Abundance, fill = Taxonomy)) + 
    facet_grid(~Primer) + labs(title = i)+
    geom_bar(stat = "identity",  colour = "black") +
    #scale_y_continuous(breaks=NULL)+
    scale_fill_manual(values = taxa_col) +
    #coord_polar("y", start=0)+
    theme_classic(base_size = 10)+
    theme(axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          #axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          #axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank(),
          strip.background = element_blank())
}

plot_grid(p[[1]],p[[2]],p[[3]], p[[4]], ncol = 2, align = "hv")

#####################################
#####################################
#Taxonomic enrichment test
#####################################
#agglomerate taxa on order level and remove any order with less than 100 reads in the dataset
#V3V4
V3V4_data.BAC.ra.glom <- tax_glom(V3V4_data.BAC, taxrank = "order")
V3V4_data.BAC.ra.glom <-prune_taxa(taxa_sums(V3V4_data.BAC.ra.glom)>100,V3V4_data.BAC.ra.glom)
V3V4_data.BAC.melt <- psmelt(V3V4_data.BAC.ra.glom)

#V4V5
V4V5_data.BAC.ra.glom <- tax_glom(V4V5_data.BAC, taxrank = "order")
V4V5_data.BAC.ra.glom <-prune_taxa(taxa_sums(V4V5_data.BAC.ra.glom)>100,V4V5_data.BAC.ra.glom)
V4V5_data.BAC.melt <- psmelt(V4V5_data.BAC.ra.glom)

#####################################
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
  #generate abundance table for order*Sample
  V3V4_data.BAC.agg <- aggregate(Abundance~Sample+order, V3V4_data.BAC.melt.sub,sum)
  V3V4_data.BAC.agg[,c("Sample","order","Abundance")] %>% spread(Sample, Abundance)-> V3V4_data.BAC.wide
  
  #V4V5
  V4V5_data.BAC.melt.sub$Sample <- paste("V4V5",V4V5_data.BAC.melt.sub$Station,
                                         V4V5_data.BAC.melt.sub$Type,
                                         V4V5_data.BAC.melt.sub$Depth, sep = ".")
  
  #generate abundance table for order*Sample
  V4V5_data.BAC.agg <- aggregate(Abundance~Sample+order, V4V5_data.BAC.melt.sub,sum)
  V4V5_data.BAC.agg[,c("Sample","order","Abundance")] %>% spread(Sample, Abundance)-> V4V5_data.BAC.wide
  
  #merge abundance tables
  deseq.test <- merge(V3V4_data.BAC.wide,V4V5_data.BAC.wide, by= "order")
  rownames(deseq.test) <- deseq.test$order
  deseq.test <- subset(deseq.test, select= -c(order))
  
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
  BAC_sed.DEseq.res.sig.sub <- BAC_sed.DEseq.res.sig[abs(BAC_sed.DEseq.res.sig$log2FoldChange)>1.58,]
  BAC_sed.DEseq.res.sig.sub$Env <- environment[i]
  
  #add taxonomy
  BAC_sed.DEseq.res.sig.sub$order <-rownames(BAC_sed.DEseq.res.sig.sub) 
  BAC_sed.DEseq.res.sig.sub <- as(BAC_sed.DEseq.res.sig.sub,"data.frame")
  BAC_sed.DEseq.res.sig.sub <- merge(BAC_sed.DEseq.res.sig.sub, unique(V3V4_data.BAC.melt[,c("order","class","phylum")]), by = "order")
  rownames(BAC_sed.DEseq.res.sig.sub) <- NULL
  
  #merge all results togather
  merged_enrichment_tests <- rbind(merged_enrichment_tests,BAC_sed.DEseq.res.sig.sub)
}

#some adjustments for plotting
merged_enrichment_tests$class <- as.character(merged_enrichment_tests$class)
merged_enrichment_tests_other <- merged_enrichment_tests[!merged_enrichment_tests$class %in% selected_taxa,]

merged_enrichment_tests$class[!merged_enrichment_tests$class %in% selected_taxa] <- "Other"

merged_enrichment_tests$class <- factor(merged_enrichment_tests$class, ordered = TRUE,
                                        levels = unique(merged_enrichment_tests$class[order(merged_enrichment_tests$class)]))
merged_enrichment_tests$order <- factor(merged_enrichment_tests$order, ordered = TRUE,
                                        levels = unique(merged_enrichment_tests$order[order(merged_enrichment_tests$class)]))
merged_enrichment_tests$Env <- factor(merged_enrichment_tests$Env, ordered = TRUE,
                                      levels = c("Sea-ice","Surface-water","Deep-water","Sediment"))
merged_enrichment_tests_other$Env <- factor(merged_enrichment_tests_other$Env, ordered = TRUE,
                                            levels = c("Sea-ice","Surface-water","Deep-water","Sediment"))

#plot the major taxonomic groups
BAC_da_taxa.p <- ggplot(data= merged_enrichment_tests[merged_enrichment_tests$class!="Other",], 
                      aes(x = log2FoldChange, y = order))+ 
                geom_point(data= merged_enrichment_tests[merged_enrichment_tests$class!="Other",], 
                                      aes(x = log2FoldChange, y = order), colour = "black",size = 5)+
                geom_point(data= merged_enrichment_tests[merged_enrichment_tests$class!="Other",], 
                                       aes(x = log2FoldChange, y = order, colour = class), size = 4)+
                ylab("log2foldchange")+
                geom_vline(aes(xintercept=0), linetype="dashed")+
                scale_colour_manual(values = taxa_col)+
                facet_grid(~Env)+
                theme_classic(base_size = 10)+
                theme(axis.text.x =element_text(angle=90) , legend.position = "bottom",
                      panel.grid.major = element_line(colour = "red", linetype = "dotted"))


#plot the other taxonomic groups
colourCount = length(unique(merged_enrichment_tests_other$phylum))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
getPalette = colorRampPalette(tol21rainbow)

BAC_daOTU.other <- ggplot(data= merged_enrichment_tests_other, 
                          aes(x = log2FoldChange, y = order))+ 
                   geom_point(data= merged_enrichment_tests_other, 
                              aes(x = log2FoldChange, y = order), colour = "black",size = 5)+
                   geom_point(data= merged_enrichment_tests_other, 
                              aes(x = log2FoldChange, y = order, colour = phylum), size = 4)+
                    ylab("log2foldchange")+
                    geom_vline(aes(xintercept=0), linetype="dashed")+
                    scale_colour_manual(values = getPalette(colourCount))+
                    facet_grid(~Env)+
                    theme_classic(base_size = 10)+
                    theme(axis.text.x =element_text(angle=90) , legend.position = "bottom",
                          panel.grid.major = element_line(colour = "red", linetype = "dotted"))

#####################################
#some summary numbers
merged_enrichment_tests.sub <- unique(subset(merged_enrichment_tests,select=-c(Env)), by = "order")


overview_order <- as.data.frame(as.list(aggregate(log2FoldChange~class+order,merged_enrichment_tests.sub,
                                            FUN = function(x) c(sum = sum(x),  count=length(x)))))

overview_class <- as.data.frame(as.list(aggregate(log2FoldChange.sum~class,overview,
                                            FUN = function(x) c(sum = sum(x),  count=length(x)))))

#####################################
#Compare to CARD-FISH
#####################################
#import raw counts
raw.counts <- read.csv("../CARD-FISH/card-fish-counts.csv", dec = ".", stringsAsFactors = FALSE)
raw.counts$conc.FISH <- as.numeric(raw.counts$conc.FISH)

#calculate mean and se for counts of ech taxonomic group
counts.agg <- as.data.frame(as.list(aggregate(conc.FISH~StationName+Depth+Domain,raw.counts,FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))))

#exclude the total bacterial counts for the plot
counts.agg.abs <- counts.agg[counts.agg$Domain!="EUB",]

levels(counts.agg.abs$Depth) <- c("ice","DCM","MESO","SED")
counts.agg.abs$Domain <-  factor(counts.agg.abs$Domain,
                                 levels= c("SAR11","BACT","POL","GAM","ALT"))

#plot absolute abundance
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

#####################################
#calculate count proportions 
subset(counts.agg, select = -conc.FISH.se)%>% spread(Domain, conc.FISH.mean)-> counts.agg.wide

counts.agg.wide.ra <- counts.agg.wide
counts.agg.wide.ra$SAR11 <- counts.agg.wide.ra$SAR11/counts.agg.wide.ra$EUB
counts.agg.wide.ra$GAM <- counts.agg.wide.ra$GAM/counts.agg.wide.ra$EUB
counts.agg.wide.ra$ALT <- counts.agg.wide.ra$ALT/counts.agg.wide.ra$EUB
counts.agg.wide.ra$BACT <- counts.agg.wide.ra$BACT/counts.agg.wide.ra$EUB
counts.agg.wide.ra$POL <- counts.agg.wide.ra$POL/counts.agg.wide.ra$EUB

#melt the table for plotting
counts.agg.ra <- melt(subset(counts.agg.wide.ra, select = -EUB),
                      id.vars = c("StationName","Depth"))

#adjust some factors for plottting
levels(counts.agg.ra$Depth) <- c("ice","DCM","MESO","SED")
counts.agg.ra$variable <-  factor(counts.agg.ra$variable,
                                  levels= c("SAR11","BACT","POL","GAM","ALT"))

#plot relative abundance 
p <- list()
for (i in levels(counts.agg.ra$Depth)){
  p[[i]] <- ggplot()+
    geom_bar(aes(y = value, x = variable, fill =variable), colour = "black", data = counts.agg.ra[counts.agg.ra$Depth==i,], stat="identity")+
    facet_grid(.~StationName)+
    #scale_y_log10()+
    ggtitle(label = i)+
    #geom_errorbar(aes(x = Domain, ymin = conc.FISH.mean-conc.FISH.se, ymax = conc.FISH.mean+conc.FISH.se), data = counts.agg.abs[counts.agg.abs$Depth==i,], width=.2,position=position_dodge(.9)) +
    scale_fill_manual(values = cbbPalette)+
    theme_classic(base_size = 12)+
    theme(axis.title.x = element_blank(), axis.text.x=element_blank(),
          legend.position="bottom",
          panel.grid.major=element_blank(),
          strip.background = element_blank())
  
}
ml.ra <- marrangeGrob(p, nrow=4, ncol=1)

#####################################
#get the relevant 16S samples
#V3V4
#extract the relevant samples
V3V4_data.BAC.fish <- subset_samples(V3V4_data.BAC.pruned, FISH =="1")
V3V4_data.BAC.fish <- prune_taxa(taxa_sums(V3V4_data.BAC.fish)>0,V3V4_data.BAC.fish)

#transform to relative sequence abudnance and melt to a data frame
V3V4_data.melt.fish <-transform_and_melt(V3V4_data.BAC.fish, taxrank= "genus", prune= 0.00)

#aggregate on different taxonomic levels
V3V4_class.agg <- as.data.frame(as.list(aggregate(Abundance~Station+merging+class,V3V4_data.melt.fish, FUN = function(x) c(sum = sum(x),  count=length(x)))))
names(V3V4_class.agg) <- c("Station","Env","Taxa","Abundance","No.OTU")
V3V4_order.agg <- as.data.frame(as.list(aggregate(Abundance~Station+merging+order,V3V4_data.melt.fish, FUN = function(x) c(sum = sum(x),  count=length(x)))))
names(V3V4_order.agg) <- c("Station","Env","Taxa","Abundance","No.OTU")
V3V4_family.agg <- as.data.frame(as.list(aggregate(Abundance~Station+merging+family,V3V4_data.melt.fish, FUN = function(x) c(sum = sum(x), count=length(x)))))
names(V3V4_family.agg) <- c("Station","Env","Taxa","Abundance","No.OTU")
V3V4_data.melt.fish$Taxonomy <- gsub("Polaribacter.*|\\[Polaribacter.*","Polaribacter",V3V4_data.melt.fish$Taxonomy)
V3V4_genus.agg <- as.data.frame(as.list(aggregate(Abundance~Station+merging+Taxonomy,V3V4_data.melt.fish, FUN = function(x) c(sum = sum(x),  count=length(x)))))
names(V3V4_genus.agg) <- c("Station","Env","Taxa","Abundance","No.OTU")

#merge all taxa to one data frame
V3V4_summary.agg <-rbind(V3V4_class.agg,V3V4_order.agg,V3V4_family.agg,V3V4_genus.agg)
V3V4_summary.agg$Method <- "V3V4"

#V4V5
#extract the relevant samples
V4V5_data.BAC.fish <- subset_samples(V4V5_data.BAC.pruned, FISH =="1")
V4V5_data.BAC.fish <- prune_taxa(taxa_sums(V4V5_data.BAC.fish)>0,V4V5_data.BAC.fish)

#transform to relative sequence abudnance and melt to a data frame
V4V5_data.melt.fish <-transform_and_melt(V4V5_data.BAC.fish, taxrank= "genus", prune= 0.00)

#aggregate on different taxonomic levels
V4V5_class.agg <- as.data.frame(as.list(aggregate(Abundance~Station+merging+class,V4V5_data.melt.fish, FUN = function(x) c(sum = sum(x),  count=length(x)))))
names(V4V5_class.agg) <- c("Station","Env","Taxa","Abundance","No.OTU")
V4V5_order.agg <- as.data.frame(as.list(aggregate(Abundance~Station+merging+order,V4V5_data.melt.fish, FUN = function(x) c(sum = sum(x),  count=length(x)))))
names(V4V5_order.agg) <- c("Station","Env","Taxa","Abundance","No.OTU")
V4V5_family.agg <- as.data.frame(as.list(aggregate(Abundance~Station+merging+family,V4V5_data.melt.fish, FUN = function(x) c(sum = sum(x), count=length(x)))))
names(V4V5_family.agg) <- c("Station","Env","Taxa","Abundance","No.OTU")
V4V5_data.melt.fish$Taxonomy <- gsub("Polaribacter.*|\\[Polaribacter.*","Polaribacter",V4V5_data.melt.fish$Taxonomy)
V4V5_genus.agg <- as.data.frame(as.list(aggregate(Abundance~Station+merging+Taxonomy,V4V5_data.melt.fish, FUN = function(x) c(sum = sum(x),  count=length(x)))))
names(V4V5_genus.agg) <- c("Station","Env","Taxa","Abundance","No.OTU")

#merge all taxa to one data frame
V4V5_summary.agg <-rbind(V4V5_class.agg,V4V5_order.agg,V4V5_family.agg,V4V5_genus.agg)
V4V5_summary.agg$Method <- "V4V5"


#extract only the taxonomic groups with counts
FISH_taxa <- c("SAR11 clade", "Alteromonadales","Polaribacter",
               "Bacteroidia",
               "Gammaproteobacteria")

#FISH_taxa <- c("SAR11 clade", "Alteromonadales","Polaribacter")
V3V4_summary.agg.sub <- V3V4_summary.agg[V3V4_summary.agg$Taxa %in% FISH_taxa,]
V4V5_summary.agg.sub <- V4V5_summary.agg[V4V5_summary.agg$Taxa %in% FISH_taxa,]


#merge both datasets togather
summary.agg.sub <- rbind(V3V4_summary.agg.sub,V4V5_summary.agg.sub)
summary.agg.sub$Abundance <- as.numeric(summary.agg.sub$Abundance)
summary.agg.sub$Taxa <- factor(summary.agg.sub$Taxa,
                               levels = c("SAR11 clade", "Bacteroidia",
                                          "Gammaproteobacteria", 
                                          "Polaribacter",
                                          "Alteromonadales"))

#plot the 16S relative abundance
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
          legend.position="none",
          panel.grid.major=element_blank(),
          strip.background = element_blank())
}

ml.16s <- marrangeGrob(q, nrow=4, ncol=1)
ml.16s
#####################################
test.fish <- counts.agg.abs
test.fish$conc.FISH.se <- "FISH"
names(test.fish) <- c("Station","Env","Taxa","Abundance","Method")

revalue(test.fish$Taxa, c("ALT" = "Alteromonadales", 
                          "GAM"= "Gammaproteobacteria",
                          "SAR11"= "SAR11 clade",
                          "BACT" = "Bacteroidia",
                          "POL"="Polaribacter")) -> test.fish$Taxa


revalue(test.fish$Env, c("ice" = "Sea-ice", 
                         "DCM"= "Surface-water",
                         "MESO"= "Deep-water",
                         "SED" = "Sediment")) -> test.fish$Env


test <- rbind(summary.agg.sub[,c("Station","Env","Taxa","Abundance","Method")], test.fish)


test %>% spread(Method, Abundance) ->test.wide

ddply(test.wide, .(Taxa), summarise, 
      "corr" = cor.test(FISH, V3V4, method = "spearman")$estimate,
      "p.value"= cor.test(FISH, V3V4, method = "spearman")$p.value)

ddply(test.wide, .(Taxa), summarise, 
      "corr" = cor.test(FISH, V4V5, method = "spearman")$estimate,
      "p.value"= cor.test(FISH, V4V5, method = "spearman")$p.value)

#####################################
#Session info
####################################
sessionInfo()
