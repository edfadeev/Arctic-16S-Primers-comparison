#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

#load functions
source('./scripts/ReadAmplicon_with_cl.R')
source("./scripts/lm_eqn.R")
source("./scripts/phyloseq_to_vegan.R")

#####################################
#Parse V3V4 OTU for Phyloseq
#####################################
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

#Remove the fractionated samples
V3V4_data.raw <- subset_samples(V3V4_data.raw, Type != "Water-FL")
V3V4_data.raw <- subset_samples(V3V4_data.raw, Type != "Water-PA")
V3V4_data.raw <- subset_samples(V3V4_data.raw, merging != "Exclude")
V3V4_data.raw <- prune_taxa(taxa_sums(V3V4_data.raw)>0, V3V4_data.raw)
V3V4_data.BAC <- V3V4_data.raw

#####################################
#Parse V4V5 OTU for Phyloseq
#####################################
Y <- ReadAmplicon_with_cl(otu = "./Data/Swarm-output/V4V5/OTU_contingency_table.csv", tax = "./Data/Swarm-output/V4V5/amplicons_seeds_taxonomy.txt", 
                  silva = "./silva_tax/silva132_tax_ssu_curated.tsv", 
                  domain = "Bacteria", silva.sep = "\t", singletons = F, unclassified = T, write.files = F)

OTU<- Y$OTU
TAX<- Y$TAX

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

#Remove the fractionated samples
V4V5_data.raw <- subset_samples(V4V5_data.raw, Type != "Water-FL")
V4V5_data.raw <- subset_samples(V4V5_data.raw, Type != "Water-PA")
V4V5_data.raw <- subset_samples(V4V5_data.raw, merging != "Exclude")
V4V5_data.raw <- prune_taxa(taxa_sums(V4V5_data.raw)>0, V4V5_data.raw)
V4V5_data.BAC <- V4V5_data.raw

#remove mitochondria and Chloroplast sequences
V4V5_data.BAC <- subset_taxa(V4V5_data.BAC, family !="Mitochondria")
V3V4_data.BAC <- subset_taxa(V3V4_data.BAC, family !="Mitochondria")

V4V5_data.BAC <- subset_taxa(V4V5_data.BAC, order !="Chloroplast")
V3V4_data.BAC <- subset_taxa(V3V4_data.BAC, order !="Chloroplast")

#save the phyloseq objects
saveRDS(V3V4_data.BAC, "./Data/V3V4_data_BAC.rds")
saveRDS(V4V5_data.BAC, "./Data/V4V5_data_BAC.rds")

saveRDS(V3V4_data.raw, "./Data/V3V4_data_raw.rds")
saveRDS(V4V5_data.raw, "./Data/V4V5_data_raw.rds")

#####################################
#remove all temporary datasets 
####################################
rm(list=ls(all=TRUE))

sessionInfo()

