#load libraries
library(phyloseq); packageVersion("phyloseq")

#load scripts
source('./scripts/ReadAmplicon_with_cl_EF.R')
# source("./scripts/phyloseq_to_vegan.R")
# source("./scripts/color_palettes.R")
# source("./scripts/summarize_taxa.R")
# source("./scripts/miseqR.R")
# source("./scripts/lm_eqn.R")

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


#add samples IDs to the metadata
sample_data(V3V4_data.raw)$SampleID <- rownames(sample_data(V3V4_data.raw))
sample_data(V4V5_data.raw)$SampleID <- rownames(sample_data(V4V5_data.raw))

#remove mitochondria and Chloroplast sequences
V4V5_data.BAC <- subset_taxa(V4V5_data.raw, family !="Mitochondria")
V3V4_data.BAC <- subset_taxa(V3V4_data.raw, family !="Mitochondria")

V4V5_data.BAC <- subset_taxa(V4V5_data.raw, order !="Chloroplast")
V3V4_data.BAC <- subset_taxa(V3V4_data.raw, order !="Chloroplast")

#####################################
# export R-native serialized RDS file
#####################################
saveRDS(V3V4_data.BAC, "./Data/V3V4_data_BAC.rds")
saveRDS(V4V5_data.BAC, "./Data/V4V5_data_BAC.rds")
