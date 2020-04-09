#set working directory
wd <- dirname("D:/Postdoc-Vienna/MPI-projects/Primers_comparison/new_analysis_030420/")
setwd(wd)

#load libraries
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(iNEXT); packageVersion("iNEXT")

theme_set(theme_bw())

#####################################
#Import dada2 output into phyloseq
#####################################
#load V3V4 dataset
load("./Primers_comparison/new_analysis_030420/data/V3V4/V3V4_dada2.Rdata")

sample_df <- read.csv("./Primers_comparison/new_analysis_030420/data/V3V4/V3V4_samples_list.csv", header = TRUE, row.names = "SampleID")


V3V4_ps <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE), 
               sample_data(sample_df), 
               tax_table(taxa)#, phy_tree(fitGTR$tree)
               )

#add reference sequence and replace variants with ASVs
dna <- Biostrings::DNAStringSet(taxa_names(V3V4_ps))
names(dna) <- taxa_names(V3V4_ps)
V3V4_ps <- merge_phyloseq(V3V4_ps, dna)
taxa_names(V3V4_ps) <- paste0("ASV", seq(ntaxa(V3V4_ps)))

#remove the unsuccesful sample 9
V3V4_ps<- subset_samples(V3V4_ps, !sample_title %in% c("FRAM_2963_Fevi_32_oben_2"))

#check how many ASVs were unclassified on phylum level, or assigned to chloroplast and Mitochondria
V3V4_ps_uncl <- subset_taxa(V3V4_ps, is.na(Phylum))

V3V4_ps_chl <- subset_taxa(V3V4_ps, Order %in% c("Chloroplast"))

V3V4_ps_mit <- subset_taxa(V3V4_ps, Family %in% c("Mitochondria") )

#remove unclassified on phylum level, chloroplast and Mitochondrial sequence variants
V3V4_ps0 <- subset_taxa(V3V4_ps, !is.na(Phylum) & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") )

#calculate proportion of excludes sequences
1-sum(sample_sums(V3V4_ps0))/sum(sample_sums(V3V4_ps))

dim(otu_table(V3V4_ps0))[2]/dim(otu_table(V3V4_ps))[2]

#calculate retained reads along the workflow
dadaFs <- append(dadaFs_hi, dadaFs_mi)
dadaRs <- append(dadaRs_hi, dadaRs_mi)
mergers <-  append(mergers_hi, mergers_mi)
out <- rbind(out_hi,out_mi)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "tabled")
rownames(track) <- sample.names

write.csv(track, file.path("./Primers_comparison/new_analysis_030420/data/V3V4/", "Report","dada2_reads_output.csv"))


rm(list=ls()[! ls() %in% c("V3V4_ps0")])


#load V4V5 dataset
load("./Primers_comparison/new_analysis_030420/data/V4V5/V4V5_dada2.Rdata")

sample_df <- read.csv("./Primers_comparison/new_analysis_030420/data/V4V5/V4V5_samples_list.csv", header = TRUE, row.names = "SampleID")


V4V5_ps <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE), 
                    sample_data(sample_df), 
                    tax_table(taxa)#, phy_tree(fitGTR$tree)
                    )

#add reference sequence and replace variants with ASVs
dna <- Biostrings::DNAStringSet(taxa_names(V4V5_ps))
names(dna) <- taxa_names(V4V5_ps)
V4V5_ps <- merge_phyloseq(V4V5_ps, dna)
taxa_names(V4V5_ps) <- paste0("ASV", seq(ntaxa(V4V5_ps)))

#remove the unsuccesful sample 9
V4V5_ps<- subset_samples(V4V5_ps, !sample_title %in% c("FRAM_2963_Fevi_32_oben_2"))

#check how many ASVs were unclassified on phylum level, or assigned to chloroplast and Mitochondria
V4V5_ps_uncl <- subset_taxa(V4V5_ps, is.na(Phylum))

V4V5_ps_chl <- subset_taxa(V4V5_ps, Order %in% c("Chloroplast"))

V4V5_ps_mit <- subset_taxa(V4V5_ps, Family %in% c("Mitochondria") )

#remove them
V4V5_ps0 <- subset_taxa(V4V5_ps, !is.na(Phylum) & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") )

#calculate proportion of excludes sequences
1-sum(sample_sums(V4V5_ps0))/sum(sample_sums(V4V5_ps))

dim(otu_table(V4V5_ps0))[2]/dim(otu_table(V4V5_ps))[2]

#calculate retained reads along the workflow
dadaFs <- append(dadaFs_hi, dadaFs_mi)
dadaRs <- append(dadaRs_hi, dadaRs_mi)
mergers <-  append(mergers_hi, mergers_mi)
out <- rbind(out_hi,out_mi)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "tabled")
rownames(track) <- sample.names

write.csv(track, file.path("./Primers_comparison/new_analysis_030420/data/V4V5/", "Report","dada2_reads_output.csv"))


rm(list=ls()[! ls() %in% c("V3V4_ps0","V4V5_ps0")])


#####################################
#Alpha diversity overview table
####################################
#V3V4
V3V4_reads.tab <- read.csv2("./Primers_comparison/new_analysis_030420/data/V3V4/Report/dada2_reads_output.csv",
                            header = TRUE, sep = ",", row.names = "X")

row.names(V3V4_reads.tab) <-paste("X",row.names(V3V4_reads.tab), sep ="")

V3V4_meta <- sample_data(V3V4_ps0)
row.names(V3V4_meta) <-paste("X",row.names(V3V4_meta), sep ="")


V3V4_alpha <- estimate_richness(V3V4_ps0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
V3V4_alpha$Row.names <- row.names(V3V4_alpha)


V3V4_comm.char <- merge(V3V4_meta,V3V4_reads.tab,by =0) %>%
    mutate_if(is.numeric, round, 0) %>%
      merge(V3V4_alpha, by ="Row.names") %>%
      mutate_if(is.numeric, round, 1) %>%
      mutate(Seq.prop = tabled/input) %>%
      mutate(Seq.prop = round(Seq.prop,2)) %>%
      select("sample_title","Primer_set","Type","Depth","input","filtered","merged","tabled",
             "Seq.prop","Observed","Chao1","Shannon","InvSimpson")%>%
      rename("Sample name" = "sample_title",
             "Primer set" = "Primer_set",
             "Microbiome" = "Type",
             "Sampling depth [m]" = "Depth",
             "Raw seq." = "input",
             "Seq. after qual. filt" = "filtered",
             "Merged seq." = "merged",
             "Final seq." = "tabled",
             "Seq. proportions" = "Seq.prop",
             "Observed ASVs" = "Observed",
             "Chao1 richness est." = "Chao1",
             "Shannon Index" = "Shannon",
             "Inverse Simpson Index" = "InvSimpson")


#V4V5
V4V5_reads.tab <- read.csv2("./Primers_comparison/new_analysis_030420/data/V4V5/Report/dada2_reads_output.csv",
                            header = TRUE, sep = ",", row.names = "X")

row.names(V4V5_reads.tab) <-paste("X",row.names(V4V5_reads.tab), sep ="")

V4V5_meta <- sample_data(V4V5_ps0)
row.names(V4V5_meta) <-paste("X",row.names(V4V5_meta), sep ="")


V4V5_alpha <- estimate_richness(V4V5_ps0, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
V4V5_alpha$Row.names <- row.names(V4V5_alpha)


V4V5_comm.char <- merge(V4V5_meta,V4V5_reads.tab,by =0) %>%
  mutate_if(is.numeric, round, 0) %>%
  merge(V4V5_alpha, by ="Row.names") %>%
  mutate_if(is.numeric, round, 1) %>%
  mutate(Seq.prop = tabled/input) %>%
  mutate(Seq.prop = round(Seq.prop,2)) %>%
  select("sample_title","Primer_set","Type","Depth","input","filtered","merged","tabled",
         "Seq.prop","Observed","Chao1","Shannon","InvSimpson")%>%
  rename("Sample name" = "sample_title",
         "Primer set" = "Primer_set",
         "Microbiome" = "Type",
         "Sampling depth [m]" = "Depth",
         "Raw seq." = "input",
         "Seq. after qual. filt" = "filtered",
         "Merged seq." = "merged",
         "Final seq." = "tabled",
         "Seq. proportions" = "Seq.prop",
         "Observed ASVs" = "Observed",
         "Chao1 richness est." = "Chao1",
         "Shannon Index" = "Shannon",
         "Inverse Simpson Index" = "InvSimpson")


overview_table <- rbind(V3V4_comm.char, V4V5_comm.char)

as.data.frame(as.list(aggregate(`Seq. proportions`~`Primer set`,
                                overview_table, 
                                FUN = function(x) c(mean = mean(x), sd = sd(x), count=length(x)))))


#####################################
#Plot rarefaction
####################################
#V3V4
V3V4_ps0.iNEXT <- iNEXT(as.data.frame(t(otu_table(V3V4_ps0, taxa_are_rows=FALSE))), q=0, datatype="abundance")

V3V4_ps0.iNEXT.rare <-fortify(V3V4_ps0.iNEXT, type=1)

V3V4_ps0.iNEXT.meta <- as(sample_data(V3V4_ps0), "data.frame")
V3V4_ps0.iNEXT.meta$site <- rownames(V3V4_ps0.iNEXT.meta)

V3V4_ps0.iNEXT.rare$Type <- V3V4_ps0.iNEXT.meta$Type[match(V3V4_ps0.iNEXT.rare$site, V3V4_ps0.iNEXT.meta$site)] 
V3V4_ps0.iNEXT.rare$SampleName <- V3V4_ps0.iNEXT.meta$sample_title[match(V3V4_ps0.iNEXT.rare$site, V3V4_ps0.iNEXT.meta$site)] 

V3V4_ps0.iNEXT.rare.point <- V3V4_ps0.iNEXT.rare[which(V3V4_ps0.iNEXT.rare$method == "observed"),]
V3V4_ps0.iNEXT.rare.line <- V3V4_ps0.iNEXT.rare[which(V3V4_ps0.iNEXT.rare$method != "observed"),]
V3V4_ps0.iNEXT.rare.line$method <- factor (V3V4_ps0.iNEXT.rare.line$method,
                                           c("interpolated", "extrapolated"),
                                           c("interpolation", "extrapolation"))

V3V4_ps0.iNEXT.rare.line$Primer <- "V3V5"
V3V4_ps0.iNEXT.rare.point$Primer <- "V3V5"

#V4V5
V4V5_ps0.iNEXT <- iNEXT(as.data.frame(t(otu_table(V4V5_ps0, taxa_are_rows=FALSE))), q=0, datatype="abundance")

V4V5_ps0.iNEXT.rare <-fortify(V4V5_ps0.iNEXT, type=1)

V4V5_ps0.iNEXT.meta <- as(sample_data(V4V5_ps0), "data.frame")
V4V5_ps0.iNEXT.meta$site <- rownames(V4V5_ps0.iNEXT.meta)

V4V5_ps0.iNEXT.rare$Type <- V4V5_ps0.iNEXT.meta$Type[match(V4V5_ps0.iNEXT.rare$site, V4V5_ps0.iNEXT.meta$site)] 
V4V5_ps0.iNEXT.rare$SampleName <- V4V5_ps0.iNEXT.meta$sample_title[match(V4V5_ps0.iNEXT.rare$site, V4V5_ps0.iNEXT.meta$site)] 

V4V5_ps0.iNEXT.rare.point <- V4V5_ps0.iNEXT.rare[which(V4V5_ps0.iNEXT.rare$method == "observed"),]
V4V5_ps0.iNEXT.rare.line <- V4V5_ps0.iNEXT.rare[which(V4V5_ps0.iNEXT.rare$method != "observed"),]
V4V5_ps0.iNEXT.rare.line$method <- factor (V4V5_ps0.iNEXT.rare.line$method,
                            c("interpolated", "extrapolated"),
                            c("interpolation", "extrapolation"))

V4V5_ps0.iNEXT.rare.line$Primer <- "V4V5"
V4V5_ps0.iNEXT.rare.point$Primer <- "V4V5"

iNEXT.rare.line <- rbind(V3V4_ps0.iNEXT.rare.line,V4V5_ps0.iNEXT.rare.line)
iNEXT.rare.point <- rbind(V3V4_ps0.iNEXT.rare.point,V4V5_ps0.iNEXT.rare.point)

levels(iNEXT.rare.line$Type) <- c("Sea ice", "Water", "Sediment.trap","Sediment")

rare.p <- ggplot(V4V5_ps0.iNEXT.rare, aes(x=x, y=y, shape = site))+
  geom_line(aes(linetype = method, colour = Primer), lwd = 0.5, data= iNEXT.rare.line)+
  geom_point(shape = 21, size =3, colour = "black", data= iNEXT.rare.point)+
  labs(x = "Sample size", y = "Species richness")+
  xlim(0,2e5)+
  facet_wrap(Type~Primer, ncol = 2,scales = "free")+
  theme_classic(base_size = 12)+theme(legend.position="none")


#####################################
#Filter ASVs of each taxa by prevalence 
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

#explore results
prevdf_V3V4_sum <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
  arrange(desc(Abundance))%>%
  mutate(Proportion = 100*Abundance / sum(Abundance))

top_phyla_V3V4 <- prevdf_V3V4_sum$Phylum[prevdf_V3V4_sum$Proportion >1]

#explore visualy 
ggplot(prevdf_V3V4[prevdf_V3V4$Phylum %in% top_phyla_V3V4,], aes(TotalAbundance, Prevalence ,color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 3, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 100, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Define prevalence threshold as 10% of total samples
prevalenceThreshold = 0.1 * nsamples(V3V4_ps0)

# Execute prevalence filter, using `prune_taxa()` function
V3V4_ps0.prev <-  prune_taxa((prev0_V3V4 > prevalenceThreshold), V3V4_ps0)

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

#explore results
prevdf_V4V5_sum <- plyr::ddply(prevdf_V4V5, "Phylum", function(df1){cbind(Samples=mean(df1$Prevalence),Abundance=sum(df1$TotalAbundance))})%>%
                arrange(desc(Abundance))%>%
                mutate(Proportion = 100*Abundance / sum(Abundance))

top_phyla_V4V5 <- prevdf_V4V5_sum$Phylum[prevdf_V4V5_sum$Proportion >1]
  

prevdf <- rbind(prevdf_V3V4,prevdf_V4V5)


#explore visualy 
ggplot(prevdf[prevdf$Phylum %in% top_phyla_V4V5,], aes(TotalAbundance, Prevalence ,color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 3, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 100, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum+Primer, ncol = 4) + theme(legend.position="none")

# Define prevalence threshold as 10% of total samples
prevalenceThreshold = 0.1 * nsamples(V4V5_ps0)

# Execute prevalence filter, using `prune_taxa()` function
V4V5_ps0.prev <-  prune_taxa((prev0_V4V5 > prevalenceThreshold), V4V5_ps0)


vegan::adonis(TotalAbundance ~ Phylum , prevdf) 

test <- prevdf[prevdf$Phylum %in% top_phyla_V3V4,]%>% 
group_by(Primer, Phylum) %>% 
  summarise(TotalAbundance = list(TotalAbundance))%>% 
spread(Primer, TotalAbundance) %>% 
  group_by(Phylum)%>% 
  mutate(p_value = t.test(unlist("V3V4"), unlist("V4V5"))$p.value,
         t_value = t.test(unlist("V3V4"), unlist("V4V5"))$statistic)


prevdf[prevdf$Phylum %in% top_phyla_V3V4,]%>% 
  select("TotalAbundance","Phylum","Primer")%>% 
  gather(key = Phylum, value = TotalAbundance, -Primer) %>%
  group_by(Primer, TotalAbundance) %>% 
  summarise(TotalAbundance = list(TotalAbundance)) %>% 
  spread(Primer, TotalAbundance) %>%
  group_by(Phylum) %>% 
  mutate(p_value = t.test(unlist("V3V4"), unlist("V4V5"))$p.value,
         t_value = t.test(unlist("V3V4"), unlist("V4V5"))$statistic)



summarise_each(funs(t.test(.[Primer == "V3V4"], .[Primer == "V4V5"])$p.value), vars = Phylum:TotalAbundance)



#####################################
#Bar plots
#####################################


top50 <- names(sort(taxa_sums(V3V4_ps0.prev), decreasing=TRUE))[1:100]
V3V4_ps0.prop <- transform_sample_counts(V3V4_ps0.prev, function(otu) otu/sum(otu))
ps.top50 <- prune_taxa(top50, V3V4_ps0.prop)
plot_bar(ps.top50,  fill="Class") + facet_wrap(~environmental.package, scales="free_x")


ord.nmds.bray <- ordinate(V3V4_ps0.prop, method="NMDS", distance="bray")

plot_ordination(V3V4_ps0.prop, ord.nmds.bray, label = "geographic.location..depth.", shape="environmental.package", title="Bray NMDS")+ geom_point(size = 5)



#####################################
#Import dada2 output into phyloseq
#####################################
#load dataset
load("./new_analysis_030420/data/V4V5/V4V5_dada2.Rdata")

sample_df <- read.csv("new_analysis_030420/data/V4V5/V4V5_samples_list.csv", header = TRUE, row.names = "SampleID")


V4V5_ps <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE), 
                    sample_data(sample_df), 
                    tax_table(taxa), phy_tree(fitGTR$tree))

#add reference sequence and replace variants with ASVs
dna <- Biostrings::DNAStringSet(taxa_names(V4V5_ps))
names(dna) <- taxa_names(V4V5_ps)
V4V5_ps <- merge_phyloseq(V4V5_ps, dna)
taxa_names(V4V5_ps) <- paste0("ASV", seq(ntaxa(V4V5_ps)))

#remove unclassified on phylum level, chloroplast and Mitochondrial sequence variants
V4V5_ps0 <- subset_taxa(V4V5_ps, !is.na(Phylum) & !Order %in% c("Chloroplast") & !Family %in% c("Mitochondria") )

rm(list=ls()[! ls() %in% c("V4V5_ps0")])

#####################################
#Preprocess bacterial OTU by prevalence of each taxa
#####################################
#in how many samples did each taxa appear at least once
prev0 <- apply(X = otu_table(V4V5_ps0),
               MARGIN = ifelse(taxa_are_rows(V4V5_ps0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prev0,
                     TotalAbundance = taxa_sums(V4V5_ps0),
                     tax_table(V4V5_ps0))

#explore results
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#explore visualy 
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(V4V5_ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept = 100, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Define prevalence threshold as 10% of total samples
prevalenceThreshold = 0.1 * nsamples(V4V5_ps0)

# Execute prevalence filter, using `prune_taxa()` function
V4V5_ps0.prev <-  prune_taxa((prev0 > prevalenceThreshold), V4V5_ps0)


#####################################
#Bar plots
#####################################


top50 <- names(sort(taxa_sums(V4V5_ps0.prev), decreasing=TRUE))[1:100]
V4V5_ps0.prop <- transform_sample_counts(V4V5_ps0.prev, function(otu) otu/sum(otu))
ps.top50 <- prune_taxa(top50, V4V5_ps0.prop)
plot_bar(ps.top50,  fill="Class") + facet_wrap(~environmental.package, scales="free_x")


ord.nmds.bray <- ordinate(V4V5_ps0.prop, method="NMDS", distance="bray")

plot_ordination(V4V5_ps0.prop, ord.nmds.bray, label = "geographic.location..depth.", shape="environmental.package", title="Bray NMDS")+ geom_point(size = 5)
