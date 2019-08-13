library("dplyr"); packageVersion("dplyr")
library("tidyr"); packageVersion("tidyr")
library("iNEXT"); packageVersion("iNEXT")
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")

#define function for standard error calc. 
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))}


#load phyloseq objects
V3V4_data.raw <- readRDS("./Data/V3V4_data_raw.rds")
V4V5_data.raw <- readRDS("./Data/V4V5_data_raw.rds")
V3V4_data.BAC <- readRDS("./Data/V3V4_data_BAC.rds")
V4V5_data.BAC <- readRDS("./Data/V4V5_data_BAC.rds")

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
V3V4_swarm_output$Swarm <- sample_sums(V3V4_data.raw)
V3V4_swarm_output$Final <- sample_sums(V3V4_data.BAC)
V3V4_swarm_output$Primer <- "V3V4"

V3V4_swarm_output <- left_join(as(sample_data(V3V4_data.BAC),"data.frame")[,c("Type","merging","SampleID")],V3V4_swarm_output, by = "SampleID")

#V4V5
#Read counts from the QC
V4V5_swarm_output <- read.table("./Data/Swarm-output/V4V5/nSeqs_V4V5_all.txt", sep="\t", header = FALSE)

#remove paths of files
V4V5_swarm_output$SampleID <- V4V5_swarm_output$V1
V4V5_swarm_output$SampleID <- gsub("./Renamed/","",V4V5_swarm_output$SampleID)
V4V5_swarm_output$SampleID <- gsub("_R1.*","",V4V5_swarm_output$SampleID)
V4V5_swarm_output$SampleID <- paste("X",V4V5_swarm_output$SampleID, sep = "")

names(V4V5_swarm_output) <- c("Raw","Clipped","Trimmed","Merged", "SampleID")
for (n in c("Raw","Clipped","Trimmed","Merged")){
  V4V5_swarm_output[,c(n)] <- as.numeric(gsub(".*fastq:","",V4V5_swarm_output[,c(n)]))
}

#add final number of sequences
V4V5_swarm_output <- V4V5_swarm_output[V4V5_swarm_output$SampleID %in% sample_names(V4V5_data.BAC),]
V4V5_swarm_output$Swarm <- sample_sums(V4V5_data.raw)
V4V5_swarm_output$Final <- sample_sums(V4V5_data.BAC)
V4V5_swarm_output$Primer <- "V4V5"

V4V5_swarm_output <- left_join(as(sample_data(V4V5_data.BAC),"data.frame")[,c("Type","merging","SampleID")],V4V5_swarm_output, by = "SampleID")

#merge into one dataframe
swarm_output.both <- rbind(V3V4_swarm_output,V4V5_swarm_output)

swarm_output.prop <- swarm_output.both[,c("Type","merging","Primer","SampleID","Raw","Clipped","Trimmed","Merged","Swarm","Final")] %>%
                          mutate_at(vars(c("Raw","Clipped","Trimmed","Merged","Swarm","Final")), funs(. / swarm_output.both$Raw))

#####################################
#summarize the workflow 
#####################################
swarm_output.agg.Type <- swarm_output.prop%>% group_by(Primer) %>%
                          summarize_at(c("Raw","Clipped","Trimmed","Merged","Swarm", "Final"), c("mean","se"))



#Counts of mitochondrial and chloroplast-related sequences
V3V4_mit_chl_counts <- data.frame(Primer = "V3V4",Station = sample_data(V3V4_data.raw)$Station,Type = sample_data(V3V4_data.raw)$Type,Depth = sample_data(V3V4_data.raw)$Depth, 
                                  Mitochondria = 1-sample_sums(subset_taxa(V3V4_data.raw, family !="Mitochondria"))/sample_sums(V3V4_data.raw), 
                                  Chloroplast = 1-sample_sums(subset_taxa(V3V4_data.raw, order !="Chloroplast"))/sample_sums(V3V4_data.raw))

V4V5_mit_chl_counts <- data.frame(Primer = "V4V5", Station = sample_data(V4V5_data.raw)$Station,Type = sample_data(V4V5_data.raw)$Type,Depth = sample_data(V4V5_data.raw)$Depth, 
                                  Mitochondria = 1-sample_sums(subset_taxa(V4V5_data.raw, family !="Mitochondria"))/sample_sums(V4V5_data.raw), 
                                  Chloroplast = 1-sample_sums(subset_taxa(V4V5_data.raw, order !="Chloroplast"))/sample_sums(V4V5_data.raw))

#summarize for both primers
mit.agg <- as.data.frame(as.list(aggregate(Mitochondria~Primer,rbind(V3V4_mit_chl_counts[V3V4_mit_chl_counts$Depth<50,],V4V5_mit_chl_counts[V4V5_mit_chl_counts$Depth<50,]),FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))))
chl.agg <- as.data.frame(as.list(aggregate(Chloroplast~Primer,rbind(V3V4_mit_chl_counts[V3V4_mit_chl_counts$Depth<50,],V4V5_mit_chl_counts[V4V5_mit_chl_counts$Depth<50,]),FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))))


#####################################
# Compute prevalence of each feature
#####################################
#V3V4
prevdf <-  apply(X = otu_table(V3V4_data.BAC),
                 MARGIN = ifelse(taxa_are_rows(V3V4_data.BAC), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# # Add taxonomy and total read counts to this data.frame
prevdf.tax  <-  data.frame(Prevalence = prevdf,
                           TotalAbundance = taxa_sums(V3V4_data.BAC),
                           tax_table(V3V4_data.BAC))
#summarize
prevdf.tax.summary <- plyr::ddply(prevdf.tax, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# 
# #plot
prev_plot_phyl <- ggplot(prevdf.tax, aes(TotalAbundance, Prevalence / nsamples(V3V4_data.BAC),color=phylum)) +
  # # Include a guess for parameter
  geom_hline(yintercept = 0.15, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold <- round(0.15 * nsamples(V3V4_data.BAC))
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
V3V4_data.BAC.prev <-  prune_taxa((prevdf > prevalenceThreshold), V3V4_data.BAC)

#V4V5
prevdf <-  apply(X = otu_table(V4V5_data.BAC),
                 MARGIN = ifelse(taxa_are_rows(V4V5_data.BAC), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# # Add taxonomy and total read counts to this data.frame
prevdf.tax  <-  data.frame(Prevalence = prevdf,
                           TotalAbundance = taxa_sums(V4V5_data.BAC),
                           tax_table(V4V5_data.BAC))
#summarize
prevdf.tax.summary <- plyr::ddply(prevdf.tax, "phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# 
# #plot
prev_plot_phyl <- ggplot(prevdf.tax, aes(TotalAbundance, Prevalence / nsamples(V4V5_data.BAC),color=phylum)) +
  # # Include a guess for parameter
  geom_hline(yintercept = 0.15, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold <- round(0.15* nsamples(V4V5_data.BAC))
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
V4V5_data.BAC.prev <-  prune_taxa((prevdf > prevalenceThreshold), V4V5_data.BAC)

#####################################
#Plot rarefaction
####################################
iNEXT.V3V4 <- iNEXT(as.data.frame(otu_table(V3V4_data.BAC.prev)), q=0, datatype="abundance")

rare.V3V4<-fortify(iNEXT.V3V4, type=1)
meta.V3V4<- as(sample_data(V3V4_data.BAC.prev), "data.frame")
meta.V3V4$site <- rownames(meta.V3V4)
rare.V3V4$Community <- meta.V3V4$Community[match(rare.V3V4$site, meta.V3V4$site)] 
rare.V3V4$StationName <- meta.V3V4$StationName[match(rare.V3V4$site, meta.V3V4$site)] 
rare.V3V4$Type <- meta.V3V4$merging[match(rare.V3V4$site, meta.V3V4$site)]
rare.V3V4$Type <- factor(rare.V3V4$Type, levels = c("Sea-ice", "Surface-water", "Deep-water", "Sediment" ))
  
rare.V3V4$label <- paste(rare.V3V4$StationName,rare.V3V4$Type, sep = "-")

rare.point.V3V4<- rare.V3V4[which(rare.V3V4$method == "observed"),]
rare.line.V3V4 <- rare.V3V4[which(rare.V3V4$method != "observed"),]
rare.line.V3V4$method <- factor(rare.line.V3V4$method,
                            c("interpolated", "extrapolated"),
                            c("interpolation", "extrapolation"))


V3V4_rare.p <- ggplot(rare.V3V4, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= rare.line.V3V4)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(aes(shape=Type), size =3, data= rare.point.V3V4, colour = "black")+
  #geom_text(aes(label=label), size =2, data= rare.point, colour = "black", nudge_y = -100)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  #xlim(0,2e5)+
  ylim(0,15000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")



iNEXT.V4V5<- iNEXT(as.data.frame(otu_table(V4V5_data.BAC.prev)), q=0, datatype="abundance")

rare.V4V5<-fortify(iNEXT.V4V5, type=1)

meta.V4V5<- as(sample_data(V4V5_data.BAC.prev), "data.frame")
meta.V4V5$site <- rownames(meta.V4V5)
rare.V4V5$Community <- meta.V4V5$Community[match(rare.V4V5$site, meta.V4V5$site)] 
rare.V4V5$StationName <- meta.V4V5$StationName[match(rare.V4V5$site, meta.V4V5$site)] 
rare.V4V5$Type <- meta.V4V5$merging[match(rare.V4V5$site, meta.V4V5$site)]
rare.V4V5$Type <- factor(rare.V4V5$Type, levels = c("Sea-ice", "Surface-water", "Deep-water", "Sediment" ))

rare.V4V5$label <- paste(rare.V4V5$StationName,rare.V4V5$Type, sep = "-")

rare.point.V4V5<- rare.V4V5[which(rare.V4V5$method == "observed"),]
rare.line.V4V5 <- rare.V4V5[which(rare.V4V5$method != "observed"),]
rare.line.V4V5$method <- factor(rare.line.V4V5$method,
                                c("interpolated", "extrapolated"),
                                c("interpolation", "extrapolation"))


V4V5_rare.p <- ggplot(rare.V4V5, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= rare.line.V4V5)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(aes(shape=Type), size =3, data= rare.point.V4V5, colour = "black")+
  #geom_text(aes(label=label), size =2, data= rare.point, colour = "black", nudge_y = -100)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  #xlim(0,2e5)+
  ylim(0,15000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

plot_grid(V3V4_rare.p,V4V5_rare.p, ncol = 2)

ggsave("./figures/rarefaction-prev.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, 
       height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#Alpha diversity table
####################################
# Calculate for V3V4
V3V4_data.BAC.prev.div <- estimate_richness(V3V4_data.BAC.prev, split = TRUE, measures = NULL)

#generate data set with all bacterial community characteristics
V3V4_data.BAC.prev_comm.char<- data.frame(Expedition = sample_data(V3V4_data.BAC.prev)$Expedition,
                                     Station = sample_data(V3V4_data.BAC.prev)$Station,
                                     Environment = sample_data(V3V4_data.BAC.prev)$merging,
                                     Depth = sample_data(V3V4_data.BAC.prev)$Depth,
                                     SampleID = sample_data(V3V4_data.BAC.prev)$SampleID,
                                     Primer = "V3V4")

V3V4_data.BAC.prev_comm.char <- left_join(V3V4_data.BAC.prev_comm.char, V3V4_swarm_output[,c("SampleID","Raw","Clipped", "Trimmed", "Merged", "Swarm", "Final")], by = "SampleID")

V3V4_data.BAC.prev_comm.char <- cbind(V3V4_data.BAC.prev_comm.char,data.frame(
                             Observed = V3V4_data.BAC.prev.div$Observed,
                             Chao1 = round(V3V4_data.BAC.prev.div$Chao1,digits=0),
                             Completness = round(100*V3V4_data.BAC.prev.div$Observed/V3V4_data.BAC.prev.div$Chao1, digits=2),
                             Shanonn = round(V3V4_data.BAC.prev.div$Shannon,digits=2),
                             Simpson = round(V3V4_data.BAC.prev.div$Simpson,digits=2),
                             Evenness = round(V3V4_data.BAC.prev.div$Shannon/log(V3V4_data.BAC.prev.div$Observed),digits=2)))


# Calculate for V4V5
V4V5_data.BAC.prev.div <- estimate_richness(V4V5_data.BAC.prev, split = TRUE, measures = NULL)

#generate data set with all bacterial community characteristics
V4V5_data.BAC.prev_comm.char<- data.frame(Expedition = sample_data(V4V5_data.BAC.prev)$Expedition,
                                     Station = sample_data(V4V5_data.BAC.prev)$Station,
                                     Environment = sample_data(V4V5_data.BAC.prev)$merging,
                                     Depth = sample_data(V4V5_data.BAC.prev)$Depth,
                                     SampleID = sample_data(V4V5_data.BAC.prev)$SampleID,
                                     Primer = "V4V5")

V4V5_data.BAC.prev_comm.char <- left_join(V4V5_data.BAC.prev_comm.char, V4V5_swarm_output[,c("SampleID","Raw","Clipped", "Trimmed", "Merged", "Swarm", "Final")], by = "SampleID")

V4V5_data.BAC.prev_comm.char <- cbind(V4V5_data.BAC.prev_comm.char,data.frame(
  Observed = V4V5_data.BAC.prev.div$Observed,
  Chao1 = round(V4V5_data.BAC.prev.div$Chao1,digits=0),
  Completness = round(100*V4V5_data.BAC.prev.div$Observed/V4V5_data.BAC.prev.div$Chao1, digits=2),
  Shanonn = round(V4V5_data.BAC.prev.div$Shannon,digits=2),
  Simpson = round(V4V5_data.BAC.prev.div$Simpson,digits=2),
  Evenness = round(V4V5_data.BAC.prev.div$Shannon/log(V4V5_data.BAC.prev.div$Observed),digits=2)))

#final overview table of all samples
BAC_comm.char <- rbind(V3V4_data.BAC.prev_comm.char,V4V5_data.BAC.prev_comm.char)

BAC_comm.char <- select(BAC_comm.char,-c("SampleID"))

write.csv(BAC_comm.char, "./Data/alpha_table.csv")


#coverage based on Chao1
BAC_comm.char.cov.agg <- group_by(BAC_comm.char, Primer) %>%
                          summarize(mean.cov = mean(Completness),
                                    se.cov = se(Completness))


#####################################
#Richness statistics
####################################
BAC_comm.char$sample <- paste(BAC_comm.char$Environment,BAC_comm.char$Station,BAC_comm.char$Depth, sep = " ")
              

select(BAC_comm.char,c("sample", "Primer","Chao1")) %>% 
group_by(sample, Primer, Chao1) %>% 
  spread(Primer, Chao1) %>% 
  separate(sample, c("Environment", "Station", "Depth"), " ") %>% 
  group_by(Environment) %>% 
  do(tidy(t.test(.$V3V4, 
                 .$V4V5, 
                 mu = 0, 
                 alt = "two.sided", 
                 paired = TRUE, 
                 conf.level = 0.99)))
  

    
select(BAC_comm.char,c("sample", "Primer","Observed")) %>% 
  group_by(sample, Primer, Observed) %>% 
  spread(Primer, Observed) %>% 
  separate(sample, c("Environment", "Station", "Depth"), " ") %>% 
  group_by(Environment) %>% 
  do(tidy(t.test(.$V3V4, 
                 .$V4V5, 
                 mu = 0, 
                 alt = "two.sided", 
                 paired = TRUE, 
                 conf.level = 0.99)))

