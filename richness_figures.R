#set working directory
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#load libraries 
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(iNEXT); packageVersion("iNEXT")
library(cowplot); packageVersion("cowplot")
library(dplyr); packageVersion("dplyr")

#load scripts


#load datasets 
V3V4_data.raw <- readRDS("./Data/V3V4_data_raw.rds")
V4V5_data.raw <- readRDS("./Data/V4V5_data_raw.rds")

#####################################
#Figure - 1 - Workflow reads retainment
#####################################
#V3V4
V3V4_swarm_output <- read.table("./Data/Swarm-output/V3V4/nSeqs_V3V4_all.txt", sep="\t", header = FALSE)
V3V4_swarm_output$SampleID <- V3V4_swarm_output$V1
V3V4_swarm_output$SampleID <- gsub("./Renamed/","",V3V4_swarm_output$SampleID)
V3V4_swarm_output$SampleID <- gsub("_R1.*","",V3V4_swarm_output$SampleID)
V3V4_swarm_output$SampleID <- paste("X",V3V4_swarm_output$SampleID, sep = "")

names(V3V4_swarm_output) <- c("Raw","Clipped","Trimmed","Merged","SampleID")
for (n in c("Raw","Clipped","Trimmed","Merged")){
  V3V4_swarm_output[,c(n)] <- as.numeric(gsub(".*fastq:","",V3V4_swarm_output[,c(n)]))
}
V3V4_swarm_output$Primer <- "V3V4"
V3V4_swarm_output <- V3V4_swarm_output[V3V4_swarm_output$SampleID %in% sample_names(V3V4_data.BAC),]

#calculate proportions for each sample
V3V4_swarm_output.prop <- sweep(V3V4_swarm_output[, -c(5:6)], 1, V3V4_swarm_output[, "Raw"], "/")
V3V4_swarm_output.prop <- cbind(V3V4_swarm_output[, c(5:6)], V3V4_swarm_output.prop)

#V4V5
V4V5_swarm_output <- read.table("./Data/Swarm-output/V4V5/nSeqs_V4V5_all.txt", sep="\t", header = FALSE)
V4V5_swarm_output$SampleID <- V4V5_swarm_output$V1
V4V5_swarm_output$SampleID <- gsub("./Renamed/","",V4V5_swarm_output$SampleID)
V4V5_swarm_output$SampleID <- gsub("_R1.*","",V4V5_swarm_output$SampleID)
V4V5_swarm_output$SampleID <- paste("X",V4V5_swarm_output$SampleID, sep = "")

names(V4V5_swarm_output) <- c("Raw","Clipped","Trimmed","Merged","SampleID")
for (n in c("Raw","Clipped","Trimmed","Merged")){
  V4V5_swarm_output[,c(n)] <- as.numeric(gsub(".*fastq:","",V4V5_swarm_output[,c(n)]))
}
V4V5_swarm_output$Primer <- "V4V5"
V4V5_swarm_output <- V4V5_swarm_output[V4V5_swarm_output$SampleID %in% sample_names(V4V5_data.BAC),]

#calculate proportions for each sample
V4V5_swarm_output.prop <- sweep(V4V5_swarm_output[, -c(5:6)], 1, V4V5_swarm_output[, "Raw"], "/")
V4V5_swarm_output.prop <- cbind(V4V5_swarm_output[, c(5:6)], V4V5_swarm_output.prop)

#merge
swarm_output <- rbind(V3V4_swarm_output.prop,V4V5_swarm_output.prop)
swarm_output <- melt(subset(swarm_output, select = -c(SampleID)),"Primer")
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

#export 
ggsave("./figures/workflow_reads.pdf",
       plot = overview.p,
       scale = 1,
       units = "cm",
       width = 8.6,
       #height = 17.4,
       dpi = 300)
#####################################

#####################################
#Amount of mitochondrial and chloroplast-related sequences
#####################################
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
#Rarefaction curves and Alpha diversity
#####################################
#load datasets 
V3V4_data.BAC <- readRDS("./Data/V3V4_data_BAC.rds")
V4V5_data.BAC <- readRDS("./Data/V4V5_data_BAC.rds")

#calculate coverage
V3V4_data.BAC.iNEXT.out <- iNEXT(as.data.frame(otu_table(V3V4_data.BAC)), q=0, datatype="abundance")
V4V5_data.BAC.iNEXT.out <- iNEXT(as.data.frame(otu_table(V4V5_data.BAC)), q=0, datatype="abundance")

#rarefaction V3V4
V3V4_rare <-fortify(V3V4_data.BAC.iNEXT.out, type=1)
V3V4_meta <- as(sample_data(V3V4_data.BAC), "data.frame")
V3V4_meta$site <- rownames(V3V4_meta)
V3V4_rare$Depth <- V3V4_meta$Depth[match(V3V4_rare$site, V3V4_meta$site)] 
V3V4_rare$SampleType <- V3V4_meta$merging[match(V3V4_rare$site, V3V4_meta$site)] 

V3V4_rare.point <- V3V4_rare[which(V3V4_rare$method == "observed"),]
V3V4_rare.line <- V3V4_rare[which(V3V4_rare$method != "observed"),]
V3V4_rare.line$method <- factor (V3V4_rare.line$method,
                                 c("interpolated", "extrapolated"),
                                 c("interpolation", "extrapolation"))

#plot
V3V4_rare.p <- ggplot()+
  geom_line(aes(x=x, y=y, colour = SampleType, group =site), lwd = 0.5, data= subset(V3V4_rare.line, method == "interpolation"))+
  geom_line(aes(x=x, y=y, colour = SampleType, group =site), linetype =2, lwd = 0.5, data= subset(V3V4_rare.line, method == "extrapolation"))+
  geom_point(aes(x=x, y=y, shape=SampleType), size =2, data= V3V4_rare.point)+
  #geom_text(aes(x=x, y=y, label = Depth), size =3, data= V3V4_rare.point, nudge_y = 500)+
  labs(x = "Sample size", y = "Species richness")+ ylim(0,10000)+xlim(0,2e5)+
  theme_classic(base_size=10)+
  theme(legend.position="none")

#coverage
V3V4_cov <-fortify(V3V4_data.BAC.iNEXT.out, type=3)
V3V4_meta <- as(sample_data(V3V4_data.BAC), "data.frame")
V3V4_meta$site <- rownames(V3V4_meta)
V3V4_cov$Depth <- V3V4_meta$Depth[match(V3V4_cov$site, V3V4_meta$site)] 
V3V4_cov$SampleType <- V3V4_meta$merging[match(V3V4_cov$site, V3V4_meta$site)] 

V3V4_cov.point <- V3V4_cov[which(V3V4_cov$method == "observed"),]
V3V4_cov.line <- V3V4_cov[which(V3V4_cov$method != "observed"),]
V3V4_cov.line$method <- factor (V3V4_cov.line$method,
                                c("interpolated", "extrapolated"),
                                c("interpolation", "extrapolation"))

V3V4_cov.p <- ggplot()+
  geom_line(aes(x=x, y=y, colour = SampleType, group =site), lwd = 0.5, data= subset(V3V4_cov.line, method == "interpolation"))+
  geom_line(aes(x=x, y=y, colour = SampleType, group =site), linetype =2, lwd = 0.5, data= subset(V3V4_cov.line, method == "extrapolation"))+
  geom_point(aes(x=x, y=y, shape = SampleType), size =2, data= V3V4_cov.point)+
  #geom_text(aes(x=x, y=y, label = Depth),size =3, data= V3V4_cov.point)+
  scale_colour_discrete(guide = FALSE)+
  geom_vline(xintercept = 0.95, linetype = "dashed")+ylim(0,10000)+
  labs(x = "Sample size", y = "Species richness")+
  theme_classic(base_size=10)+theme(legend.position="none")

#rarefaction V4V5
V4V5_rare <-fortify(V4V5_data.BAC.iNEXT.out, type=1)
V4V5_meta <- as(sample_data(V4V5_data.BAC), "data.frame")
V4V5_meta$site <- rownames(V4V5_meta)
V4V5_rare$Depth <- V4V5_meta$Depth[match(V4V5_rare$site, V4V5_meta$site)] 
V4V5_rare$SampleType <- V4V5_meta$merging[match(V4V5_rare$site, V4V5_meta$site)] 

V4V5_rare.point <- V4V5_rare[which(V4V5_rare$method == "observed"),]
V4V5_rare.line <- V4V5_rare[which(V4V5_rare$method != "observed"),]
V4V5_rare.line$method <- factor (V4V5_rare.line$method,
                                 c("interpolated", "extrapolated"),
                                 c("interpolation", "extrapolation"))

#plot
V4V5_rare.p <- ggplot()+
  geom_line(aes(x=x, y=y, colour = SampleType, group =site), lwd = 0.5, data= subset(V4V5_rare.line, method == "interpolation"))+
  geom_line(aes(x=x, y=y, colour = SampleType, group =site), linetype =2, lwd = 0.5, data= subset(V4V5_rare.line, method == "extrapolation"))+
  geom_point(aes(x=x, y=y, shape=SampleType), size =2, data= V4V5_rare.point)+
  #geom_text(aes(x=x, y=y, label = Depth), size =3, data= V4V5_rare.point, nudge_y = 500)+
  labs(x = "Sample size", y = "Species richness")+ ylim(0,10000)+xlim(0,2e5)+
  theme_classic(base_size=10)+
  theme(legend.position="none")

#coverage
V4V5_cov <-fortify(V4V5_data.BAC.iNEXT.out, type=3)
V4V5_meta <- as(sample_data(V4V5_data.BAC), "data.frame")
V4V5_meta$site <- rownames(V4V5_meta)
V4V5_cov$Depth <- V4V5_meta$Depth[match(V4V5_cov$site, V4V5_meta$site)] 
V4V5_cov$SampleType <- V4V5_meta$merging[match(V4V5_cov$site, V4V5_meta$site)] 

V4V5_cov.point <- V4V5_cov[which(V4V5_cov$method == "observed"),]
V4V5_cov.line <- V4V5_cov[which(V4V5_cov$method != "observed"),]
V4V5_cov.line$method <- factor (V4V5_cov.line$method,
                                c("interpolated", "extrapolated"),
                                c("interpolation", "extrapolation"))

V4V5_cov.p <- ggplot()+
  geom_line(aes(x=x, y=y, colour = SampleType, group =site), lwd = 0.5, data= subset(V4V5_cov.line, method == "interpolation"))+
  geom_line(aes(x=x, y=y, colour = SampleType, group =site), linetype =2, lwd = 0.5, data= subset(V4V5_cov.line, method == "extrapolation"))+
  geom_point(aes(x=x, y=y, shape = SampleType), size =2, data= V4V5_cov.point)+
  #geom_text(aes(x=x, y=y, label = Depth),size =3, data= V4V5_cov.point)+
  #scale_colour_discrete(guide = FALSE)+
  geom_vline(xintercept = 0.95, linetype = "dashed")+ylim(0,10000)+
  labs(x = "Sample size", y = "Species richness")+
  theme_classic(base_size=10)+theme(legend.position="none")



#combined plot
plot_grid(V3V4_rare.p, V4V5_rare.p, ncol = 2, align = "hv")

ggsave("./figures/rare_curve.pdf",
       plot = ggplot2::last_plot(),
       scale = 1,
       units = "cm",
       width = 17.8,
       #height = 17.4,
       dpi = 300)


#####################################
#Alpha diversity indices
#####################################
#V3V4
#summary table
V3V4_richness <- V3V4_data.BAC.iNEXT.out$AsyEst[V3V4_data.BAC.iNEXT.out$AsyEst$Diversity == "Species richness",]
V3V4_shannon <- V3V4_data.BAC.iNEXT.out$AsyEst[V3V4_data.BAC.iNEXT.out$AsyEst$Diversity == "Shannon diversity",]
V3V4_simpson <- V3V4_data.BAC.iNEXT.out$AsyEst[V3V4_data.BAC.iNEXT.out$AsyEst$Diversity == "Simpson diversity",]

V3V4_comm.char<- data.frame(PANGAEA_sampleID = paste(sample_data(V3V4_data.BAC)$Expedition,sample_data(V3V4_data.BAC)$Station, sep = "/"),
                            #StationName = sample_data(V3V4_data.BAC)$Station,
                            Type = sample_data(V3V4_data.BAC)$Type,
                            Environment = sample_data(V3V4_data.BAC)$merging,
                            Depth = sample_data(V3V4_data.BAC)$Depth,
                            #Sample_sum = V3V4_data.BAC.iNEXT.out$DataInfo$n,
                            Observed = V3V4_data.BAC.iNEXT.out$DataInfo$S.obs,
                            Chao = V3V4_richness$Estimator,
                            Chao.cov = V3V4_richness$Observed/V3V4_richness$Estimator,
                            Shannon = V3V4_shannon$Observed,
                            Shannon.ext = V3V4_shannon$Estimator,
                            Simpson = V3V4_simpson$Observed,
                            Simpson.ext = V3V4_simpson$Estimator,
                            Sample.cov = 100*V3V4_data.BAC.iNEXT.out$DataInfo$SC,
                            Primer="V3V4")


#V4V5
#summary table
V4V5_richness <- V4V5_data.BAC.iNEXT.out$AsyEst[V4V5_data.BAC.iNEXT.out$AsyEst$Diversity == "Species richness",]
V4V5_shannon <- V4V5_data.BAC.iNEXT.out$AsyEst[V4V5_data.BAC.iNEXT.out$AsyEst$Diversity == "Shannon diversity",]
V4V5_simpson <- V4V5_data.BAC.iNEXT.out$AsyEst[V4V5_data.BAC.iNEXT.out$AsyEst$Diversity == "Simpson diversity",]

V4V5_comm.char<- data.frame(PANGAEA_sampleID = paste(sample_data(V4V5_data.BAC)$Expedition,sample_data(V4V5_data.BAC)$Station, sep = "/"),
                            #StationName = sample_data(V4V5_data.BAC)$Station,
                            Type = sample_data(V4V5_data.BAC)$Type,
                            Environment = sample_data(V4V5_data.BAC)$merging,
                            Depth = sample_data(V4V5_data.BAC)$Depth,
                            #Sample_sum = V4V5_data.BAC.iNEXT.out$DataInfo$n,
                            Observed = V4V5_data.BAC.iNEXT.out$DataInfo$S.obs,
                            Chao = V4V5_richness$Estimator,
                            Chao.cov = V4V5_richness$Observed/V4V5_richness$Estimator,
                            Shannon = V4V5_shannon$Observed,
                            Shannon.ext = V4V5_shannon$Estimator,
                            Simpson = V4V5_simpson$Observed,
                            Simpson.ext = V4V5_simpson$Estimator,
                            Sample.cov = 100*V4V5_data.BAC.iNEXT.out$DataInfo$SC,
                            Primer="V4V5")

#####################################
#Alpha diversity significance tests
#####################################
alpha_together <- rbind(V3V4_comm.char,V4V5_comm.char)

df <- alpha_together %>%
  select(Primer, Environment,Chao.cov)

summary(df %>% filter(Primer == "V3V4") %>% .$Chao.cov)
summary(df %>% filter(Primer == "V4V5") %>% .$Chao.cov)

#Chao
t.test(Chao.cov~ Primer,data = alpha_together %>% filter(Environment == "Sea-ice"))
t.test(Chao.cov~ Primer,data = alpha_together %>% filter(Environment == "Sediment"))
t.test(Chao.cov~ Primer,data = alpha_together %>% filter(Environment == "Surface-water"))
t.test(Chao.cov~ Primer,data = alpha_together %>% filter(Environment == "Deep-water"))


t.test(Shannon~ Primer,data = alpha_together %>% filter(Environment == "Sea-ice"))
t.test(Shannon~ Primer,data = alpha_together %>% filter(Environment == "Surface-water"))
t.test(Shannon~ Primer,data = alpha_together %>% filter(Environment == "Deep-water"))
t.test(Shannon~ Primer,data = alpha_together %>% filter(Environment == "Sediment"))

t.test(Simpson~ Primer,data = alpha_together %>% filter(Environment == "Sea-ice"))
t.test(Simpson~ Primer,data = alpha_together %>% filter(Environment == "Surface-water"))
t.test(Simpson~ Primer,data = alpha_together %>% filter(Environment == "Deep-water"))
t.test(Simpson~ Primer,data = alpha_together %>% filter(Environment == "Sediment"))

write.table(alpha_together,"./Data/alpha_div.csv",col.names = TRUE,quote = FALSE, sep = ",")

#####################################
#OTU overlaps
#####################################
library(VennDiagram); packageVersion("VennDiagram")
library(gridExtra); packageVersion("gridExtra")

#calculate overlaps
V3V4_over <- list()
V4V5_over <- list()


for (d in levels(sample_data(V3V4_data.BAC)$merging)){
  sub <- subset_samples(V3V4_data.BAC, merging == d)
  #sub <- prune_taxa(taxa_sums(sub)>0, sub) #remove unobserved
  sub = filter_taxa(sub, function(x) sum(x >= 1) >= (nsamples(sub)/2), TRUE)
  V3V4_over[[d]] <- as.character(row.names(otu_table(sub)))
}

for (d in levels(sample_data(V4V5_data.BAC)$merging)){
  sub <- subset_samples(V4V5_data.BAC, merging == d)
  #sub <- prune_taxa(taxa_sums(sub)>0, sub) #remove unobserved
  sub = filter_taxa(sub, function(x) sum(x >= 1) >= (nsamples(sub)/2), TRUE)
  V4V5_over[[d]] <- as.character(row.names(otu_table(sub)))
}

venn_V3V4 <- venn.diagram(V3V4_over, lwd = 1, fill = c("red","blue","darkblue","gray"),
                          col = c("red", "blue", "green", "orange"), alpha = 0.5, filename = NULL, scaled = TRUE,
                          reverse = TRUE, print.mode  = c("raw","percent"),
                          sigdigs = 1)
venn_V4V5 <- venn.diagram(V4V5_over, lwd = 1, fill = c("red","blue","darkblue","gray"),
                          col = c("red", "blue", "green", "orange"), alpha = 0.5, filename = NULL, scaled = TRUE,
                          reverse = TRUE, print.mode  = c("raw","percent"),
                          sigdigs = 1)
#plot venn diagrams
grid.arrange(grobTree(venn_V3V4), grobTree(venn_V4V5), ncol= 2, nrow= 1)

