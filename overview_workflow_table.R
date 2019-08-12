library("dplyr"); packageVersion("dplyr")
library("tidyr"); packageVersion("tidyr")
library("iNEXT"); packageVersion("iNEXT")
library("ggplot2"); packageVersion("ggplot2")
library("cowplot"); packageVersion("cowplot")

#define function for standard error calc. 
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))}

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
swarm_output.agg.Type <- swarm_output.prop%>% group_by(Primer, merging) %>%
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
#Plot rarefaction
####################################
iNEXT.V3V4 <- iNEXT(as.data.frame(otu_table(V3V4_data.BAC)), q=0, datatype="abundance")

rare.V3V4<-fortify(iNEXT.V3V4, type=1)
meta.V3V4<- as(sample_data(V3V4_data.BAC), "data.frame")
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
  geom_point(aes(shape=Type), size =3, data= rare.point.V3V4)+
  #geom_text(aes(label=label), size =2, data= rare.point, colour = "black", nudge_y = -100)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  #xlim(0,2e5)+
  ylim(0,15000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")



iNEXT.V4V5<- iNEXT(as.data.frame(otu_table(V4V5_data.BAC)), q=0, datatype="abundance")

rare.V4V5<-fortify(iNEXT.V4V5, type=1)

meta.V4V5<- as(sample_data(V4V5_data.BAC), "data.frame")
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
  geom_point(aes(shape=Type), size =3, data= rare.point.V4V5)+
  #geom_text(aes(label=label), size =2, data= rare.point, colour = "black", nudge_y = -100)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  #xlim(0,2e5)+
  ylim(0,15000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

plot_grid(V3V4_rare.p,V4V5_rare.p, ncol = 2)

ggsave("./figures/rarefaction.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, 
       height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#Alpha diversity statistical tests
####################################
# Calculate for V3V4
V3V4_data.BAC.div <- estimate_richness(V3V4_data.BAC, split = TRUE, measures = NULL)

#generate data set with all bacterial community characteristics
V3V4_data.BAC_comm.char<- data.frame(Expedition = sample_data(V3V4_data.BAC)$Expedition,
                                     Environment = sample_data(V3V4_data.BAC)$merging,
                                     Depth = sample_data(V3V4_data.BAC)$Depth,
                                     SampleID = sample_data(V3V4_data.BAC)$SampleID,
                                     Primer = "V3V4")

V3V4_data.BAC_comm.char <- left_join(V3V4_data.BAC_comm.char, V3V4_swarm_output[,c("SampleID","Raw","Clipped", "Trimmed", "Merged", "Swarm", "Final")], by = "SampleID")

V3V4_data.BAC_comm.char <- cbind(V3V4_data.BAC_comm.char,data.frame(
                             Observed = V3V4_data.BAC.div$Observed,
                             Chao1 = round(V3V4_data.BAC.div$Chao1,digits=0),
                             Completness = round(100*V3V4_data.BAC.div$Observed/V3V4_data.BAC.div$Chao1, digits=2),
                             Shanonn = round(V3V4_data.BAC.div$Shannon,digits=2),
                             Simpson = round(V3V4_data.BAC.div$Simpson,digits=2),
                             Evenness = round(V3V4_data.BAC.div$Shannon/log(V3V4_data.BAC.div$Observed),digits=2)))


# Calculate for V4V5
V4V5_data.BAC.div <- estimate_richness(V4V5_data.BAC, split = TRUE, measures = NULL)

#generate data set with all bacterial community characteristics
V4V5_data.BAC_comm.char<- data.frame(Expedition = sample_data(V4V5_data.BAC)$Expedition,
                                     Environment = sample_data(V4V5_data.BAC)$merging,
                                     Depth = sample_data(V4V5_data.BAC)$Depth,
                                     SampleID = sample_data(V4V5_data.BAC)$SampleID,
                                     Primer = "V4V5")

V4V5_data.BAC_comm.char <- left_join(V4V5_data.BAC_comm.char, V4V5_swarm_output[,c("SampleID","Raw","Clipped", "Trimmed", "Merged", "Swarm", "Final")], by = "SampleID")

V4V5_data.BAC_comm.char <- cbind(V4V5_data.BAC_comm.char,data.frame(
  Observed = V4V5_data.BAC.div$Observed,
  Chao1 = round(V4V5_data.BAC.div$Chao1,digits=0),
  Completness = round(100*V4V5_data.BAC.div$Observed/V4V5_data.BAC.div$Chao1, digits=2),
  Shanonn = round(V4V5_data.BAC.div$Shannon,digits=2),
  Simpson = round(V4V5_data.BAC.div$Simpson,digits=2),
  Evenness = round(V4V5_data.BAC.div$Shannon/log(V4V5_data.BAC.div$Observed),digits=2)))



#final overview table of all samples
BAC_comm.char <- rbind(V3V4_data.BAC_comm.char,V4V5_data.BAC_comm.char)

BAC_comm.char <- select(BAC_comm.char,-c("SampleID"))

write.csv(BAC_comm.char, "./Data/alpha_table.csv")



