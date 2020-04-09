#load libraries and set random seed
library(ggplot2); packageVersion("ggplot2")
library(dada2); packageVersion("dada2")
library(gridExtra); packageVersion("gridExtra")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
set.seed(123)


# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files("Clipped", pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files("Clipped", pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_clip_R1.fastq"), `[`, 1)

filt_path <- file.path("Report")
if(!file_test("-d", filt_path)) dir.create(filt_path)

# quality check
QualityProfileFs <- list()
for(i in 1:length(fnFs)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(fnFs[i])
}
pdf(file.path("Report","RawProfileForward.pdf"))
for(i in 1:length(fnFs)) {
  do.call("grid.arrange", QualityProfileFs[[i]])  
}
dev.off()
rm(QualityProfileFs)

QualityProfileRs <- list()
for(i in 1:length(fnRs)) {
  QualityProfileRs[[i]] <- list()
  QualityProfileRs[[i]][[1]] <- plotQualityProfile(fnRs[i])
}
pdf(file.path("Report","RawProfileReverse.pdf"))
for(i in 1:length(fnRs)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)

# Make directory and filenames for the filtered fastqs
filt_path <- file.path("Filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path("Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#separate the different runs
#Hiseq
fnFs_hi<- sort(file.path("Clipped",paste(c(22,23,26,27), "_clip_R1.fastq", sep = "")))
fnRs_hi<- sort(file.path("Clipped",paste(c(22,23,26,27), "_clip_R2.fastq", sep = ""))) 
filtFs_hi<- sort(file.path("Filtered",paste(c(22,23,26,27), "_F_filt.fastq.gz", sep = "")))
filtRs_hi<- sort(file.path("Filtered",paste(c(22,23,26,27), "_R_filt.fastq.gz", sep = ""))) 
#Filter and trim
out_hi <- filterAndTrim(fnFs_hi, filtFs_hi, fnRs_hi, filtRs_hi,
                        maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE)
# Learn errors 
errF_hi <- learnErrors(filtFs_hi, multithread = TRUE)
errR_hi <- learnErrors(filtRs_hi, multithread = TRUE)
# Sample Inference 
dadaFs_hi <- dada(filtFs_hi, err=errF_hi, multithread=TRUE)
dadaRs_hi <- dada(filtRs_hi, err=errR_hi, multithread=TRUE)
#Merge paired reads
mergers_hi <- mergePairs(dadaFs_hi, filtFs_hi, dadaRs_hi, filtRs_hi, verbose=TRUE, minOverlap = 10)

#MiSeq 1 - 493:307
fnFs_mi1 <- sort(file.path("Clipped",paste(c(11,12,14:21), "_clip_R1.fastq", sep = "")))
fnRs_mi1 <- sort(file.path("Clipped",paste(c(11,12,14:21), "_clip_R2.fastq", sep = "")))
filtFs_mi1 <- sort(file.path("Filtered",paste(c(11,12,14:21), "_F_filt.fastq.gz", sep = "")))
filtRs_mi1 <- sort(file.path("Filtered",paste(c(11,12,14:21), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi1 <- filterAndTrim(fnFs_mi1, filtFs_mi1, fnRs_mi1, filtRs_mi1, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi1 <- learnErrors(filtFs_mi1, multithread = TRUE)
errR_mi1 <- learnErrors(filtRs_mi1, multithread = TRUE)
# Sample Inference 
dadaFs_mi1 <- dada(filtFs_mi1, err=errF_mi1, multithread=TRUE)
dadaRs_mi1 <- dada(filtRs_mi1, err=errR_mi1, multithread=TRUE)
#Merge paired reads
mergers_mi1 <- mergePairs(dadaFs_mi1, filtFs_mi1, dadaRs_mi1, filtRs_mi1, verbose=TRUE, minOverlap = 10)

#MiSeq 2 - 493:227
fnFs_mi2 <- sort(file.path("Clipped",paste(c(24,25), "_clip_R1.fastq", sep = "")))
fnRs_mi2 <- sort(file.path("Clipped",paste(c(24,25), "_clip_R2.fastq", sep = "")))
filtFs_mi2 <- sort(file.path("Filtered",paste(c(24,25), "_F_filt.fastq.gz", sep = "")))
filtRs_mi2 <- sort(file.path("Filtered",paste(c(24,25), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi2 <- filterAndTrim(fnFs_mi2, filtFs_mi2, fnRs_mi2, filtRs_mi2, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi2 <- learnErrors(filtFs_mi2, multithread = TRUE)
errR_mi2 <- learnErrors(filtRs_mi2, multithread = TRUE)
# Sample Inference 
dadaFs_mi2 <- dada(filtFs_mi2, err=errF_mi2, multithread=TRUE)
dadaRs_mi2 <- dada(filtRs_mi2, err=errR_mi2, multithread=TRUE)
#Merge paired reads
mergers_mi2 <- mergePairs(dadaFs_mi2, filtFs_mi2, dadaRs_mi2, filtRs_mi2, verbose=TRUE, minOverlap = 10)

#MiSeq 3 - 493:291
fnFs_mi3 <- sort(file.path("Clipped",paste(c(28:31), "_clip_R1.fastq", sep = "")))
fnRs_mi3 <- sort(file.path("Clipped",paste(c(28:31), "_clip_R2.fastq", sep = "")))
filtFs_mi3 <- sort(file.path("Filtered",paste(c(28:31), "_F_filt.fastq.gz", sep = "")))
filtRs_mi3 <- sort(file.path("Filtered",paste(c(28:31), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi3 <- filterAndTrim(fnFs_mi3, filtFs_mi3, fnRs_mi3, filtRs_mi3, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi3 <- learnErrors(filtFs_mi3, multithread = TRUE)
errR_mi3 <- learnErrors(filtRs_mi3, multithread = TRUE)
# Sample Inference 
dadaFs_mi3 <- dada(filtFs_mi3, err=errF_mi3, multithread=TRUE)
dadaRs_mi3 <- dada(filtRs_mi3, err=errR_mi3, multithread=TRUE)
#Merge paired reads
mergers_mi3 <- mergePairs(dadaFs_mi3, filtFs_mi3, dadaRs_mi3, filtRs_mi3, verbose=TRUE, minOverlap = 10)

#MiSeq 4 - 493:255
fnFs_mi4 <- sort(file.path("Clipped",paste(c(1,2,9), "_clip_R1.fastq", sep = "")))
fnRs_mi4 <- sort(file.path("Clipped",paste(c(1,2,9), "_clip_R2.fastq", sep = "")))
filtFs_mi4 <- sort(file.path("Filtered",paste(c(1,2,9), "_F_filt.fastq.gz", sep = "")))
filtRs_mi4 <- sort(file.path("Filtered",paste(c(1,2,9), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi4 <- filterAndTrim(fnFs_mi4, filtFs_mi4, fnRs_mi4, filtRs_mi4, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi4 <- learnErrors(filtFs_mi4, multithread = TRUE)
errR_mi4 <- learnErrors(filtRs_mi4, multithread = TRUE)
# Sample Inference 
dadaFs_mi4 <- dada(filtFs_mi4, err=errF_mi4, multithread=TRUE)
dadaRs_mi4 <- dada(filtRs_mi4, err=errR_mi4, multithread=TRUE)
#Merge paired reads
mergers_mi4 <- mergePairs(dadaFs_mi4, filtFs_mi4, dadaRs_mi4, filtRs_mi4, verbose=TRUE, minOverlap = 10)

#MiSeq 5 - 493:295
fnFs_mi5 <- sort(file.path("Clipped",paste(c(32), "_clip_R1.fastq", sep = "")))
fnRs_mi5 <- sort(file.path("Clipped",paste(c(32), "_clip_R2.fastq", sep = "")))
filtFs_mi5 <- sort(file.path("Filtered",paste(c(32), "_F_filt.fastq.gz", sep = "")))
filtRs_mi5 <- sort(file.path("Filtered",paste(c(32), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi5 <- filterAndTrim(fnFs_mi5, filtFs_mi5, fnRs_mi5, filtRs_mi5, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi5 <- learnErrors(filtFs_mi5, multithread = TRUE)
errR_mi5 <- learnErrors(filtRs_mi5, multithread = TRUE)
# Sample Inference 
dadaFs_mi5 <- dada(filtFs_mi5, err=errF_mi5, multithread=TRUE)
dadaRs_mi5 <- dada(filtRs_mi5, err=errR_mi5, multithread=TRUE)
#Merge paired reads
mergers_mi5 <- mergePairs(dadaFs_mi5, filtFs_mi5, dadaRs_mi5, filtRs_mi5, verbose=TRUE, minOverlap = 10)

#MiSeq 6 - 493:293
fnFs_mi6 <- sort(file.path("Clipped",paste(c(33:35), "_clip_R1.fastq", sep = "")))
fnRs_mi6 <- sort(file.path("Clipped",paste(c(33:35), "_clip_R2.fastq", sep = "")))
filtFs_mi6 <- sort(file.path("Filtered",paste(c(33:35), "_F_filt.fastq.gz", sep = "")))
filtRs_mi6 <- sort(file.path("Filtered",paste(c(33:35), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi6 <- filterAndTrim(fnFs_mi6, filtFs_mi6, fnRs_mi6, filtRs_mi6, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi6 <- learnErrors(filtFs_mi6, multithread = TRUE)
errR_mi6 <- learnErrors(filtRs_mi6, multithread = TRUE)
# Sample Inference 
dadaFs_mi6 <- dada(filtFs_mi6, err=errF_mi6, multithread=TRUE)
dadaRs_mi6 <- dada(filtRs_mi6, err=errR_mi6, multithread=TRUE)
#Merge paired reads
mergers_mi6 <- mergePairs(dadaFs_mi6, filtFs_mi6, dadaRs_mi6, filtRs_mi6, verbose=TRUE, minOverlap = 10)

#MiSeq 7 - MISEQ:226
fnFs_mi7 <- sort(file.path("Clipped",paste(c(36,38), "_clip_R1.fastq", sep = "")))
fnRs_mi7 <- sort(file.path("Clipped",paste(c(36,38), "_clip_R2.fastq", sep = "")))
filtFs_mi7 <- sort(file.path("Filtered",paste(c(36,38), "_F_filt.fastq.gz", sep = "")))
filtRs_mi7 <- sort(file.path("Filtered",paste(c(36,38), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi7 <- filterAndTrim(fnFs_mi7, filtFs_mi7, fnRs_mi7, filtRs_mi7, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi7 <- learnErrors(filtFs_mi7, multithread = TRUE)
errR_mi7 <- learnErrors(filtRs_mi7, multithread = TRUE)
# Sample Inference 
dadaFs_mi7 <- dada(filtFs_mi7, err=errF_mi7, multithread=TRUE)
dadaRs_mi7 <- dada(filtRs_mi7, err=errR_mi7, multithread=TRUE)
#Merge paired reads
mergers_mi7 <- mergePairs(dadaFs_mi7, filtFs_mi7, dadaRs_mi7, filtRs_mi7, verbose=TRUE, minOverlap = 10)

#MiSeq 8 - MISEQ:221
fnFs_mi8 <- sort(file.path("Clipped",paste(c(37,39:42), "_clip_R1.fastq", sep = "")))
fnRs_mi8 <- sort(file.path("Clipped",paste(c(37,39:42), "_clip_R2.fastq", sep = "")))
filtFs_mi8 <- sort(file.path("Filtered",paste(c(37,39:42), "_F_filt.fastq.gz", sep = "")))
filtRs_mi8 <- sort(file.path("Filtered",paste(c(37,39:42), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi8 <- filterAndTrim(fnFs_mi8, filtFs_mi8, fnRs_mi8, filtRs_mi8, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi8 <- learnErrors(filtFs_mi8, multithread = TRUE)
errR_mi8 <- learnErrors(filtRs_mi8, multithread = TRUE)
# Sample Inference 
dadaFs_mi8 <- dada(filtFs_mi8, err=errF_mi8, multithread=TRUE)
dadaRs_mi8 <- dada(filtRs_mi8, err=errR_mi8, multithread=TRUE)
#Merge paired reads
mergers_mi8 <- mergePairs(dadaFs_mi8, filtFs_mi8, dadaRs_mi8, filtRs_mi8, verbose=TRUE, minOverlap = 10)

#MiSeq 9 - 493:307
fnFs_mi9 <- sort(file.path("Clipped",paste(c(7,8,10), "_clip_R1.fastq", sep = "")))
fnRs_mi9 <- sort(file.path("Clipped",paste(c(7,8,10), "_clip_R2.fastq", sep = "")))
filtFs_mi9 <- sort(file.path("Filtered",paste(c(7,8,10), "_F_filt.fastq.gz", sep = "")))
filtRs_mi9 <- sort(file.path("Filtered",paste(c(7,8,10), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi9 <- filterAndTrim(fnFs_mi9, filtFs_mi9, fnRs_mi9, filtRs_mi9, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi9 <- learnErrors(filtFs_mi9, multithread = TRUE)
errR_mi9 <- learnErrors(filtRs_mi9, multithread = TRUE)
# Sample Inference 
dadaFs_mi9 <- dada(filtFs_mi9, err=errF_mi9, multithread=TRUE)
dadaRs_mi9 <- dada(filtRs_mi9, err=errR_mi9, multithread=TRUE)
#Merge paired reads
mergers_mi9 <- mergePairs(dadaFs_mi9, filtFs_mi9, dadaRs_mi9, filtRs_mi9, verbose=TRUE, minOverlap = 10)

#MiSeq 10 - 493:471
fnFs_mi10 <- sort(file.path("Clipped",paste(c(3:6), "_clip_R1.fastq", sep = "")))
fnRs_mi10 <- sort(file.path("Clipped",paste(c(3:6), "_clip_R2.fastq", sep = "")))
filtFs_mi10 <- sort(file.path("Filtered",paste(c(3:6), "_F_filt.fastq.gz", sep = "")))
filtRs_mi10 <- sort(file.path("Filtered",paste(c(3:6), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi10 <- filterAndTrim(fnFs_mi10, filtFs_mi10, fnRs_mi10, filtRs_mi10, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi10 <- learnErrors(filtFs_mi10, multithread = TRUE)
errR_mi10 <- learnErrors(filtRs_mi10, multithread = TRUE)
# Sample Inference 
dadaFs_mi10 <- dada(filtFs_mi10, err=errF_mi10, multithread=TRUE)
dadaRs_mi10 <- dada(filtRs_mi10, err=errR_mi10, multithread=TRUE)
#Merge paired reads
mergers_mi10 <- mergePairs(dadaFs_mi10, filtFs_mi10, dadaRs_mi10, filtRs_mi10, verbose=TRUE, minOverlap = 10)

##summary of filtering

# quality check
QualityProfileFs <- list()
for(i in 1:length(filtFs)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(filtFs[i])
}
pdf(file.path("Report","FiltProfileForward.pdf"))
for(i in 1:length(filtFs)) {
  do.call("grid.arrange", QualityProfileFs[[i]])  
}
dev.off()
rm(QualityProfileFs)

QualityProfileRs <- list()
for(i in 1:length(filtRs)) {
  QualityProfileRs[[i]] <- list()
  QualityProfileRs[[i]][[1]] <- plotQualityProfile(filtRs[i])
}
pdf(file.path("Report","FiltProfileReverse.pdf"))
for(i in 1:length(filtRs)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)

# Plot error profiles
pdf(file.path("Report","ErrorProfiles.pdf"))
plotErrors(errF_hi, nominalQ = TRUE)+ggtitle("HiSeq- 168:111")
plotErrors(errF_mi1, nominalQ = TRUE)+ggtitle("MiSeq- 493:307")
plotErrors(errF_mi2, nominalQ = TRUE)+ggtitle("MiSeq- 493:227")
plotErrors(errF_mi3, nominalQ = TRUE)+ggtitle("MiSeq- 493:291")
plotErrors(errF_mi4, nominalQ = TRUE)+ggtitle("MiSeq- 493:255")
plotErrors(errF_mi5, nominalQ = TRUE)+ggtitle("MiSeq- 493:295")
plotErrors(errF_mi6, nominalQ = TRUE)+ggtitle("MiSeq- 493:293")
plotErrors(errF_mi7, nominalQ = TRUE)+ggtitle("MiSeq- MISEQ:226")
plotErrors(errF_mi8, nominalQ = TRUE)+ggtitle("MiSeq- MISEQ:221")
plotErrors(errF_mi9, nominalQ = TRUE)+ggtitle("MiSeq- 493:307")
plotErrors(errF_mi10, nominalQ = TRUE)+ggtitle("MiSeq- MISEQ:471")
plotErrors(errR_hi, nominalQ = TRUE)+ggtitle("HiSeq- 168:111")
plotErrors(errR_mi1, nominalQ = TRUE)+ggtitle("MiSeq- 493:307")
plotErrors(errR_mi2, nominalQ = TRUE)+ggtitle("MiSeq- 493:227")
plotErrors(errR_mi3, nominalQ = TRUE)+ggtitle("MiSeq- 493:291")
plotErrors(errR_mi4, nominalQ = TRUE)+ggtitle("MiSeq- 493:255")
plotErrors(errR_mi5, nominalQ = TRUE)+ggtitle("MiSeq- 493:295")
plotErrors(errR_mi6, nominalQ = TRUE)+ggtitle("MiSeq- 493:293")
plotErrors(errR_mi7, nominalQ = TRUE)+ggtitle("MiSeq- MISEQ:226")
plotErrors(errR_mi8, nominalQ = TRUE)+ggtitle("MiSeq- MISEQ:221")
plotErrors(errR_mi9, nominalQ = TRUE)+ggtitle("MiSeq- 493:307")
plotErrors(errR_mi10, nominalQ = TRUE)+ggtitle("MiSeq- MISEQ:471")
dev.off()


#write out filtered read counts
write.csv(rbind(out_hi,out_mi1,out_mi2,out_mi3,out_mi4,out_mi5,out_mi6,
                out_mi7,out_mi8,
                out_mi9,out_mi10), 
          file= file.path("Report","dada2_filterAndTrim_output.csv"))



#merge the hiseq and miseq into a single sequence table
seqtab<- mergeSequenceTables(table1= makeSequenceTable(mergers_hi), 
                             table2 = makeSequenceTable(mergers_mi1),
                             table3 = makeSequenceTable(mergers_mi2),
                             table4 = makeSequenceTable(mergers_mi3),
                             table5 = makeSequenceTable(mergers_mi4),
                             table6 = makeSequenceTable(mergers_mi5),
                             table7 = makeSequenceTable(mergers_mi6),
                             table8 = makeSequenceTable(mergers_mi7),
                             table9 = makeSequenceTable(mergers_mi8),
                             table10 = makeSequenceTable(mergers_mi9),
                             table11 = makeSequenceTable(mergers_mi10))

save.image("V3V4_dada2_sep_runs.Rdata")

#Combine together sequences that are identical 
seqtab1 <- collapseNoMismatch(seqtab, verbose = TRUE)

dim(seqtab1)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab1)))

seqtab.nochim <- removeBimeraDenovo(seqtab1, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#proportion of chimeras
sum(seqtab.nochim)/sum(seqtab1)

# inspect output: remove singletons and 'junk' sequences
# read lengths modified for V34 amplicons / based upon output table where majority of reads occurs
seqtab.nochim2 <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% c(370:430) & colSums(seqtab.nochim) > 1]
dim(seqtab.nochim2)
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "tabled")
rownames(track) <- sample.names

write.csv(track, file.path(path, "Report","dada2_reads_output.csv"))

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim2, "../tax/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE, verbose = TRUE)
taxa <- addSpecies(taxa, "../tax/silva_species_assignment_v138.fa.gz", tryRC = TRUE, verbose = TRUE)

save.image("V3V4_dada2_sep_runs.Rdata")
