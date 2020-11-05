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
filt_path <- file.path("Seq.Tables")
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
                        compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors 
errF_hi <- learnErrors(filtFs_hi, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_hi <- learnErrors(filtRs_hi, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_hi <- dada(filtFs_hi, err=errF_hi, multithread=TRUE, verbose = TRUE)
dadaRs_hi <- dada(filtRs_hi, err=errR_hi, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_hi <- mergePairs(dadaFs_hi, filtFs_hi, dadaRs_hi, filtRs_hi, verbose=TRUE, minOverlap = 10)
#generate sequence table and save it
seqtab_hi <- makeSequenceTable(mergers_hi)
saveRDS(seqtab_hi, file.path("Seq.Tables","seqtab_hi.rds"))

#MiSeq 1 - 493:307
fnFs_mi1 <- sort(file.path("Clipped",paste(c(1:21,24:25,28:37), "_clip_R1.fastq", sep = "")))
fnRs_mi1 <- sort(file.path("Clipped",paste(c(1:21,24:25,28:37), "_clip_R2.fastq", sep = "")))
filtFs_mi1 <- sort(file.path("Filtered",paste(c(1:21,24:25,28:37), "_F_filt.fastq.gz", sep = "")))
filtRs_mi1 <- sort(file.path("Filtered",paste(c(1:21,24:25,28:37), "_R_filt.fastq.gz", sep = "")))
#Filter and trim
out_mi1 <- filterAndTrim(fnFs_mi1, filtFs_mi1, fnRs_mi1, filtRs_mi1, truncLen=c(255,200),
                         maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE, verbose = TRUE)
# Learn errors
errF_mi1 <- learnErrors(filtFs_mi1, nbases = 2e+08, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
errR_mi1 <- learnErrors(filtRs_mi1, nbases = 2e+08, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 30, verbose = TRUE)
# Sample Inference 
dadaFs_mi1 <- dada(filtFs_mi1, err=errF_mi1, multithread=TRUE, verbose = TRUE)
dadaRs_mi1 <- dada(filtRs_mi1, err=errR_mi1, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi1 <- mergePairs(dadaFs_mi1, filtFs_mi1, dadaRs_mi1, filtRs_mi1, verbose=TRUE, minOverlap = 10)
#generate sequence table and save it
seqtab_mi1 <- makeSequenceTable(mergers_mi1)
saveRDS(seqtab_mi1, file.path("Seq.Tables","seqtab_mi1.rds"))


##summary of filtering
save.image("V3V4_dada2_sep_runs.Rdata")
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
plotErrors(errF_hi, nominalQ = TRUE)+ggtitle("HiSeq")
plotErrors(errF_mi1, nominalQ = TRUE)+ggtitle("MiSeq")
plotErrors(errR_hi, nominalQ = TRUE)+ggtitle("HiSeq")
plotErrors(errR_mi1, nominalQ = TRUE)+ggtitle("MiSeq")
dev.off()

#write out filtered read counts
write.csv(rbind(out_hi,out_mi1), 
          file= file.path("Report","dada2_filterAndTrim_output.csv"))



#merge the hiseq and miseq into a single sequence table
seqtab<- mergeSequenceTables(table1= makeSequenceTable(mergers_hi), 
                             table2 = makeSequenceTable(mergers_mi1))

save.image("V3V4_dada2_merged_runs.Rdata")

#Combine together sequences that are identical 
seqtab1 <- collapseNoMismatch(seqtab, verbose = TRUE)

dim(seqtab1)

save.image("V3V4_dada2_sep_runs.Rdata")

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

#correct sample name in the table
rownames(seqtab.nochim)[25]<- "32_F_filt.fastq.gz"
rownames(seqtab.nochim2)[25]<- "32_F_filt.fastq.gz"

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim2, "../tax/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE, verbose = TRUE)
taxa <- addSpecies(taxa, "../tax/silva_species_assignment_v138.fa.gz", tryRC = TRUE, verbose = TRUE)

save.image("V3V4_dada2_sep_runs.Rdata")

# get summary tables 
sample.order <- names(c(dadaFs_mi1,dadaFs_hi))

getN <- function(x) sum(getUniques(x))
track <- cbind(rbind(out_mi1,out_hi), 
               sapply(c(dadaFs_mi1,dadaFs_hi), getN),
               sapply(c(mergers_mi1,mergers_hi), getN),
               rowSums(seqtab.nochim[sample.order,]),
               rowSums(seqtab.nochim2[sample.order,]))
colnames(track) <- c("input", "filtered", "denoised", "merged", "nochim", "tabled")
rownames(track) <- gsub("_F_filt.fastq.gz","",sample.order)
track <- data.frame(track)

#add unclassified levels of taxonomy 
TAX <- taxa
k <- ncol(TAX) - 1
for (i in 2:k) {
  if (sum(is.na(TAX[, i])) > 1) {
    test <- TAX[is.na(TAX[, i]), ]
    for (j in 1:nrow(test)) {
      if (sum(is.na(test[j, i:(k + 1)])) == length(test[j, i:(k + 1)])) {
        test[j, i] <- paste(test[j, (i - 1)], "_uncl", sep = "")
        test[j, (i + 1):(k + 1)] <- test[j, i]
      }
    }
    TAX[is.na(TAX[, i]), ] <- test
  }
  if (sum(is.na(TAX[, i])) == 1) {
    test <- TAX[is.na(TAX[, i]), ]
    if (sum(is.na(test[i:(k + 1)])) == length(test[i:(k + 1)])) {
      test[i] <- paste(test[(i - 1)], "_uncl", sep = "")
      test[(i + 1):(k + 1)] <- test[i]
    }
    TAX[is.na(TAX[, i]),] <- test
  }
}
TAX[is.na(TAX[, (k + 1)]), (k + 1)] <- paste(TAX[is.na(TAX[, (k + 1)]), k], "_uncl", sep = "")

# write output
write.table(t(seqtab.nochim2), "dada2_seqtab_nochim2.txt", quote = F, sep = "\t")
write.table(TAX, "dada2_taxonomy_table.txt", sep = "\t", quote = F)
write.table(track, "libs_summary_table.txt", sep = "\t", quote = F)
