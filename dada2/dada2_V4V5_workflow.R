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

fnFs_hi<- sort(c("Clipped/42_clip_R1.fastq","Clipped/41_clip_R1.fastq", "Clipped/40_clip_R1.fastq","Clipped/39_clip_R1.fastq","Clipped/38_clip_R1.fastq","Clipped/37_clip_R1.fastq"))
fnRs_hi<- sort(c("Clipped/42_clip_R2.fastq","Clipped/41_clip_R2.fastq", "Clipped/40_clip_R2.fastq","Clipped/39_clip_R2.fastq","Clipped/38_clip_R2.fastq","Clipped/37_clip_R2.fastq"))

fnFs_mi <- sort(fnFs[!fnFs %in%fnFs_hi])
fnRs_mi <- sort(fnRs[!fnRs %in%fnRs_hi]) 
  

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

#separate miseq and hiseq data
filtFs_hi <- filtFs[names(filtFs) %in% c("37","38","39","40","41","42")]
filtRs_hi <- filtRs[names(filtRs) %in% c("37","38","39","40","41","42")]

filtFs_mi <- filtFs[!names(filtFs) %in% c("37","38","39","40","41","42")]
filtRs_mi <- filtRs[!names(filtRs) %in% c("37","38","39","40","41","42")]


#Filter and trim
out_hi <- filterAndTrim(fnFs_hi, filtFs_hi, fnRs_hi, filtRs_hi,
                     maxN=0, maxEE=c(2,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

out_mi <- filterAndTrim(fnFs_mi, filtFs_mi, fnRs_mi, filtRs_mi, truncLen=c(255,200),
                        maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE, verbose = TRUE)

write.csv(rbind(out_hi,out_mi), file= file.path("Report","dada2_filterAndTrim_output.csv"))

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
# Learn errors 
errF_hi <- learnErrors(filtFs_hi, multithread = TRUE)
errR_hi <- learnErrors(filtRs_hi, multithread = TRUE)

errF_mi <- learnErrors(filtFs_mi, multithread = TRUE)
errR_mi <- learnErrors(filtRs_mi, multithread = TRUE)

# Plot error profiles
pdf(file.path("Report","ErrorProfiles.pdf"))
plotErrors(errF_hi, nominalQ = TRUE)
plotErrors(errF_mi, nominalQ = TRUE)
plotErrors(errR_hi, nominalQ = TRUE)
plotErrors(errR_mi, nominalQ = TRUE)
dev.off()

# Sample Inference 
dadaFs_hi <- dada(filtFs_hi, err=errF_hi, multithread=TRUE)
dadaRs_hi <- dada(filtRs_hi, err=errR_hi, multithread=TRUE)

dadaFs_hi[[1]]

dadaFs_mi <- dada(filtFs_mi, err=errF_mi, multithread=TRUE)
dadaRs_mi <- dada(filtRs_mi, err=errR_mi, multithread=TRUE)

dadaFs_mi[[1]]

#Merge paired reads
mergers_hi <- mergePairs(dadaFs_hi, filtFs_hi, dadaRs_hi, filtRs_hi, verbose=TRUE, minOverlap = 10)

mergers_mi <- mergePairs(dadaFs_mi, filtFs_mi, dadaRs_mi, filtRs_mi, verbose=TRUE, minOverlap = 10)

#merge the hiseq and miseq into a single sequence table
seqtab<- mergeSequenceTables(table1= makeSequenceTable(mergers_hi), table2 = makeSequenceTable(mergers_mi))

save.image(file.path("V4V5_dada2.Rdata"))

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

save.image(file.path("V4V5_dada2.Rdata"))


seqs <- getSequences(seqtab.nochim2)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

save.image(file.path("V4V5_dada2.Rdata"))
