#install missing packages
.cran_packages <- c("ggplot2", "gridExtra","BiocManager")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(.bioc_packages[!.inst], ask = F)    
}

#load libraries and set random seed
library(ggplot2); packageVersion("ggplot2")
library(dada2); packageVersion("dada2")
library(gridExtra); packageVersion("gridExtra")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
set.seed(123)

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(file.path("Clipped"), pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(file.path("Clipped"), pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_clip_R1.fastq"), `[`, 1)

filt_path <- file.path(path, "Report")
if(!file_test("-d", filt_path)) dir.create(filt_path)

# quality check
QualityProfileFs <- list()
for(i in 1:length(fnFs)) {
  QualityProfileFs[[i]] <- list()
  QualityProfileFs[[i]][[1]] <- plotQualityProfile(fnFs[i])
}
pdf(file.path("Report","QualityProfileForward.pdf"))
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
pdf(file.path("Report","QualityProfileReverse.pdf"))
for(i in 1:length(fnRs)) {
  do.call("grid.arrange", QualityProfileRs[[i]])  
}
dev.off()
rm(QualityProfileRs)

# Place filtered files in filtered/ subdirectory
# Make directory and filenames for the filtered fastqs
filt_path <- file.path("Filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path("Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=10, truncLen=c(240,180),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)


write.csv(out, file= file.path("Report","dada2_filterAndTrim_output.csv"))

# Quality check 
QualityProfileFs.filt <- list()
for(i in 1:length(filtFs)) {
  QualityProfileFs.filt[[i]] <- list()
  QualityProfileFs.filt[[i]][[1]] <- plotQualityProfile(filtFs[i])
}
pdf(file.path("Report","QualityProfileForwardFiltered.pdf"))
for(i in 1:length(filtFs)) {
  do.call("grid.arrange", QualityProfileFs.filt[[i]])  
}
dev.off()
rm(QualityProfileFs.filt)

QualityProfileRs.filt <- list()
for(i in 1:length(filtRs)) {
  QualityProfileRs.filt[[i]] <- list()
  QualityProfileRs.filt[[i]][[1]] <- plotQualityProfile(filtRs[i])
}
pdf(file.path("Report","QualityProfileReverseFiltered.pdf"))
for(i in 1:length(filtRs)) {
  do.call("grid.arrange", QualityProfileRs.filt[[i]])  
}
dev.off()
rm(QualityProfileRs.filt)
# Learn errors 
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Plot error profiles
pdf(file.path("Report","ErrorProfiles.pdf"))
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Sample Inference 
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 10)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#proportion of chimeras
sum(seqtab.nochim)/sum(seqtab)

# inspect output: remove singletons and 'junk' sequences
# read lengths modified for V34 amplicons / based upon output table where majority of reads occurs
seqtab.nochim2 <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% c(360:400) & colSums(seqtab.nochim) > 1]
dim(seqtab.nochim2)
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "tabled")
rownames(track) <- sample.names

write.csv(track, file.path("Report","dada2_reads_output.csv"))

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim2, "../tax/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)
taxa <- addSpecies(taxa, "../tax/silva_species_assignment_v138.fa.gz")

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

# format output
seqtab.nochim.print <- t(seqtab.nochim2)
tax.print.arch <- taxa
all.equal(rownames(seqtab.nochim.print), rownames(taxa))  #TRUE
rownames(seqtab.nochim.print) <- paste("sq", 1:ncol(seqtab.nochim.print), sep = "")
rownames(taxa) <- rownames(seqtab.nochim.print)

# write output
write.table(seqtab.nochim.print, file.path("V3V4_seqtab_nochim.txt"), quote = F, sep = "\t")
write.table(taxa, file.path("V3V4_taxonomy_table.txt"), sep = "\t", quote = F)

save.image("dada2_V3V4.Rdata")
