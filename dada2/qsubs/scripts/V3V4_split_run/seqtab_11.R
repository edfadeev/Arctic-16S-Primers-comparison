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

# Make directory and filenames for the filtered fastqs
filt_path <- file.path("Filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filt_path <- file.path("Seq.Tables")
if(!file_test("-d", filt_path)) dir.create(filt_path)

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
errF_mi10 <- learnErrors(filtFs_mi10, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 20, verbose = TRUE)
errR_mi10 <- learnErrors(filtRs_mi10, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 20, verbose = TRUE)
# Sample Inference 
dadaFs_mi10 <- dada(filtFs_mi10, err=errF_mi10, multithread=TRUE,verbose = TRUE)
dadaRs_mi10 <- dada(filtRs_mi10, err=errR_mi10, multithread=TRUE,verbose = TRUE)
#Merge paired reads
mergers_mi10 <- mergePairs(dadaFs_mi10, filtFs_mi10, dadaRs_mi10, filtRs_mi10, verbose=TRUE, minOverlap = 10)
#generate sequence table and save it
seqtab_mi10 <- makeSequenceTable(mergers_mi10)
saveRDS(seqtab_mi10, file.path("Seq.Tables","seqtab_mi10.rds"))

