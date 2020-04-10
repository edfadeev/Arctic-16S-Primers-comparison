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
errF_mi5 <- learnErrors(filtFs_mi5, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 20, verbose = TRUE)
errR_mi5 <- learnErrors(filtRs_mi5, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 20, verbose = TRUE)
# Sample Inference 
dadaFs_mi5 <- dada(filtFs_mi5, err=errF_mi5, multithread=TRUE, verbose = TRUE)
dadaRs_mi5 <- dada(filtRs_mi5, err=errR_mi5, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi5 <- mergePairs(dadaFs_mi5, filtFs_mi5, dadaRs_mi5, filtRs_mi5, verbose=TRUE, minOverlap = 10)
#generate sequence table and save it
seqtab_mi5 <- makeSequenceTable(list(sample1=mergers_mi5))
saveRDS(seqtab_mi5, file.path("Seq.Tables","seqtab_mi5.rds"))

