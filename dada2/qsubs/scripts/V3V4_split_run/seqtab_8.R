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
errF_mi7 <- learnErrors(filtFs_mi7, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 20, verbose = TRUE)
errR_mi7 <- learnErrors(filtRs_mi7, multithread = TRUE, randomize = TRUE, MAX_CONSIST = 20, verbose = TRUE)
# Sample Inference 
dadaFs_mi7 <- dada(filtFs_mi7, err=errF_mi7, multithread=TRUE, verbose = TRUE)
dadaRs_mi7 <- dada(filtRs_mi7, err=errR_mi7, multithread=TRUE, verbose = TRUE)
#Merge paired reads
mergers_mi7 <- mergePairs(dadaFs_mi7, filtFs_mi7, dadaRs_mi7, filtRs_mi7, verbose=TRUE, minOverlap = 10)
#generate sequence table and save it
seqtab_mi7 <- makeSequenceTable(mergers_mi7)
saveRDS(seqtab_mi7, file.path("Seq.Tables","seqtab_mi7.rds"))
