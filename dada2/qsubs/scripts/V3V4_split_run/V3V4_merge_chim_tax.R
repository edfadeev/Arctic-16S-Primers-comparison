#load libraries and set random seed
library(ggplot2); packageVersion("ggplot2")
library(dada2); packageVersion("dada2")
library(gridExtra); packageVersion("gridExtra")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
set.seed(123)


seq_tabs <- c(file.path("Seq.Tables","seqtab_mi1.rds"),file.path("Seq.Tables","seqtab_mi2.rds"),
              file.path("Seq.Tables","seqtab_mi3.rds"),file.path("Seq.Tables","seqtab_mi4.rds"),
              file.path("Seq.Tables","seqtab_mi5.rds"),file.path("Seq.Tables","seqtab_mi6.rds"),
              file.path("Seq.Tables","seqtab_mi7.rds"),file.path("Seq.Tables","seqtab_mi8.rds"),
              file.path("Seq.Tables","seqtab_mi9.rds"),file.path("Seq.Tables","seqtab_mi10.rds"),
              file.path("Seq.Tables","seqtab_hi.rds"))

#merge the hiseq and miseq into a single sequence table
seqtab<- mergeSequenceTables(tables= seq_tabs)

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