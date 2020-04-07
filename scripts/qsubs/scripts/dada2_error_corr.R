#load libraries and set random seed
library(ggplot2); packageVersion("ggplot2")
library(dada2); packageVersion("dada2")
library(gridExtra); packageVersion("gridExtra")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
set.seed(123)

load(file.path("V3V4_dada2.Rdata"))


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


save.image(file.path("V3V4_dada2.Rdata"))


