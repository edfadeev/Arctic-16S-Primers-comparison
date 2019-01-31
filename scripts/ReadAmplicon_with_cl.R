#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements:
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# wrapper to read the output of AmpliconNGSworkflow.txt and to output curated sample by OTU and OTU by taxonomy tables
#
# input: 
# otu - file name of OTU_contingency_table.csv
#       columns: OTU (OTU number), amplicon (seed sequence header), sequence counts per sample, total (total sequence abundance per OTU over all samples)
# tax - file name of amplicons_seeds_taxonomy.txt
#       columns: accnos (seed sequence header), align (alignment quality), path (taxonomic path)
# silva - file name of Silva taxmap
#         columns: path (taxonomic path), node (name of last element of path), rank (corresponding taxonomic rank)
# domain - Which domain was sequenced? (Bacteria, Archaea, Eukaryota)
# silva.sep - separator of Silva taxmap
# singletons - Does data set include singeltons (default TRUE)
# unclassified - In case of unclassified taxonomic levels, 
#                should "_unclassified" be appended to the last classified taxonomic level? (default TRUE)
# write.files - Should output tables be written to file? (default TRUE)
#
# output:
# Sample by OTU and OTU by taxonomy table (may need further manual curation), 
# alignment quality for representative sequence per OTU, original seed sequence accession number and original Silva taxonomic path
#
# dependencies:
# function optimized for output of AmpliconNGSworkflow.txt
#=============================================================================
#documentation end


ReadAmplicon_with_cl <- function(otu, tax, silva, domain,
                         silva.sep = "\t", singletons = T, unclassified = T, write.files = T) {
  
  #reading otu table
  Data0 <- read.table(otu, h = T, sep = "\t")
  colnames(Data0) <- sapply(colnames(Data0), function(x) {
    strsplit(x, "_der", fixed = T)[[1]][1]
  }
  )
  rownames(Data0) <- paste("otu", Data0$OTU, sep = "")
  
  if (singletons == FALSE) {
    Data0 <- Data0[Data0$total > 1, ]
  }
  
  #reading taxonomic information
  Tax0 <- read.table(tax, h = F, sep = "\t")
  colnames(Tax0) <- c("accnos", "align", "path")
  
  #check order of sequences
  accnos <- sapply(as.character(Tax0$accnos), function(x) {
    strsplit(x, "_", fixed = T)[[1]][1]
  }
  )
  names(accnos) <- NULL
  if (all.equal(as.character(accnos), as.character(Data0$amplicon)) != TRUE) {
    print("OTU_contigency_table and taxonomy do not match\n")
  } else {
    Tax0$accnos <- accnos
    rownames(Tax0) <- rownames(Data0)
    
    #removing unwanted lineages
    if (domain == "Bacteria") {
      Data0a <- Data0[grep("Bacteria", Tax0$path), ]
      Tax0a <- Tax0[grep("Bacteria", Tax0$path), ]
      Data0b <- Data0a#[-c(grep("chloroplast", Tax0a$path), grep("mitochondria", Tax0a$path)), ]
      Tax0b <- Tax0a#[-c(grep("chloroplast", Tax0a$path), grep("mitochondria", Tax0a$path)), ]
    }
    if (domain == "Archaea") {
      Data0b <- Data0[grep("Archaea", Tax0$path), ]
      Tax0b <- Tax0[grep("Archaea", Tax0$path), ]
    }
    if (domain == "Eukaryota") {
      Data0b <- Data0[grep("Eukaryota", Tax0$path), ]
      Tax0b <- Tax0[grep("Eukaryota", Tax0$path), ]
    }
    
    #parsing taxonomic path
    if (domain %in% c("Bacteria", "Archaea")) {
      SILVAtaxopath <- function(tax, SILVA) {
        # create matrix with nrow=number of OTUs and ncol=number of possible taxonomic levels
        output <- matrix(NA, nrow = length(tax), ncol = length(levels(SILVA$rank)))
        colnames(output) <- levels(SILVA$rank)
        
        for (i in 1:length(tax)) {
          if (length(tax[[i]]) == 6){
            output[i, c("domain","phylum","class","order","family","genus")] <- tax[[i]]
          } else {
            # this loop goes through the taxonomic path of an OTU one level at a time and assigns the name of each level to its correct rank according to SILVA
            for (j in 1:length(levels(SILVA$rank))) {
              if (paste(tax[[i]][1:j], collapse=";") %in% SILVA$path) {
                output[i, as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "rank"])] <- as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "node"])
              }
            }
          }
        }
        return(output)
      }  
    }
    if (domain == "Eukaryota") {
      SILVAtaxopath <- function(tax, SILVA){
        output <- matrix(NA, nrow = length(tax), ncol = length(levels(SILVA$rank)))
        colnames(output) <- levels(SILVA$rank)
        for (i in 1:length(tax)) {
          for (j in 1:length(levels(SILVA$rank))) {
            if (paste(tax[[i]][1:j], collapse = ";") %in% SILVA$path) {
              output[i, as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "rank"])] <- as.character(SILVA[SILVA$path == paste(tax[[i]][1:j], collapse = ";"), "node"])
            }
          }
        }
        return(output)
      }
    }
    
    SILVA0 <- read.table(silva, sep = silva.sep, h=T)
    SILVA <- SILVA0[grep(domain, SILVA0$path), ]
    Tax0b_split <- strsplit(as.character(Tax0b$path), ";")
    Tax <- SILVAtaxopath(Tax0b_split, SILVA)
    if (domain %in% c("Bacteria", "Archaea")) {
      #select only the ranks of interest
      Taxa <- Tax[,c("domain","phylum","class","order","family","genus")] 
    }
    if (domain == "Eukaroyta") {
      Taxa <- Tax[,c("domain","major_clade","kingdom","phylum","class","order","family","genus")] 
    }
    rownames(Taxa) <- rownames(Data0b)
    Taxb <- Taxa
    
    if (unclassified == TRUE) {
      Taxb[Taxb == "uncultured"] <- NA
      k <- ncol(Taxb) - 1
      for (i in 2:k) {
        if (sum(is.na(Taxb[, i])) > 1) {
          test <- Taxb[is.na(Taxb[, i]), ]
          for (j in 1:nrow(test)) {
            if (sum(is.na(test[j, i:(k + 1)])) == length(test[j, i:(k + 1)])) {
              test[j, i] <- paste(test[j, (i - 1)], "_unclassified", sep = "")
              test[j, (i + 1):(k + 1)] <- test[j, i]
            }
          }
          Taxb[is.na(Taxb[, i]), ] <- test
        }
        if (sum(is.na(Taxb[, i])) == 1) {
          test <- Taxb[is.na(Taxb[, i]), ]
          if (sum(is.na(test[i:(k + 1)])) == length(test[i:(k + 1)])) {
            test[i] <- paste(test[(i - 1)], "_unclassified", sep = "")
            test[(i + 1):(k + 1)] <- test[i]
          }
          Taxb[is.na(Taxb[, i]),] <- test
        }
      }
      Taxb[is.na(Taxb[, (k + 1)]), (k + 1)] <- paste(Taxb[is.na(Taxb[, (k + 1)]), k], "_unclassified", sep = "")
    }
    
    align_qual <- Tax0b[rownames(Data0b), c("accnos","align","path")]
    
    output <- list(OTU = Data0b[, 3:(ncol(Data0b) - 1)], TAX = Taxb, ALIGN = align_qual)
    if (write.files == TRUE) {
      write.table(output$OTU, "OTU_table.txt", sep = "\t", quote = F)
      write.table(output$TAX, "Taxonomy_table.txt", sep = "\t", quote = F)
      write.table(output$ALIGN, "Align_qual.txt", sep = "\t", quote = F)
    }
    return(output)
  } 
}
