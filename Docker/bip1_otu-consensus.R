# v0.3 unified OTU consensus script

args <- commandArgs(trailingOnly = TRUE)

sampleid <- args[1]
wkdir <- args[2]
stage <- args[3]   # primary or secondary

setwd(wkdir)

library(Biostrings)
library(muscle)
library(DECIPHER)

output <- DNAStringSet()

# ---------------------------
# Configure stage-specific settings
# ---------------------------

if (stage == "primary") {
  
  pattern <- paste0(sampleid, "_OTU.*_Reads.fasta")
  otu_pattern <- "_OTU"
  outfile <- sprintf("%s_OTUs.tmp", sampleid)

} else if (stage == "secondary") {
  
  pattern <- paste0(sampleid, "_FinalOTU.*_Reads.fasta")
  otu_pattern <- "_FinalOTU"
  outfile <- sprintf("%s_finalOTUs.tmp", sampleid)

} else {
  
  stop("stage must be 'primary' or 'secondary'")
  
}

file_list <- grep(pattern, list.files(), value = TRUE)

# ---------------------------
# Process OTU files
# ---------------------------

for (file in file_list) {

  OTUname <- gsub("_Reads.fasta", "", file)
  OTUname <- gsub(otu_pattern, "|OTU", OTUname)

  input <- readDNAStringSet(file)

  # ---------------------------
  # Expand sequences by per-sequence read count (weighted consensus)
  # Handles size=N (primary, from vsearch dereplication) and reads-N (secondary)
  # ---------------------------

  expanded <- DNAStringSet()
  for (i in seq_along(input)) {
    n <- suppressWarnings(
      as.numeric(sub(".*(?:size=|reads-)([0-9]+)$", "\\1", names(input[i]), perl = TRUE))
    )
    if (is.na(n) || n < 1) n <- 1
    expanded <- c(expanded, rep(input[i], n))
  }
  input <- expanded
  len <- length(input)
  readcount <- len

  # ---------------------------
  # Multi-read OTUs
  # ---------------------------

  if (len >= 2) {

    # subsample huge OTUs
    if (len > 1000) {

      indices <- sample(seq_len(len), 1000, replace = FALSE)
      input <- input[indices]
    }
    
    align <- DNAStringSet(muscle(input))
    
    cons <- DNAStringSet(
      ConsensusSequence(
        align,
        includeNonLetters = TRUE,
        includeTerminalGaps = TRUE,
        ambiguity = FALSE,
        noConsensusChar = "N",
        threshold = 0.75,
        minInformation = 0.25
      )
    )
    
    names(cons) <- sprintf("%s|reads-%s", OTUname, readcount)
    
    output <- c(output, cons)
  }
  
  # ---------------------------
  # Single-sequence OTUs
  # ---------------------------
  
  if (len == 1) {
    
    cons <- input
    
    names(cons) <- sprintf("%s|reads-%s", OTUname, readcount)
    
    output <- c(output, cons)
  }
}

# ---------------------------
# Cleanup consensus sequences
# ---------------------------

output <- DNAStringSet(
  gsub("-", "", output)
)

degen <- "[^ATCGN]"

output <- DNAStringSet(
  gsub(degen, "N", output)
)

# ---------------------------
# Write output
# ---------------------------

writeXStringSet(output, outfile)
