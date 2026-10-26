# v0.4 unified OTU consensus script (MAFFT, distinct reads, abundance-weighted consensus)

args <- commandArgs(trailingOnly = TRUE)

sampleid <- args[1]
wkdir <- args[2]
stage <- args[3]   # primary or secondary

setwd(wkdir)

library(Biostrings)
library(DECIPHER)

# Align with MAFFT (default FFT-NS-2 strategy). One thread: BIP already runs one bip1 per
# sample in parallel. Returns the alignment in the same order as the input.
mafft_align <- function(seqs) {
  names(seqs) <- paste0("s", seq_along(seqs))
  infile <- tempfile(fileext = ".fasta"); outfile <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(infile, outfile)))
  writeXStringSet(seqs, infile)
  status <- system2("mafft", c("--quiet", "--thread", "1", infile), stdout = outfile)
  if (!identical(status, 0L) || file.size(outfile) == 0) {
    stop("MAFFT failed (exit status ", status, "). Is mafft installed and on PATH?")
  }
  aln <- readDNAStringSet(outfile)
  aln <- DNAStringSet(toupper(as.character(aln)))
  aln[names(seqs)]
}

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
  # Read count behind each distinct sequence
  # Handles size=N (primary, from vsearch dereplication) and reads-N (secondary)
  # ---------------------------

  weights <- suppressWarnings(
    as.numeric(sub(".*(?:size=|reads-)([0-9]+)$", "\\1", names(input), perl = TRUE))
  )
  weights[is.na(weights) | weights < 1] <- 1
  readcount <- sum(weights)

  # ---------------------------
  # Multi-read OTUs
  # ---------------------------

  if (readcount >= 2) {

    # Subsample huge OTUs to 1000 reads, drawn in proportion to abundance as before, then
    # keep each distinct sequence once with the number of draws it received. Seeded so the
    # same input always gives the same consensus.
    if (readcount > 1000) {
      set.seed(1)
      draws <- tabulate(sample(rep(seq_along(input), weights), 1000), nbins = length(input))
      keep <- draws > 0
      input <- input[keep]
      weights <- draws[keep]
    }

    # Align each distinct sequence once, then weight its aligned row by its read count, so
    # the consensus is the same as aligning every read but far cheaper when reads repeat.
    if (length(input) == 1) {
      align <- input
    } else {
      align <- mafft_align(input)
    }
    align <- align[rep(seq_along(align), weights)]

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
  
  if (readcount == 1) {
    
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
