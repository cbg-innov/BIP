# v0.1

#get argument from bash script (plate name)
args <- commandArgs(trailingOnly = TRUE)
runid <- args[1]
wkdir <- args[2]
minreads <- args[3]
minreads <- as.numeric(minreads)
minreads <- round(minreads, digits = 0)
marker <- args[4]
if(marker == "COI-5P"){
  marker <- TRUE
}
primers <- args[5]

#set working directory based on plate's folder, load required libraries
setwd(wkdir)
library("Biostrings")
library("readxl")
library(muscle)
library("dplyr")
library("msa")
library(stringr)
library(DECIPHER)
library(readr)

# input OTU consensus sequences and tax ID results
input.fasta <- readDNAStringSet(sprintf("%s.corrected.fasta", runid))
input.table <- read.table(sprintf("%s.table", runid), header = TRUE, fill = TRUE, sep = "\t", row.names = NULL)
input.table <- read_tsv(sprintf("%s.table", runid), col_names = TRUE)
if(nrow(input.table) == 0) {
  input.table <- data.frame("SeqName" = names(input.fasta),
                            "TaxAssign" = "unknown",
                            "Strand" = "unknown",
                            "TaxAssignFinal" = "k:unknown;p:unknown;c:unknown;o:unknown;f:unknown;g:unknown;s:unknown")
}

##### make OTU table #####
names <- names(input.fasta)
seqs <- as.character(input.fasta, use.names = FALSE)
df <- data.frame(Sequence = seqs, SeqName = names, stringsAsFactors = FALSE)
parsed <- strsplit(df$SeqName, "\\|")
df$Sample <- sub("__.*$", "", sapply(parsed, "[[", 1))
df$OTUName <- sapply(parsed, "[[", 2)
df$ReadCount <- sapply(parsed, "[[", 3)
df$STOPCODON <- sapply(parsed, function(x) "STOP" %in% x)
df$ReadCount <- as.numeric(as.character(gsub("reads-","",df$ReadCount)))
df <- df[order(df$Sample, -df$ReadCount),]

output <- merge(df, input.table, by = "SeqName", all.x = TRUE)
# Map SINTAX rank prefixes to column names — shared by both parse functions below
rank_prefix_map <- c(k = "Kingdom", p = "Phylum", c = "Class",
                     o = "Order",   f = "Family", g = "Genus", s = "Species")

# Determine which ranks exist in the reference library by scanning raw SINTAX results.
# Ranks absent from the library get "N/A"; ranks present but unassigned get "unknown".
all_prefixes_seen <- unique(unlist(regmatches(
  output$TaxAssign[!is.na(output$TaxAssign)],
  gregexpr("[a-z]+(?=:)", output$TaxAssign[!is.na(output$TaxAssign)], perl = TRUE)
)))
ranks_in_library <- unname(rank_prefix_map[names(rank_prefix_map) %in% all_prefixes_seen])
all_rank_names   <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_defaults     <- setNames(
  ifelse(all_rank_names %in% ranks_in_library, "unknown", "N/A"),
  all_rank_names
)

parse_tax_final <- function(tax_string) {
  result <- tax_defaults
  if (is.na(tax_string) || trimws(tax_string) == "") return(result)
  for (part in unlist(strsplit(tax_string, "[,;]"))) {
    part   <- trimws(part)
    prefix <- sub(":.*", "", part)
    value  <- sub("^[^:]+:", "", part)
    rank   <- rank_prefix_map[prefix]
    if (!is.na(rank)) result[rank] <- value
  }
  result
}

tax_cols <- as.data.frame(
  t(sapply(output$TaxAssignFinal, parse_tax_final)),
  stringsAsFactors = FALSE
)
output <- cbind(output, tax_cols)
output$RunID <- runid

#add sequence length to output table
output$SeqLength <- nchar(output$Sequence)

# add number Ns to output table
output$Ns <- nchar(output$Sequence) - nchar(gsub("N", "", output$Sequence))

#generate tax table
input.tax <- read.delim("../parameters.tsv",
                         header = TRUE,
                         sep = "\t",
                         na.strings = "NA",
                         stringsAsFactors = FALSE,
                         check.names = FALSE)
# parameters.tsv holds one row per WELL and covers EVERY primer pair, so a sample
# run under two markers appears more than once. Merging that straight onto the OTU
# table fans out - every OTU is emitted once per matching row - which double-counts
# reads and makes the duplicated rows tie on ReadCount, demoting genuine Dominant
# OTUs to "Tied" and losing their barcodes. Restrict to this primer pair's rows and
# keep one row per sample before merging.
input.tax <- input.tax[
  paste(input.tax[["Forward Primer Name"]],
        input.tax[["Reverse Primer Name"]], sep = "_") == primers, ]
input.tax <- input.tax[!duplicated(input.tax$Sample), ]

input.tax.table <- input.tax[,c("Sample","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
names(input.tax.table) <- c("Sample","USER.Kingdom","USER.Phylum","USER.Class","USER.Order","USER.Family","USER.Genus","USER.Species")
input.tax.table$USER.Kingdom[is.na(input.tax.table$USER.Kingdom) | input.tax.table$USER.Kingdom == ""] <- "unknown"
input.tax.table$USER.Phylum[is.na(input.tax.table$USER.Phylum)   | input.tax.table$USER.Phylum   == ""] <- "unknown"
input.tax.table$USER.Class[is.na(input.tax.table$USER.Class)     | input.tax.table$USER.Class     == ""] <- "unknown"
input.tax.table$USER.Order[is.na(input.tax.table$USER.Order)     | input.tax.table$USER.Order     == ""] <- "unknown"
input.tax.table$USER.Family[is.na(input.tax.table$USER.Family)   | input.tax.table$USER.Family   == ""] <- "unknown"
input.tax.table$USER.Genus[is.na(input.tax.table$USER.Genus)     | input.tax.table$USER.Genus     == ""] <- "unknown"
input.tax.table$USER.Species[is.na(input.tax.table$USER.Species) | input.tax.table$USER.Species   == ""] <- "unknown"

#add USER Phylum, Class, and Order ID to OTU table
output$merge_key <- output$Sample
output <- merge(output, input.tax.table, by.x = "merge_key", by.y = "Sample", all.x = TRUE)
output$merge_key <- NULL  # drop the temp column after merging

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

known.non.targets <- c("Bacteria", "Fungi", "Protista", "Nematoda", "Heterokontophyta", "Rhodophyta")

#flag OTUs that have taxonomic mismatches or are unknown
compute_rank_match <- function(seq_rank, user_rank) {
  ifelse(
    is.na(seq_rank) | is.na(user_rank) | seq_rank == "unknown" | user_rank == "unknown" | seq_rank == "N/A",
    "UNKNOWN",
    ifelse(seq_rank == user_rank, "OK", "MISMATCH")
  )
}

# Apply to each rank
output$Kingdom_Match <- compute_rank_match(output$Kingdom, output$USER.Kingdom)
output$Phylum_Match  <- compute_rank_match(output$Phylum, output$USER.Phylum)
output$Class_Match   <- compute_rank_match(output$Class, output$USER.Class)
output$Order_Match   <- compute_rank_match(output$Order, output$USER.Order)
output$Family_Match  <- compute_rank_match(output$Family, output$USER.Family)
output$Genus_Match   <- compute_rank_match(output$Genus, output$USER.Genus)
output$Species_Match <- compute_rank_match(output$Species, output$USER.Species)

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_match_cols <- paste0(ranks, "_Match")

# Cascade: if any higher rank is MISMATCH, force all lower ranks to MISMATCH
for (i in seq_along(rank_match_cols)) {
  if (i == 1) next
  higher_cols <- rank_match_cols[1:(i - 1)]
  output[[rank_match_cols[i]]] <- ifelse(
    apply(output[, higher_cols, drop = FALSE] == "MISMATCH", 1, any),
    "MISMATCH",
    output[[rank_match_cols[i]]]
  )
}

output[, rank_match_cols] <- lapply(output[, rank_match_cols], function(col) {
  replace(col, is.na(col), "UNKNOWN")
})

# Tax_Match based on high ranks only (Kingdom–Order); Family/Genus/Species → LowRankScore
high_cols <- c("Kingdom_Match", "Phylum_Match", "Class_Match", "Order_Match")
output$Tax_Match <- ifelse(
  apply(output[, high_cols] == "MISMATCH", 1, any), "MISMATCH",
  ifelse(
    apply(output[, high_cols] == "UNKNOWN", 1, any), "UNKNOWN",
    "OK"
  )
)

low_cols <- c("Family_Match", "Genus_Match", "Species_Match")
output$LowRankScore <- rowSums(
  sapply(output[, low_cols], function(x) ifelse(x == "OK", 1L, 0L)),
  na.rm = TRUE
)

# add dominant tag to output table
output <- output[order(output$Sample, -output$ReadCount),]
output$OTURank <- "Secondary"
topotu <- which(!duplicated(output$Sample))
output$OTURank[topotu] <- "Dominant"

# add tied tag to output table
for(i in 1:(nrow(output)-1)){
  if(output$OTURank[i] == "Dominant" && output$Sample[i] == output$Sample[i+1] && output$ReadCount[i] == output$ReadCount[i+1]) {
    output$OTURank[i] <- "Tied"
    output$OTURank[i+1] <- "Tied"
  }
}

# add INDEL tag to table (default to FALSE if not COI-658)
expected_lengths <- c(640, 643, 646, 649, 652, 655, 658, 661, 664, 667, 670)
ifelse(marker=="TRUE",
       output$INDEL <- ifelse(output$SeqLength %in% expected_lengths, FALSE, TRUE),
       output$INDEL <- FALSE)
output$HMM <- ifelse(output$INDEL == FALSE & output$STOPCODON == FALSE, "OK", "HMM_ISSUE")

# split SINTAX results (with probabilities) into their own columns by rank prefix
parse_sintax <- function(tax_string) {
  ranks_all <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  out_names <- c(rbind(paste0("SINTAX_", ranks_all), paste0("SINTAX_", ranks_all, "_Prob")))
  result <- setNames(rep(NA_character_, length(out_names)), out_names)
  if (is.na(tax_string) || trimws(tax_string) == "") return(result)
  for (part in unlist(strsplit(tax_string, ","))) {
    part   <- trimws(part)
    prefix <- sub(":.*", "", part)
    rank   <- rank_prefix_map[prefix]
    if (is.na(rank)) next
    inner  <- sub("^[^:]+:", "", part)
    value  <- sub("\\(.*", "", inner)
    prob_m <- regmatches(inner, regexpr("(?<=\\().*(?=\\))", inner, perl = TRUE))
    prob   <- if (length(prob_m) == 1) prob_m else NA_character_
    result[paste0("SINTAX_", rank)]          <- value
    result[paste0("SINTAX_", rank, "_Prob")] <- prob
  }
  result
}

df_split <- as.data.frame(
  t(sapply(output$TaxAssign, parse_sintax)),
  stringsAsFactors = FALSE
)

# Combine the split data back into the original data frame
output <- cbind(output, df_split)

# add OTU category to table
output$OTUCategory <- if(marker == "TRUE"){
                            ifelse(output$SeqLength < 640, "NUMT",
                              ifelse(output$SeqLength > 670, "NUMT",
                                ifelse(output$INDEL == FALSE & output$STOPCODON == TRUE, "NUMT",
                                    ifelse(output$HMM != "HMM_ISSUE", "mtCOI", "unknown"))))
} else{
  output$OTUCategory <- "unknown"
}
                            
output$OTUCategory[is.na(output$OTUCategory)] <- "unknown"

# remove secondary OTUs with fewer than minreads that also match dominant OTU taxonomy
output <- data.frame(output %>%
                       group_by(Sample) %>%
                       mutate(
                         should_remove = OTURank == "Secondary" & ReadCount < minreads &
                           Order %in% unique(Order[OTURank %in% c("Dominant", "Tied")])
                       ) %>%
                       filter(!(OTURank == "Secondary" & should_remove)) %>%
                       select(-should_remove))

# select the target barcode sequence (if possible)
if(marker == "TRUE"){
  output$OTUDestination <- ifelse(output$Phylum %in% known.non.targets | output$Kingdom %in% known.non.targets, "NTS",
                                     ifelse(output$Kingdom == "unknown", "NTS",
                                              ifelse(output$OTUCategory != "mtCOI", "NTS",
                                                   ifelse(output$Tax_Match == "MISMATCH", "NTS",
                                                          ifelse(output$Ns / output$SeqLength > 0.01, "NTS",
                                                                 "target")))))
} else{
  output$OTUDestination <- ifelse(output$Phylum %in% known.non.targets | output$Kingdom %in% known.non.targets, "NTS",
                                ifelse(output$Kingdom == "unknown", "NTS",
                                       ifelse(output$Tax_Match == "MISMATCH", "NTS",
                                              ifelse(output$Ns / output$SeqLength > 0.01, "NTS",
                                                     "target"))))
}

output <- output %>%
  group_by(Sample) %>%
  group_modify(~{
    df <- .x
    candidates <- df %>%
      filter(OTUDestination == "target") %>%
      arrange(desc(LowRankScore), Ns, desc(ReadCount))

    if (nrow(candidates) == 0) return(df)

    top <- candidates[1, ]
    tied <- candidates %>%
      filter(LowRankScore == top$LowRankScore,
             Ns            == top$Ns,
             ReadCount     == top$ReadCount)

    if (nrow(tied) == 1) {
      df$OTUDestination <- if_else(df$OTUName == top$OTUName, "target", "NTS")
    } else {
      df$OTUDestination <- "NTS"
    }

    df
  }) %>%
  ungroup()

# Add detail to Tax_Match results (mismatch_at_X if confirmed mismatch, blank otherwise)
output$Tax_Match_Annotated <- apply(output, 1, function(row) {
  matches <- row[rank_match_cols]
  if (row["Tax_Match"] == "MISMATCH") {
    for (rank in ranks) {
      if (matches[[paste0(rank, "_Match")]] == "MISMATCH") {
        return(paste0("mismatch_at_", rank))
      }
    }
  }
  return("")
})

# output results table
if(marker == "TRUE"){
  output <- output %>%
    select(
      RunID,
      Sample,
      OTU_Name = OTUName,
      Read_Count = ReadCount,
      Barcode_Status = OTUDestination,
      Sequence_Length = SeqLength,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      Species,
      Ambiguous_Base_Count = Ns,
      Indel = INDEL,
      Stop_Codon = STOPCODON,
      Tax_Match = Tax_Match_Annotated,
      Sequence_Name = SeqName,
      Sequence,
      OTU_Rank = OTURank,
      OTU_Category = OTUCategory,
      Probability_Kingdom = SINTAX_Kingdom_Prob,
      Probability_Phylum = SINTAX_Phylum_Prob,
      Probability_Class = SINTAX_Class_Prob,
      Probability_Order = SINTAX_Order_Prob,
      Probability_Family = SINTAX_Family_Prob,
      Probability_Genus = SINTAX_Genus_Prob,
      Probability_Species = SINTAX_Species_Prob,
      User_Kindgom = USER.Kingdom,
      User_Phylum = USER.Phylum,
      User_Class = USER.Class,
      User_Order = USER.Order,
      User_Family = USER.Family,
      User_Genus = USER.Genus,
      User_Species = USER.Species,
      Kingdom_Match,
      Phylum_Match,
      Class_Match,
      Order_Match,
      Family_Match,
      Genus_Match,
      Species_Match)
} else{
  output <- output %>%
    select(
      RunID,
      Sample,
      OTU_Name = OTUName,
      Read_Count = ReadCount,
      Barcode_Status = OTUDestination,
      Sequence_Length = SeqLength,
      Kingdom,
      Phylum,
      Class,
      Order,
      Family,
      Genus,
      Species,
      Ambiguous_Base_Count = Ns,
      Tax_Match = Tax_Match_Annotated,
      Sequence_Name = SeqName,
      Sequence,
      OTU_Rank = OTURank,
      Probability_Kingdom = SINTAX_Kingdom_Prob,
      Probability_Phylum = SINTAX_Phylum_Prob,
      Probability_Class = SINTAX_Class_Prob,
      Probability_Order = SINTAX_Order_Prob,
      Probability_Family = SINTAX_Family_Prob,
      Probability_Genus = SINTAX_Genus_Prob,
      Probability_Species = SINTAX_Species_Prob,
      User_Kindgom = USER.Kingdom,
      User_Phylum = USER.Phylum,
      User_Class = USER.Class,
      User_Order = USER.Order,
      User_Family = USER.Family,
      User_Genus = USER.Genus,
      User_Species = USER.Species,
      Kingdom_Match,
      Phylum_Match,
      Class_Match,
      Order_Match,
      Family_Match,
      Genus_Match,
      Species_Match)
}
output <- output %>%
  group_by(Sample) %>%
  mutate(
    Sample = if_else(
      Barcode_Status == "NTS",
      paste0(Sample, ".NTS", cumsum(Barcode_Status == "NTS")),
      Sample
    )
  ) %>%
  ungroup() %>%
  mutate(parent_sample = sub("\\.NTS[0-9]+$", "", Sample)) %>%
  arrange(parent_sample,
          ifelse(Barcode_Status == "target", 0, 1),
          desc(Read_Count)) %>%
  select(-parent_sample)

write.table(output, sprintf("%s__%s__OTUDetails.tsv", runid, primers), quote = F, row.names = F, sep = "\t")

# output FASTA files for target and non-target sequences
fasta.target <- DNAStringSet(output$Sequence[output$Barcode_Status == "target"])
names(fasta.target) <- output$Sample[output$Barcode_Status == "target"]

fasta.nts <- DNAStringSet(output$Sequence[output$Barcode_Status == "NTS"])
names(fasta.nts) <- output$Sample[output$Barcode_Status == "NTS"]

writeXStringSet(fasta.target, "sequences_barcodes.fasta")
writeXStringSet(fasta.nts, "sequences_nts.fasta")






