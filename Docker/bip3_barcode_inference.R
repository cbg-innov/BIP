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
min_rank_match <- args[6]
if (is.na(min_rank_match) || !nzchar(min_rank_match)) {
  stop("min_rank_match (arg 6) is required: one of k, p, c, o, f, g, s")
}
min_rank_match <- tolower(min_rank_match)
if (!min_rank_match %in% c("k", "p", "c", "o", "f", "g", "s")) {
  stop(sprintf("invalid min_rank_match '%s': must be one of k, p, c, o, f, g, s", min_rank_match))
}
# In-frame COI-5P lengths, comma-separated; defined in the advanced parameters of BIP.sh.
coi_expected_lengths <- args[7]
if (is.na(coi_expected_lengths) || !nzchar(coi_expected_lengths)) {
  stop("coi_expected_lengths (arg 7) is required: comma-separated integers, e.g. 640,643,646")
}
# --minreads from BIP.sh: OTUs left with fewer reads than this after negative-control
# subtraction are removed.
minreads_otu <- suppressWarnings(as.numeric(args[8]))
if (is.na(minreads_otu)) {
  stop("minreads (arg 8) is required: the --minreads value passed to BIP.sh")
}
expected_lengths <- suppressWarnings(
  as.numeric(trimws(strsplit(coi_expected_lengths, ",")[[1]]))
)
if (length(expected_lengths) == 0 || anyNA(expected_lengths)) {
  stop(sprintf("coi_expected_lengths must be a comma-separated list of integers (got '%s')",
               coi_expected_lengths))
}

#set working directory based on plate's folder, load required libraries
setwd(wkdir)
library("Biostrings")
library("readxl")
library("dplyr")
library("msa")
library(stringr)
library(DECIPHER)
library(readr)

# input OTU consensus sequences and tax ID results
input.fasta <- readDNAStringSet(sprintf("%s.corrected.fasta", runid))
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
# run under two markers appears more than once. Restrict to this primer pair's rows and
# keep one row per sample before merging.
input.tax <- input.tax[
  paste(input.tax[["Forward Primer Name"]],
        input.tax[["Reverse Primer Name"]], sep = "_") == primers, ]
# A sample can occupy wells on more than one plate; keep every plate it is on (for
# plate-specific control subtraction) before reducing to one row per sample.
sample_plates <- lapply(split(as.character(input.tax$Plate), input.tax$Sample),
                        function(p) sort(unique(p[!is.na(p) & p != ""])))
input.tax <- input.tax[!duplicated(input.tax$Sample), ]

# Target amplicon length drives the sequence-length tie-break during barcode inference.
# It comes from the dictionary via parameters.tsv and is required.
amplen_col <- "Target amplicon length"
if (!amplen_col %in% names(input.tax)) {
  stop(sprintf("required column '%s' not found in parameters.tsv", amplen_col))
}
target_lengths <- unique(suppressWarnings(as.numeric(input.tax[[amplen_col]])))
target_lengths <- target_lengths[!is.na(target_lengths)]
if (length(target_lengths) == 0) {
  stop(sprintf("required column '%s' in parameters.tsv has no numeric value for primer pair %s",
               amplen_col, primers))
}
if (length(target_lengths) > 1) {
  warning(sprintf("primer pair %s has %d different '%s' values in parameters.tsv; using %s",
                  primers, length(target_lengths), amplen_col, target_lengths[1]))
}
target_length <- target_lengths[1]

# ---------------------------------------------------------------------------------------
# Negative-control subtraction (plate-specific), before barcode inference
# ---------------------------------------------------------------------------------------
# A well is a control if its Negative Control value is yes, y or 1 (any case). Every OTU
# from a control well is a contaminant: it is written to the ControlContaminants file and
# removed from everything downstream (OTU table, target and NTS FASTA files). Every other
# OTU on the SAME plate that matches a contaminant (vsearch, identity >= ctrl_match_id over
# the whole alignment) loses the read count of its highest-read matching contaminant. OTUs
# left with fewer than --minreads reads are removed. Read counts from here on are
# post-subtraction.
# Matching uses --iddef 1, which counts terminal gaps, so a partial fragment in a control
# cannot match full-length OTUs. The 99% threshold is deliberately stricter than OTU
# clustering, so conserved markers such as 18S, where related species can share >97%
# identity, are not over-subtracted.
ctrl_match_id <- 0.99
is_control_value <- function(x) {
  v <- tolower(trimws(as.character(x)))
  !is.na(x) & (v %in% c("yes", "y") | grepl("^1(\\.0+)?$", v))
}
ctrl_samples <- input.tax$Sample[is_control_value(input.tax[["Negative Control"]])]
plates_of <- function(s) if (s %in% names(sample_plates)) sample_plates[[s]] else character(0)
output$Plate <- vapply(output$Sample, function(s) {
  p <- plates_of(s)
  if (length(p) == 0) NA_character_ else paste(p, collapse = ";")
}, character(1), USE.NAMES = FALSE)
output$Reads_Subtracted <- rep(0, nrow(output))
is_ctrl_otu  <- output$Sample %in% ctrl_samples
contaminants <- output[is_ctrl_otu, ]
# rep(), not a bare 0: a run with no control-well OTUs gives a zero-row table
contaminants$OTUs_Matched           <- rep(0, nrow(contaminants))
contaminants$Reads_Subtracted_Total <- rep(0, nrow(contaminants))
contaminants$OTUs_Removed           <- rep(0, nrow(contaminants))
subtraction_log <- data.frame(RunID = character(0), Plate = character(0), Sample = character(0),
                              OTU_Name = character(0), Reads_Before = numeric(0),
                              Reads_Subtracted = numeric(0), Reads_After = numeric(0),
                              Removed = logical(0), Matched_Control = character(0),
                              Matched_Control_OTU = character(0), Pct_Identity = numeric(0),
                              Sequence_Name = character(0), stringsAsFactors = FALSE)
removed_by_sub <- rep(FALSE, nrow(output))

# plates that at least one control OTU came from
ctrl_plates <- unique(unlist(lapply(unique(contaminants$Sample), plates_of)))
cand <- which(!is_ctrl_otu &
              vapply(output$Sample, function(s) any(plates_of(s) %in% ctrl_plates), logical(1),
                     USE.NAMES = FALSE))
if (nrow(contaminants) > 0 && length(cand) > 0) {
  q <- DNAStringSet(output$Sequence[cand]);  names(q) <- paste0("q", seq_along(cand))
  d <- DNAStringSet(contaminants$Sequence);  names(d) <- paste0("c", seq_len(nrow(contaminants)))
  qf <- tempfile(fileext = ".fasta"); dbf <- tempfile(fileext = ".fasta"); hf <- tempfile(fileext = ".tsv")
  writeXStringSet(q, qf); writeXStringSet(d, dbf)
  # Exhaustive search (maxaccepts/maxrejects 0) so every matching contaminant is found.
  status <- system2("vsearch", c("--usearch_global", qf, "--db", dbf, "--id", ctrl_match_id,
                                 "--iddef", "1", "--maxaccepts", "0", "--maxrejects", "0",
                                 "--strand", "plus", "--userout", hf,
                                 "--userfields", "query+target+id", "--threads", "4", "--quiet"))
  if (!identical(status, 0L)) stop("vsearch failed during negative-control matching")
  hits <- if (file.exists(hf) && file.size(hf) > 0) {
    read.delim(hf, header = FALSE, col.names = c("q", "c", "id"), stringsAsFactors = FALSE)
  } else data.frame(q = character(0), c = character(0), id = numeric(0))
  unlink(c(qf, dbf, hf))
  if (nrow(hits) > 0) {
    hits$oi <- cand[as.integer(sub("^q", "", hits$q))]          # row in output
    hits$ci <- as.integer(sub("^c", "", hits$c))                # row in contaminants
    same_plate <- mapply(function(s, cs) any(plates_of(s) %in% plates_of(cs)),
                         output$Sample[hits$oi], contaminants$Sample[hits$ci], USE.NAMES = FALSE)
    hits <- hits[as.logical(same_plate), ]                                  # same plate only
  }
  if (nrow(hits) > 0) {
    # each OTU is charged its highest-read matching contaminant on its plate
    hits$creads <- contaminants$ReadCount[hits$ci]
    hits <- hits[order(hits$oi, -hits$creads, -hits$id), ]
    best <- hits[!duplicated(hits$oi), ]
    before <- output$ReadCount[best$oi]
    after  <- before - best$creads
    gone   <- after < minreads_otu
    subtraction_log <- data.frame(
      RunID = runid, Plate = output$Plate[best$oi], Sample = output$Sample[best$oi],
      OTU_Name = output$OTUName[best$oi], Reads_Before = before, Reads_Subtracted = best$creads,
      Reads_After = pmax(after, 0), Removed = gone, Matched_Control = contaminants$Sample[best$ci],
      Matched_Control_OTU = contaminants$OTUName[best$ci], Pct_Identity = best$id,
      Sequence_Name = output$SeqName[best$oi], stringsAsFactors = FALSE)
    output$ReadCount[best$oi]        <- pmax(after, 0)
    output$Reads_Subtracted[best$oi] <- best$creads
    removed_by_sub[best$oi]          <- gone
    for (j in unique(best$ci)) {
      k <- best$ci == j
      contaminants$OTUs_Matched[j]           <- sum(k)
      contaminants$Reads_Subtracted_Total[j] <- sum(pmin(best$creads[k], before[k]))
      contaminants$OTUs_Removed[j]           <- sum(gone[k])
    }
    # keep the read count in each sequence name in step with the new count
    hit_rows <- best$oi[!gone]
    output$SeqName[hit_rows] <- mapply(function(n, r) sub("\\|reads-[0-9]+", paste0("|reads-", r), n),
                                       output$SeqName[hit_rows], as.integer(output$ReadCount[hit_rows]),
                                       USE.NAMES = FALSE)
  }
}
contam_out <- data.frame(
  RunID = rep(runid, nrow(contaminants)), Plate = contaminants$Plate,
  Control_Sample = contaminants$Sample, OTU_Name = contaminants$OTUName,
  Read_Count = contaminants$ReadCount, Sequence_Length = contaminants$SeqLength,
  Kingdom = contaminants$Kingdom, Phylum = contaminants$Phylum, Class = contaminants$Class,
  Order = contaminants$Order, Family = contaminants$Family, Genus = contaminants$Genus,
  Species = contaminants$Species, OTUs_Matched_On_Plate = contaminants$OTUs_Matched,
  Reads_Subtracted_Total = contaminants$Reads_Subtracted_Total,
  OTUs_Removed = contaminants$OTUs_Removed, Sequence_Name = contaminants$SeqName,
  Sequence = contaminants$Sequence, stringsAsFactors = FALSE)
write.table(contam_out, sprintf("%s__%s__ControlContaminants.tsv", runid, primers),
            quote = FALSE, row.names = FALSE, sep = "\t")
write.table(subtraction_log, sprintf("%s__%s__ControlSubtraction.tsv", runid, primers),
            quote = FALSE, row.names = FALSE, sep = "\t")
cat(sprintf("Negative controls (%s): %d control OTU(s) on %d plate(s); %d OTU(s) had reads subtracted, %d removed (< %d reads)\n",
            primers, nrow(contaminants), length(unique(contaminants$Plate)),
            nrow(subtraction_log), sum(subtraction_log$Removed), as.integer(minreads_otu)))
output <- output[!is_ctrl_otu & !removed_by_sub, ]

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

# Non-target status is derived from the taxonomy in the parameters file provided by user

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

# Tax_Match is enforced from Kingdom down to --min_rank_match: a mismatch at or above
# that rank blocks target designation. Ranks below it cannot block it, and instead
# contribute to LowRankScore, the first tie-breaker between candidate OTUs.
min_rank <- unname(rank_prefix_map[min_rank_match])
if (is.na(min_rank)) {
  stop(sprintf("invalid min_rank_match '%s': must be one of %s",
               min_rank_match, paste(names(rank_prefix_map), collapse = ", ")))
}
min_rank_i <- match(min_rank, ranks)

high_cols <- rank_match_cols[seq_len(min_rank_i)]
low_cols  <- if (min_rank_i < length(ranks)) {
  rank_match_cols[(min_rank_i + 1):length(ranks)]
} else {
  character(0)
}

output$Tax_Match <- ifelse(
  apply(output[, high_cols, drop = FALSE] == "MISMATCH", 1, any), "MISMATCH",
  ifelse(
    apply(output[, high_cols, drop = FALSE] == "UNKNOWN", 1, any), "UNKNOWN",
    "OK"
  )
)

# No ranks below the enforced one means nothing is left to score, so every OTU ties here.
output$LowRankScore <- if (length(low_cols) == 0) {
  rep(0L, nrow(output))
} else {
  Reduce(`+`, lapply(output[, low_cols, drop = FALSE],
                     function(x) as.integer(!is.na(x) & x == "OK")))
}

# add dominant tag to output table
output <- output[order(output$Sample, -output$ReadCount),]
# rep(), not a bare value: control subtraction can leave a zero-row table
output$OTURank <- rep("Secondary", nrow(output))
topotu <- which(!duplicated(output$Sample))
output$OTURank[topotu] <- "Dominant"

# add tied tag to output table: when other OTUs have the same read count as the Dominant,
# all of them, the Dominant included, become Tied, however many there are (a pairwise loop
# missed the third OTU in a three-way tie and crashed on a single-row table)
top_reads <- ave(output$ReadCount, output$Sample, FUN = max)
is_top    <- output$ReadCount == top_reads
n_top     <- ave(as.integer(is_top), output$Sample, FUN = sum)
output$OTURank[is_top & n_top > 1] <- "Tied"

# The in-frame COI checks below (indel flag, NUMT category, mtCOI required for a target)
# assume a full-length barcode, so they apply only when this primer pair's target amplicon
# length lies inside the coi_expected_lengths range - the same condition BIP.sh uses to
# decide whether bip2 runs. A shorter COI amplicon (e.g. 418 bp) is judged like any other
# marker; it still gets the COI output columns.
min_expected_length <- min(expected_lengths)
max_expected_length <- max(expected_lengths)
full_length_coi <- marker == "TRUE" &&
  target_length >= min_expected_length && target_length <= max_expected_length

# add INDEL tag to table (default to FALSE if not full-length COI)
# expected_lengths comes from coi_expected_lengths in the advanced parameters of BIP.sh
ifelse(full_length_coi,
       output$INDEL <- ifelse(output$SeqLength %in% expected_lengths, FALSE, TRUE),
       output$INDEL <- rep(FALSE, nrow(output)))
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
    # SINTAX appends the probability as "(0.51)"; split on that final bracket only,
    # because library names can contain brackets, e.g. "Homo_sapiens_(human)".
    m      <- regmatches(inner, regexec("^(.*)\\(([0-9.eE+-]+)\\)$", inner))[[1]]
    value  <- if (length(m) == 3) m[2] else inner
    prob   <- if (length(m) == 3) m[3] else NA_character_
    result[paste0("SINTAX_", rank)]          <- value
    result[paste0("SINTAX_", rank, "_Prob")] <- prob
  }
  result
}

df_split <- if (nrow(output) > 0) {
  as.data.frame(t(sapply(output$TaxAssign, parse_sintax)), stringsAsFactors = FALSE)
} else {
  # control subtraction removed every OTU: keep the columns, with no rows
  as.data.frame(t(parse_sintax(NA))[0, , drop = FALSE], stringsAsFactors = FALSE)
}

# Combine the split data back into the original data frame
output <- cbind(output, df_split)

# add OTU category to table
# A sequence outside the range spanned by coi_expected_lengths (set in the advanced
# parameters of BIP.sh) is too short or too long to be a real barcode, so it is a NUMT.
output$OTUCategory <- if(full_length_coi){
                            ifelse(output$SeqLength < min_expected_length, "NUMT",
                              ifelse(output$SeqLength > max_expected_length, "NUMT",
                                ifelse(output$INDEL == FALSE & output$STOPCODON == TRUE, "NUMT",
                                    ifelse(output$HMM != "HMM_ISSUE", "mtCOI", "unknown"))))
} else{
  output$OTUCategory <- rep("unknown", nrow(output))
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
if(full_length_coi){
  output$OTUDestination <- ifelse(output$Kingdom == "unknown", "NTS",
                                       ifelse(output$OTUCategory != "mtCOI", "NTS",
                                            ifelse(output$Tax_Match == "MISMATCH", "NTS",
                                                   ifelse(output$Ns / output$SeqLength > 0.01, "NTS",
                                                          "target"))))
} else{
  output$OTUDestination <- ifelse(output$Kingdom == "unknown", "NTS",
                                       ifelse(output$Tax_Match == "MISMATCH", "NTS",
                                              ifelse(output$Ns / output$SeqLength > 0.01, "NTS",
                                                     "target")))
}

output <- output %>%
  group_by(Sample) %>%
  group_modify(~{
    df <- .x
    # Tie-breaking, in order: matches at ranks below --min_rank_match, then distance
    # from the target amplicon length (so 655 beats 664 against a 658 target), then
    # ambiguous bases, then read count. OTUs equidistant from the target length stay
    # tied on length and are separated by the criteria that follow.
    candidates <- df %>%
      filter(OTUDestination == "target") %>%
      mutate(LengthDist = abs(SeqLength - target_length)) %>%
      arrange(desc(LowRankScore), LengthDist, Ns, desc(ReadCount))

    if (nrow(candidates) == 0) return(df)

    top <- candidates[1, ]
    tied <- candidates %>%
      filter(LowRankScore == top$LowRankScore,
             LengthDist    == top$LengthDist,
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
# (apply() calls the function once with a dummy row on a zero-row table, so skip it then)
output$Tax_Match_Annotated <- if (nrow(output) == 0) character(0) else apply(output, 1, function(row) {
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
      User_Kingdom = USER.Kingdom,
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
      Species_Match,
      Reads_Subtracted)
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
      User_Kingdom = USER.Kingdom,
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
      Species_Match,
      Reads_Subtracted)
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






