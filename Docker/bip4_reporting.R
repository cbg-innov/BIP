# Barcode Inference Pipeline Report Generator
# Usage: Rscript generate_report.R <runid> <wkdir>

args <- commandArgs(trailingOnly = TRUE)
runid <- args[1]
wkdir <- args[2]

setwd(wkdir)

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(gtable)
library(scales)

# ============================================================
# 1. IMPORT DATA
# ============================================================

# --- Parameters file ---
params <- read.table(file.path(wkdir, "parameters.tsv"),
                     header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)

# parameters.tsv has one row per WELL and covers EVERY primer pair, so a sample run
# under two markers appears twice - on two different plates. Sample alone is therefore
# NOT a unique key here: joining on it fans out and cross-contaminates the markers
# (a COI OTU would be attributed to the sample's 28S plate as well). The real key is
# Sample + primer_pair, so carry the primer pair alongside every sample.
params$primer_pair <- paste(params[["Forward Primer Name"]],
                            params[["Reverse Primer Name"]], sep = "_")

# One row per (sample, primer pair) for joining against the OTU tables.
sample_plate_map <- params %>%
  distinct(Sample, primer_pair, Plate)

# Control wells: any row with a non-empty "Negative Control" value
ctrl_samples <- params$Sample[
  !is.na(params[["Negative Control"]]) &
  nchar(trimws(as.character(params[["Negative Control"]]))) > 0
]

# --- Read counts file ---
readcounts <- read.table(file.path(wkdir, paste0(runid, "_readcounts.tsv")),
                         header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
colnames(readcounts) <- c("runid", "stage", "primer_pair", "reads")

# --- OTU tables: find all primer pair directories and import ---
primer_dirs <- list.dirs(wkdir, full.names = TRUE, recursive = FALSE)
primer_dirs <- primer_dirs[file.exists(file.path(primer_dirs,
                                                 list.files(primer_dirs[1], pattern = "_OTUDetails.tsv", full.names = FALSE)[1]))]

# More robust: find any directory containing an OTUDetails file
all_dirs <- list.dirs(wkdir, full.names = TRUE, recursive = FALSE)
otu_files <- lapply(all_dirs, function(d) {
  f <- list.files(d, pattern = "_OTUDetails\\.tsv$", full.names = TRUE)
  if (length(f) > 0) f else NULL
})
otu_files <- unlist(Filter(Negate(is.null), otu_files))

otu_all <- bind_rows(lapply(otu_files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, comment.char = "")
  # Extract primer pair from filename: {runid}__{primerpair}__OTUDetails.tsv
  bn <- basename(f)
  pp <- gsub(paste0("^", runid, "__(.+)__OTUDetails\\.tsv$"), "\\1", bn)
  df$primer_pair <- pp
  df
}))

# ============================================================
# 2. STACKED READ RETENTION HISTOGRAM
# ============================================================

stage_order <- c("Raw", "PrimerSplit", "LengthFilter", "LinkedUMI",
                 "ReadsInOTUs", "ReadsInTargetBarcodes")

stage_labels <- c(
  "Raw"                   = "Raw Reads",
  "PrimerSplit"           = "With Primers",
  "LengthFilter"          = "After Length Filter",
  "LinkedUMI"             = "After Demultiplexing",
  "ReadsInOTUs"           = "Reads in OTUs",
  "ReadsInTargetBarcodes" = "Reads in target barcodes"
)

rc <- readcounts %>% filter(stage %in% stage_order)
rc$stage <- factor(rc$stage, levels = stage_order)

# Raw reads are logged once, as primer_pair = "ALL": before primer splitting a read
# is not attributable to any marker (and here 14% of them go on to match neither).
# So the Raw bar is drawn as ONE neutral bar rather than being split across the
# primer pairs. It used to be divided equally (reads / n_pp), which produced a
# 50/50 stack that looked like a measurement but was pure fabrication.
primer_pairs <- unique(rc$primer_pair[rc$primer_pair != "ALL"])
n_pp <- length(primer_pairs)

rc_plot <- rc

# Total reads per stage for label placement
stage_totals <- rc_plot %>%
  group_by(stage) %>%
  summarise(total = sum(reads), .groups = "drop")

raw_total <- rc$reads[rc$stage == "Raw" & rc$primer_pair == "ALL"][1]

stage_totals <- stage_totals %>%
  mutate(pct = paste0("(", round(total / raw_total * 100, 0), "% remaining)"),
         label = paste0(format(total, big.mark = ",", trim = TRUE), "\n", pct))

# Colour palette for primer pairs, plus a neutral grey for the unassigned Raw bar
custom_palette <- c("#2166AC", "#D6604D", "#4DAF4A", "#FF7F00", "#984EA3", "#A65628")
pp_colours <- setNames(
  if (n_pp <= length(custom_palette)) custom_palette[seq_len(n_pp)] else hue_pal()(n_pp),
  primer_pairs
)
# Raw/unassigned bar: a dark green carried at low opacity, so it reads as present
# but recessive next to the saturated per-primer bars. Alpha is baked into the hex
# (RRGGBBAA) rather than set via alpha= so it applies only to this fill.
pp_colours["ALL"] <- "#2E7D3259"

# Order the legend so the primer pairs come first and "ALL" (raw, unassigned) last
rc_plot$primer_pair <- factor(rc_plot$primer_pair, levels = c(primer_pairs, "ALL"))

plot1 <- ggplot(rc_plot, aes(x = factor(stage, levels = stage_order),
                             y = reads, fill = primer_pair)) +
  geom_bar(stat = "identity") +
  geom_text(data = stage_totals,
            aes(x = factor(stage, levels = stage_order),
                y = total,
                label = label),
            inherit.aes = FALSE,
            vjust = -0.2, size = 2.8, fontface = "bold") +
  scale_x_discrete(labels = stage_labels) +
  # breaks = primer_pairs: the raw/unassigned bar keeps its fill but is left out of
  # the legend, which only needs to decode the per-primer colours
  scale_fill_manual(values = pp_colours, breaks = primer_pairs, name = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, max(stage_totals$total) * 1.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 9, angle = 30, hjust = 1,
                                   color = "black", face = "bold"),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# ============================================================
# 3. READS/SAMPLE AND OTUs/SAMPLE HISTOGRAMS (faceted by primer pair)
# ============================================================

# Aggregate reads and OTUs per sample per primer pair
# Strip .NTS[0-9]+ suffix so secondary OTUs are counted under their parent sample
sample_stats <- otu_all %>%
  mutate(Sample = sub("\\.NTS[0-9]+$", "", Sample)) %>%
  group_by(Sample, primer_pair) %>%
  summarise(Reads = sum(Read_Count),
            OTUs  = n_distinct(OTU_Name),
            .groups = "drop")

# Bin breaks (same as old script)
read_breaks <- c(-1, 0, 5, 10, 15, 20, 25, 50, 100, 200, 300, 400, 500, Inf)
read_labels <- c("0","1-5","6-10","11-15","16-20","21-25","26-50",
                 "51-100","101-200","201-300","301-400","401-500",">500")

otu_breaks <- c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, Inf)
otu_labels  <- c("0","1","2","3","4","5","6","7","8","9","10",">10")

sample_stats$bin_reads <- cut(sample_stats$Reads, breaks = read_breaks, labels = read_labels)
sample_stats$bin_otus  <- cut(sample_stats$OTUs,  breaks = otu_breaks,  labels = otu_labels)

reads_hist <- sample_stats %>%
  group_by(primer_pair, bin_reads) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  complete(primer_pair, bin_reads, fill = list(Freq = 0))

otus_hist <- sample_stats %>%
  group_by(primer_pair, bin_otus) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  complete(primer_pair, bin_otus, fill = list(Freq = 0))

plot2 <- ggplot(reads_hist, aes(x = factor(bin_reads, levels = read_labels), y = Freq)) +
  geom_bar(stat = "identity", fill = "blue", colour = "black") +
  facet_wrap(~ primer_pair, scales = "free_y", ncol = 1) +
  theme_bw() +
  ylab("No. Samples") + xlab("Reads/Sample") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        strip.text = element_text(face = "bold", size = 8))

plot3 <- ggplot(otus_hist, aes(x = factor(bin_otus, levels = otu_labels), y = Freq)) +
  geom_bar(stat = "identity", fill = "forestgreen", colour = "black") +
  facet_wrap(~ primer_pair, scales = "free_y", ncol = 1) +
  theme_bw() +
  ylab("No. Samples") + xlab("OTUs/Sample") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, face = "bold"),
        strip.text = element_text(face = "bold", size = 8))

# ============================================================
# 4. SUCCESS BY PLATE
# ============================================================

# Denominator: non-control samples per plate PER primer pair, taken straight from
# the parameters file (not the OTU results). Sourcing this from params is what makes
# EVERY plate assigned to a primer set appear in that primer's facet — including
# plates that produced no target barcodes, or no reads at all. Those simply show a
# 0% bar. Counting per (Plate, primer_pair) also stops a plate shared across two
# markers from lending its full sample count to both facets.
plate_primer_totals <- params %>%
  filter(!Sample %in% ctrl_samples) %>%
  group_by(Plate, primer_pair) %>%
  summarise(Total = n(), .groups = "drop")

# Numerator: unique non-control samples with a target barcode, per plate per primer pair
target_per_plate <- otu_all %>%
  filter(Barcode_Status == "target") %>%
  mutate(Sample = sub("__.*$", "", Sample)) %>%
  filter(!Sample %in% ctrl_samples) %>%
  left_join(sample_plate_map, by = c("Sample", "primer_pair")) %>%
  filter(!is.na(Plate)) %>%
  distinct(Plate, primer_pair, Sample) %>%
  group_by(Plate, primer_pair) %>%
  summarise(Targets = n(), .groups = "drop")

plate_success <- plate_primer_totals %>%
  left_join(target_per_plate, by = c("Plate", "primer_pair")) %>%
  mutate(
    Targets = ifelse(is.na(Targets), 0L, Targets),
    Pct     = Targets / Total * 100,
    Colour  = ifelse(Pct >= 75, "#66BB6A",
              ifelse(Pct >= 50, "#FFB300", "firebrick"))
  ) %>%
  arrange(primer_pair, Plate)

n_plates <- length(unique(plate_success$Plate))
plate_text_size <- max(4, min(8, 120 / n_plates))

plate_levels <- rev(sort(unique(plate_success$Plate)))

plot_success <- ggplot(plate_success,
                       aes(x = Pct,
                           y = factor(Plate, levels = plate_levels),
                           fill = Colour)) +
  geom_bar(stat = "identity") +
  # scales = "free": free_y drops each panel to only the plates belonging to that
  # primer set (so a 28S-only plate does not occupy an empty row in the COI facet);
  # x stays 0-100% via the fixed limits below.
  facet_wrap(~ primer_pair, ncol = n_pp, scales = "free") +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%")) +
  scale_y_discrete(labels = function(x) gsub("-", "\uad", x)) +
  scale_fill_identity() +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 9, face = "bold", color = "black"),
    axis.text.y  = element_text(size = plate_text_size, face = "bold", color = "black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor    = element_blank(),
    axis.ticks          = element_blank(),
    strip.text          = element_text(face = "bold", size = 8)
  ) +
  labs(x = "Barcode Success Rate")

# ============================================================
# 5. PLATE HEATMAPS
# ============================================================

# Build plate layout from parameters. primer_pair is carried through so reads are
# matched per marker: a sample run under two markers occupies a well on each of two
# plates, and joining on Sample alone would put its combined read total in both.
plate_layout <- params %>%
  select(Plate, Well, Sample, primer_pair) %>%
  mutate(WellRow = substr(Well, 1, 1),
         WellCol = sprintf("%02d", as.integer(substr(Well, 2, nchar(Well)))))

# Infer plate grid dimensions from actual well labels
well_rows <- LETTERS[1:which(LETTERS == max(plate_layout$WellRow))]
well_cols <- sprintf("%02d", 1:max(as.integer(plate_layout$WellCol)))

# Aggregate OTU reads per sample, PER MARKER (summing across markers would show a
# shared sample's combined COI+28S total on both of its wells).
# The .NTS[0-9]+ suffix is stripped as well, so secondary OTUs are counted under
# their parent sample - without it those rows never match plate_layout and their
# reads vanish from the heatmap (cf. sample_stats, which already does this).
reads_per_sample_all <- otu_all %>%
  mutate(Sample = sub("\\.NTS[0-9]+$", "", sub("__.*$", "", Sample))) %>%
  group_by(Sample, primer_pair) %>%
  summarise(TotalReads = sum(Read_Count), .groups = "drop")

# The well heatmaps are rendered by print_heatmap_pages() further down.

# ============================================================
# 5. SUMMARY TABLE
# ============================================================

df_nonctrl <- otu_all %>%
  mutate(Sample = sub("__.*$", "", Sample)) %>%
  filter(!Sample %in% ctrl_samples)

df_ctrl <- otu_all %>%
  mutate(Sample = sub("__.*$", "", Sample)) %>%
  filter(Sample %in% ctrl_samples)

total_reads        <- sum(df_nonctrl$Read_Count)
total_samples_seq  <- length(unique(sub("\\.NTS[0-9]+$", "", df_nonctrl$Sample[df_nonctrl$Read_Count > 0])))
total_otus         <- nrow(df_nonctrl)
avg_reads_sample   <- round(total_reads / max(total_samples_seq, 1), 0)
avg_otus_sample    <- round(total_otus  / max(total_samples_seq, 1), 1)
avg_reads_otu      <- round(total_reads / max(total_otus, 1), 0)

n_target           <- nrow(df_nonctrl[df_nonctrl$Barcode_Status == "target", ])
n_nontarget        <- nrow(df_nonctrl[df_nonctrl$Barcode_Status != "target", ])
pct_nontarget      <- round(n_nontarget / max(total_otus, 1) * 100, 1)

total_nonctrl_samples <- length(unique(params$Sample[!params$Sample %in% ctrl_samples]))

# Per-primer denominator: sum params-based non-control counts across plates for each primer
primer_totals <- plate_primer_totals %>%
  group_by(primer_pair) %>%
  summarise(TotalSamples = sum(Total), .groups = "drop")

primer_success_df <- otu_all %>%
  filter(Barcode_Status == "target") %>%
  mutate(Sample = sub("__.*$", "", Sample)) %>%
  filter(!Sample %in% ctrl_samples) %>%
  distinct(Sample, primer_pair) %>%
  group_by(primer_pair) %>%
  summarise(WithTarget = n_distinct(Sample), .groups = "drop") %>%
  left_join(primer_totals, by = "primer_pair") %>%
  mutate(
    Metric = paste0("Barcode success rate (", primer_pair, ")"),
    # pmax(), not max(): TotalSamples is a column here, one row per primer pair.
    # max() collapses it to a single scalar - the largest denominator across ALL
    # primer pairs - so every primer's success rate got divided by the biggest
    # primer's sample count.
    Value  = paste0(round(WithTarget / pmax(TotalSamples, 1) * 100, 0), "%")
  ) %>%
  select(Metric, Value) %>%
  as.data.frame(stringsAsFactors = FALSE)

if (nrow(df_ctrl) > 0) {
  n_control_otus        <- nrow(df_ctrl)
  plates_with_ctrl      <- length(unique(df_ctrl$Sample[df_ctrl$Read_Count > 0]))
  n_plates              <- length(unique(params$Plate))
  avg_ctrl_otus_plate   <- round(n_control_otus / max(n_plates, 1), 1)
  ctrl_otus_gt10        <- nrow(df_ctrl[df_ctrl$Read_Count > 10, ])
  max_ctrl_reads        <- max(df_ctrl$Read_Count)
} else {
  n_control_otus      <- 0
  plates_with_ctrl    <- 0
  avg_ctrl_otus_plate <- 0
  ctrl_otus_gt10      <- 0
  max_ctrl_reads      <- 0
}

summary_df <- data.frame(
  Metric = c(
    "Total reads",
    "Total samples with sequences",
    "Average reads per sample",
    "Average OTUs per sample",
    "Average reads per OTU",
    "Number of control OTUs",
    "Average number of control OTUs per plate",
    "Number of control OTUs with > 10 reads",
    "Max reads per control OTU",
    "Total OTUs",
    "Non-target OTUs",
    "Percent non-target OTUs"
  ),
  Value = c(
    format(total_reads,       big.mark = ","),
    format(total_samples_seq, big.mark = ","),
    format(avg_reads_sample,  big.mark = ","),
    avg_otus_sample,
    format(avg_reads_otu,     big.mark = ","),
    n_control_otus,
    avg_ctrl_otus_plate,
    ctrl_otus_gt10,
    format(max_ctrl_reads,    big.mark = ","),
    format(total_otus,        big.mark = ","),
    n_nontarget,
    paste0(pct_nontarget, "%")
  ),
  stringsAsFactors = FALSE
)

# Category colours:
# Rows 1-5:  Run Performance  -> grey85 / grey92 alternating
# Rows 6-9:  Neg Control      -> moccasin shades
# Rows 10-12: OTU Summary     -> palegreen3 shades
# Row 13:    Barcode success  -> steelblue1

perf_cols  <- c("grey85",   "grey92",   "grey85",   "grey92",   "grey85")
ctrl_cols  <- c("moccasin", "#FFE4A0",  "moccasin", "#FFE4A0")
otu_cols   <- c("palegreen3","#7DC87D", "palegreen3")
bcs_col    <- "steelblue1"

row_fills <- c(perf_cols, ctrl_cols, otu_cols, bcs_col)

# Section header rows get a darker shade
header_rows <- c(1, 6, 10)
row_fills[header_rows] <- c("grey60", "#CC9900", "forestgreen")

# Font colours: headers white, others black
text_cols <- rep("black", 13)
text_cols[header_rows] <- "white"

# Fontface: headers bold, barcode success rate bold, others plain
text_face <- rep("plain", 13)
text_face[c(header_rows, 13)] <- "bold"

# Insert section header rows into summary_df
n_bcs <- nrow(primer_success_df)

summary_df_display <- rbind(
  data.frame(Metric = "Run Performance Summary",  Value = "", stringsAsFactors = FALSE),
  summary_df[1:5, ],
  data.frame(Metric = "Negative Control Summary", Value = "", stringsAsFactors = FALSE),
  summary_df[6:9, ],
  data.frame(Metric = "OTU Summary",              Value = "", stringsAsFactors = FALSE),
  summary_df[10:12, ],
  primer_success_df
)

# Row positions in display df:
# 1        = Run Performance header   -> grey60,       white,  bold
# 2-6      = perf data rows           -> grey85/92 alt, black, plain
# 7        = Neg Control header       -> #CC9900,       white,  bold
# 8-11     = ctrl data rows           -> moccasin alt,  black, plain
# 12       = OTU Summary header       -> forestgreen,   white,  bold
# 13-15    = otu data rows            -> palegreen alt,  black, plain
# 16+      = Barcode success rows     -> steelblue1,    black,  bold (one per primer pair)

n_disp <- nrow(summary_df_display)

disp_fills <- c(
  "grey60",                                          # Run perf header
  "grey85", "grey92", "grey85", "grey92", "grey85",  # perf rows 1-5
  "#CC9900",                                         # ctrl header
  "moccasin", "#FFE4A0", "moccasin", "#FFE4A0",      # ctrl rows 6-9
  "forestgreen",                                     # OTU header
  "palegreen3", "#7DC87D", "palegreen3",             # OTU rows 10-12
  rep("steelblue1", n_bcs)                           # barcode success rows
)

disp_textcol <- c(
  "white",                          # Run perf header
  rep("black", 5),                  # perf rows
  "white",                          # ctrl header
  rep("black", 4),                  # ctrl rows
  "white",                          # OTU header
  rep("black", 3),                  # OTU rows
  rep("black", n_bcs)               # barcode success rows
)

disp_face <- c(
  "bold",                           # Run perf header
  rep("plain", 5),                  # perf rows
  "bold",                           # ctrl header
  rep("plain", 4),                  # ctrl rows
  "bold",                           # OTU header
  rep("plain", 3),                  # OTU rows
  rep("bold", n_bcs)                # barcode success rows
)

table_theme <- ttheme_default(
  core = list(
    bg_params = list(
      fill = cbind(disp_fills, disp_fills),
      col  = NA
    ),
    fg_params = list(
      hjust    = cbind(rep(0,    n_disp), rep(1,    n_disp)),
      x        = cbind(rep(0.02, n_disp), rep(0.98, n_disp)),
      fontface = cbind(disp_face, disp_face),
      col      = cbind(disp_textcol, disp_textcol),
      cex      = rep(0.95, n_disp)
    )
  ),
  colhead = list(
    bg_params = list(fill = NA, col = NA),
    fg_params = list(col = NA, cex = 0)   # hide column headers
  )
)

table.1 <- tableGrob(summary_df_display, rows = NULL, theme = table_theme)
table.1$widths[1] <- unit(4.5, "in")
table.1$widths[2] <- unit(1.5, "in")

# Outer border
table.1 <- gtable_add_grob(table.1,
                           grobs = rectGrob(gp = gpar(lwd = 4, fill = NA)),
                           t = 1, b = nrow(table.1), l = 1, r = ncol(table.1))

# ============================================================
# 6. HELPER FUNCTIONS FOR PDF LAYOUT
# ============================================================

# --- Histogram page printer ---
# Prints a histogram plot (already faceted by primer pair) with dynamic height.
# 1 primer  -> occupies top half of page
# 2 primers -> full page, 1 col x 2 rows
# 3+ primers -> 2 cols x ceil(n/2) rows, full page
print_hist_page <- function(plot_obj, title_label, n_primers) {
  if (n_primers == 1) {
    print(ggdraw() +
            draw_plot_label(title_label, size = 14, x = 0.5, y = 0.97, hjust = 0.5) +
            draw_plot(plot_obj, x = 0.05, y = 0.48, width = 0.90, height = 0.46))
  } else if (n_primers == 2) {
    print(ggdraw() +
            draw_plot_label(title_label, size = 14, x = 0.5, y = 0.97, hjust = 0.5) +
            draw_plot(plot_obj + facet_wrap(~ primer_pair, ncol = 1),
                      x = 0.05, y = 0.03, width = 0.90, height = 0.91))
  } else {
    print(ggdraw() +
            draw_plot_label(title_label, size = 14, x = 0.5, y = 0.97, hjust = 0.5) +
            draw_plot(plot_obj + facet_wrap(~ primer_pair, ncol = 2),
                      x = 0.05, y = 0.03, width = 0.90, height = 0.91))
  }
}

# --- Heatmap page printer ---
# Plates are fixed at 4cm wide x 3cm tall, 6 columns.
# Global colour scale shared across all pages.
# Each page prints up to 6 cols x floor((page_height_cm - margins) / 3) rows.

print_heatmap_pages <- function(reads_df_all, plate_layout, well_rows, well_cols, title_all) {
  
  margin_cm  <- 1.5
  panel_w_cm <- 3.5
  panel_h_cm <- 2.8
  n_cols     <- 4
  page_h_cm  <- 11 * 2.54
  
  # rows that fit per page given fixed panel height (+ ~1cm per row for strip/axis)
  n_rows_page     <- max(1, floor((page_h_cm - margin_cm * 2 - 1.0) / (panel_h_cm + 1.0)))
  plates_per_page <- n_cols * n_rows_page
  
  plates <- unique(plate_layout$Plate)
  
  render_heatmap_pages <- function(reads_df, page_title) {
    
    all_wells <- expand.grid(WellRow = well_rows, WellCol = well_cols,
                             stringsAsFactors = FALSE)
    # Join on Sample AND primer_pair. On Sample alone, a sample run under two markers
    # picks up the other marker's reads: its 28S well would be painted with its COI
    # read count even when the 28S reaction produced nothing.
    df_joined <- plate_layout %>%
      left_join(reads_df, by = c("Sample", "primer_pair")) %>%
      mutate(TotalReads = ifelse(is.na(TotalReads), 0, TotalReads))
    
    df_full <- bind_rows(lapply(plates, function(p) {
      p_df <- df_joined %>% filter(Plate == p)
      merge(all_wells, p_df, by = c("WellRow","WellCol"), all.x = TRUE) %>%
        mutate(Plate = p,
               TotalReads = ifelse(is.na(TotalReads), 0, TotalReads))
    }))
    df_full$WellRow <- factor(df_full$WellRow, levels = rev(well_rows))
    df_full$WellCol <- factor(df_full$WellCol, levels = well_cols)
    
    global_max <- max(df_full$TotalReads, na.rm = TRUE)
    if (is.infinite(global_max) || global_max == 0) global_max <- 1
    
    page_groups <- split(plates, ceiling(seq_along(plates) / plates_per_page))
    
    for (pg in seq_along(page_groups)) {
      pg_plates <- page_groups[[pg]]
      df_page   <- df_full %>% filter(Plate %in% pg_plates)
      df_page$Plate <- factor(df_page$Plate, levels = pg_plates)
      
      p <- ggplot(df_page, aes(x = WellCol, y = WellRow, fill = TotalReads)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        facet_wrap(~ Plate, ncol = n_cols) +
        scale_fill_gradientn(
          colours  = c("white", "#FFFF00", "#FFA500", "#FF2200", "#8B0000"),
          values   = c(0, 0.25, 0.5, 0.75, 1),
          limits   = c(0, global_max),
          na.value = "grey90",
          name     = "Reads", labels = comma) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1),
              axis.text.y = element_text(size = 5),
              axis.title  = element_blank(),
              strip.text  = element_text(face = "bold", size = 7),
              panel.grid  = element_blank(),
              panel.spacing = unit(0.3, "cm"),
              legend.position = "right",
              legend.key.height = unit(1.5, "cm"))
      
      # Convert to gtable and FORCE panel cells to fixed physical size
      g <- ggplotGrob(p)
      panel_cols <- unique(g$layout$l[grepl("^panel", g$layout$name)])
      panel_rows <- unique(g$layout$t[grepl("^panel", g$layout$name)])
      g$widths[panel_cols]  <- unit(panel_w_cm, "cm")
      g$heights[panel_rows] <- unit(panel_h_cm, "cm")
      
      grid.newpage()
      grid.text(if (pg == 1) page_title else paste(page_title, "(cont.)"),
                x = 0.5, y = unit(page_h_cm - 0.8, "cm"),
                gp = gpar(fontsize = 13, fontface = "bold"))
      
      # Draw the gtable centred near the top; it now has its own intrinsic size
      pushViewport(viewport(
        x      = unit(margin_cm, "cm"),
        y      = unit(page_h_cm - margin_cm - 1.4, "cm"),
        width  = grobWidth(g),
        height = grobHeight(g),
        just   = c("left", "top")
      ))
      grid.draw(g)
      popViewport()
    }
  }
  
  render_heatmap_pages(reads_df_all, title_all)
}

# ============================================================
# 7. OUTPUT PDF
# ============================================================

pdf(file.path(wkdir, paste0(runid, "_Report.pdf")), width = 8.5, height = 11)

# --- Page 1: Summary table ---
print(ggdraw() +
        draw_plot_label("Barcode Inference Pipeline Report", size = 22,
                        color = "blue", x = 0.5, y = 0.97, hjust = 0.5) +
        draw_plot_label(runid, size = 16, x = 0.5, y = 0.92, hjust = 0.5) +
        draw_plot_label(format(Sys.time(), "%B %d, %Y"), size = 11,
                        x = 0.5, y = 0.88, hjust = 0.5) +
        draw_plot(table.1, x = 0.05, y = 0.15, width = 0.90, height = 0.68))

# --- Page 2: Read retention stacked bar ---
print(ggdraw() +
        draw_plot_label("Read Retention by Stage", size = 14,
                        x = 0.5, y = 0.97, hjust = 0.5) +
        draw_plot(plot1, x = 0.05, y = 0.45, width = 0.90, height = 0.45))

# --- Reads/sample histogram ---
print_hist_page(plot2, "Reads per Sample", n_pp)

# --- OTUs/sample histogram ---
print_hist_page(plot3, "OTUs per Sample", n_pp)

# --- Success by plate ---
print(ggdraw() +
        draw_plot_label("Success by Plate", size = 14,
                        x = 0.5, y = 0.97, hjust = 0.5) +
        draw_plot(plot_success, x = 0.05, y = 0.03, width = 0.90, height = 0.91))

# --- Heatmaps ---
print_heatmap_pages(reads_per_sample_all, plate_layout, well_rows, well_cols, "Reads per Well")

dev.off()

message(sprintf("Report written to: %s/%s_Report.pdf", wkdir, runid))
