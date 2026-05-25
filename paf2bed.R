#!/usr/bin/env Rscript
# PAF to BED Converter - Genomics Format Conversion Suite
# Convert PAF (Pairwise mApping Format) to BED format
# Author: Genomics Conversion Suite
# Version: 1.0.0
# GitHub: https://github.com/ank-man/genodot

suppressPackageStartupMessages(library(optparse))

script_version <- "1.0.0"
script_name <- "PAF to BED Converter"

# ─────────────────────────────────────────────
#  Helper Functions
# ─────────────────────────────────────────────

#' Parse human-friendly bp strings like 400kb, 1.5mb, 2gb to numeric
parse_bp <- function(x) {
  x <- trimws(tolower(as.character(x)))
  multipliers <- c(kb = 1e3, mb = 1e6, gb = 1e9, tb = 1e12)
  for (suffix in names(multipliers)) {
    if (grepl(paste0(suffix, "$"), x)) {
      num <- as.numeric(sub(suffix, "", x))
      if (!is.na(num)) return(num * multipliers[suffix])
    }
  }
  val <- suppressWarnings(as.numeric(x))
  if (is.na(val)) stop(paste("Cannot parse bp value:", x))
  val
}

# ─────────────────────────────────────────────
#  CLI Options
# ─────────────────────────────────────────────

option_list <- list(
  make_option(c("-o", "--output"), type = "character",
              help = "Output BED filename [input.bed]",
              dest = "output_filename"),
  make_option(c("-m", "--min-alignment-length"), type = "character", default = "1kb",
              help = "Min alignment length, accepts kb/mb/gb [%default]",
              dest = "min_align"),
  make_option(c("-i", "--min-identity"), type = "numeric", default = 0.9,
              help = "Min percent identity (0-1) [%default]",
              dest = "min_identity"),
  make_option(c("-q", "--min-mapq"), type = "numeric", default = 10,
              help = "Min mapping quality [%default]",
              dest = "min_mapq"),
  make_option(c("-r", "--reference-only"), action = "store_true", default = FALSE,
              help = "Output reference coordinates only [%default]",
              dest = "ref_only"),
  make_option(c("-s", "--strand-specific"), action = "store_true", default = FALSE,
              help = "Include strand information in name field [%default]",
              dest = "strand_specific"),
  make_option(c("-v", "--version"), action = "store_true", default = FALSE,
              help = "Show version and exit")
)

parser <- OptionParser(
  usage = "%prog [options] input.paf\n\nPAF to BED Converter - Transform genome alignments to BED format\nFor more information, see https://github.com/ank-man/genodot",
  option_list = option_list
)
opts <- parse_args(parser, positional_arguments = c(0, 1))
opt  <- opts$options

if (opt$version) {
  cat(paste0(script_name, " version: ", script_version, "\n"))
  cat("GitHub: https://github.com/ank-man/genodot\n")
  quit(status = 0)
}

# ─────────────────────────────────────────────
#  Input Validation
# ─────────────────────────────────────────────

input_file <- opts$args
if (length(input_file) == 0) {
  cat("Error: missing input file.\n\n")
  print_help(parser)
  quit(status = 1)
} else if (file.access(input_file, mode = 4) == -1) {
  cat(sprintf("Error: '%s' does not exist or cannot be read.\n\n", input_file))
  print_help(parser)
  quit(status = 1)
}

if (is.null(opt$output_filename)) {
  opt$output_filename <- sub("\\.paf$", ".bed", input_file)
}

# parse bp threshold
min_align <- tryCatch(parse_bp(opt$min_align), error = function(e) { 
  cat(e$message, "\n"); quit(status=1) 
})

# ─────────────────────────────────────────────
#  Read and Process PAF
# ─────────────────────────────────────────────

cat("Reading PAF file...\n")
alignments <- read.table(input_file, stringsAsFactors = FALSE, row.names = NULL,
                         fill = TRUE, header = FALSE)[, 1:12]
alignments[, c(2:4, 7:12)] <- lapply(alignments[, c(2:4, 7:12)], as.numeric)
colnames(alignments)[1:12] <- c("queryID", "queryLen", "queryStart", "queryEnd",
                                "strand", "refID", "refLen", "refStart", "refEnd",
                                "numResidueMatches", "lenAln", "mapQ")

# Calculate percent identity
alignments$percentID <- alignments$numResidueMatches / alignments$lenAln

# Apply filters
cat("Applying filters...\n")
alignments <- alignments[alignments$lenAln >= min_align, ]
alignments <- alignments[alignments$percentID >= opt$min_identity, ]
alignments <- alignments[alignments$mapQ >= opt$min_mapq, ]

cat(sprintf("Filtered alignments: %d\n", nrow(alignments)))

if (nrow(alignments) == 0) {
  cat("Error: no alignments remain after filtering. Adjust thresholds and retry.\n")
  quit(status = 1)
}

# ─────────────────────────────────────────────
#  Generate BED Output
# ─────────────────────────────────────────────

cat("Generating BED output...\n")

if (opt$ref_only) {
  # Reference coordinates only
  bed_data <- data.frame(
    chrom = alignments$refID,
    start = alignments$refStart,
    end = alignments$refEnd,
    name = paste0(alignments$queryID, ":", alignments$queryStart, "-", alignments$queryEnd),
    score = round(alignments$percentID * 1000),
    strand = alignments$strand,
    stringsAsFactors = FALSE
  )
} else {
  # Both reference and query coordinates (standard BED6)
  bed_data <- data.frame(
    chrom = alignments$refID,
    start = alignments$refStart,
    end = alignments$refEnd,
    name = if (opt$strand_specific) {
      paste0(alignments$queryID, "_", alignments$strand, "_", 
             round(alignments$percentID * 100), "%")
    } else {
      paste0(alignments$queryID, "_", round(alignments$percentID * 100), "%")
    },
    score = round(alignments$percentID * 1000),
    strand = alignments$strand,
    stringsAsFactors = FALSE
  )
}

# Sort BED file by chromosome and position
bed_data <- bed_data[order(bed_data$chrom, bed_data$start), ]

# Write BED file
write.table(bed_data, file = opt$output_filename, sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ─────────────────────────────────────────────
#  Summary Statistics
# ─────────────────────────────────────────────

total_alignments <- nrow(bed_data)
total_bp <- sum(bed_data$end - bed_data$start)
avg_identity <- mean(alignments$percentID) * 100
unique_refs <- length(unique(bed_data$chrom))
unique_queries <- length(unique(alignments$queryID))

cat(sprintf("\n✓ %s completed successfully!\n", script_name))
cat(sprintf("Output saved: %s\n", opt$output_filename))
cat(sprintf("Total alignments: %d\n", total_alignments))
cat(sprintf("Total aligned bp: %.1f Mb\n", total_bp / 1e6))
cat(sprintf("Average identity: %.1f%%\n", avg_identity))
cat(sprintf("Reference sequences: %d\n", unique_refs))
cat(sprintf("Query sequences: %d\n", unique_queries))
cat("GitHub: https://github.com/ank-man/genodot\n")
