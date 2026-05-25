#!/usr/bin/env Rscript
# PAF to SAM Converter - Genomics Format Conversion Suite
# Convert PAF (Pairwise mApping Format) to SAM (Sequence Alignment/Map) format
# Author: Genomics Conversion Suite
# Version: 1.0.0
# GitHub: https://github.com/ank-man/genodot

suppressPackageStartupMessages(library(optparse))

script_version <- "1.0.0"
script_name <- "PAF to SAM Converter"

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

#' Convert PAF strand to SAM flag
strand_to_flag <- function(strand) {
  if (strand == "+") return(0)    # forward strand
  else if (strand == "-") return(16)  # reverse strand
  else return(0)  # default
}

#' Generate SAM CIGAR string from PAF alignment
generate_cigar <- function(lenAln, numResidueMatches) {
  # Simplified CIGAR generation - in practice, this would need more detailed parsing
  # For now, we'll use a simple approach assuming all matches
  paste0(lenAln, "M")
}

# ─────────────────────────────────────────────
#  CLI Options
# ─────────────────────────────────────────────

option_list <- list(
  make_option(c("-o", "--output"), type = "character",
              help = "Output SAM filename [input.sam]",
              dest = "output_filename"),
  make_option(c("-r", "--reference-header"), type = "character", default = NULL,
              help = "Reference FASTA file for header generation",
              dest = "reference_fasta"),
  make_option(c("-m", "--min-alignment-length"), type = "character", default = "1kb",
              help = "Min alignment length, accepts kb/mb/gb [%default]",
              dest = "min_align"),
  make_option(c("-i", "--min-identity"), type = "numeric", default = 0.9,
              help = "Min percent identity (0-1) [%default]",
              dest = "min_identity"),
  make_option(c("-q", "--min-mapq"), type = "numeric", default = 10,
              help = "Min mapping quality [%default]",
              dest = "min_mapq"),
  make_option(c("-H", "--header-only"), action = "store_true", default = FALSE,
              help = "Generate SAM header only [%default]",
              dest = "header_only"),
  make_option(c("-v", "--version"), action = "store_true", default = FALSE,
              help = "Show version and exit")
)

parser <- OptionParser(
  usage = "%prog [options] input.paf\n\nPAF to SAM Converter - Transform PAF alignments to SAM format\nFor more information, see https://github.com/ank-man/genodot",
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
  opt$output_filename <- sub("\\.paf$", ".sam", input_file)
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
if (!opt$header_only) {
  cat("Applying filters...\n")
  alignments <- alignments[alignments$lenAln >= min_align, ]
  alignments <- alignments[alignments$percentID >= opt$min_identity, ]
  alignments <- alignments[alignments$mapQ >= opt$min_mapq, ]
  
  cat(sprintf("Filtered alignments: %d\n", nrow(alignments)))
  
  if (nrow(alignments) == 0) {
    cat("Error: no alignments remain after filtering. Adjust thresholds and retry.\n")
    quit(status = 1)
  }
}

# ─────────────────────────────────────────────
#  Generate SAM Header
# ─────────────────────────────────────────────

cat("Generating SAM header...\n")

# Get unique references and queries
unique_refs <- unique(alignments$refID)
unique_queries <- unique(alignments$queryID)

# Generate SAM header lines
header_lines <- c(
  "@HD",  # Header line
  "VN:1.6",  # SAM format version
  "SO:unsorted",  # Sort order
  "@PG",  # Program line
  "ID:paf2sam",  # Program ID
  "PN:paf2sam.R",  # Program name
  "VN:1.0.0",  # Program version
  "CL:paf2sam.R"  # Command line
)

# Add reference sequences to header
for (ref in unique_refs) {
  ref_len <- unique(alignments$refLen[alignments$refID == ref])[1]
  header_lines <- c(header_lines, "@SQ", paste0("SN:", ref), paste0("LN:", ref_len))
}

# Add read group information (optional)
for (query in unique_queries) {
  query_len <- unique(alignments$queryLen[alignments$queryID == query])[1]
  header_lines <- c(header_lines, "@RG", paste0("ID:", query), paste0("SM:", query))
}

# ─────────────────────────────────────────────
#  Generate SAM Records
# ─────────────────────────────────────────────

if (!opt$header_only) {
  cat("Generating SAM records...\n")
  
  # Convert PAF to SAM format
  sam_records <- data.frame(
    QNAME = alignments$queryID,  # Query template name
    FLAG = sapply(alignments$strand, strand_to_flag),  # Bitwise flag
    RNAME = alignments$refID,  # Reference sequence name
    POS = alignments$refStart + 1,  # 1-based leftmost mapping position
    MAPQ = alignments$mapQ,  # Mapping quality
    CIGAR = sapply(1:nrow(alignments), function(i) generate_cigar(alignments$lenAln[i], alignments$numResidueMatches[i])),  # CIGAR string
    RNEXT = "*",  # Reference name of the mate/next read
    PNEXT = 0,  # Position of the mate/next read
    TLEN = 0,  # Template length
    SEQ = "*",  # Query sequence
    QUAL = "*",  # Query quality
    NM = alignments$lenAln - alignments$numResidueMatches,  # Edit distance
    AS = alignments$numResidueMatches,  # Alignment score
    stringsAsFactors = FALSE
  )
  
  # Add optional tags
  sam_records$PI <- paste0("PI:i:", alignments$percentID * 100)
} else {
  sam_records <- NULL
}

# ─────────────────────────────────────────────
#  Write SAM Output
# ─────────────────────────────────────────────

cat("Writing SAM output...\n")

# Open output file
sink(opt$output_filename)

# Write header
cat("@HD\tVN:1.6\tSO:unsorted\n")
cat("@PG\tID:paf2sam\tPN:paf2sam.R\tVN:1.0.0\tCL:paf2sam.R\n")

# Write reference sequences
for (ref in unique_refs) {
  ref_len <- unique(alignments$refLen[alignments$refID == ref])[1]
  cat(paste("@SQ", paste0("SN:", ref), paste0("LN:", ref_len), sep = "\t"), "\n")
}

# Write SAM records if not header only
if (!opt$header_only) {
  for (i in 1:nrow(sam_records)) {
    record <- paste(
      sam_records$QNAME[i],
      sam_records$FLAG[i],
      sam_records$RNAME[i],
      sam_records$POS[i],
      sam_records$MAPQ[i],
      sam_records$CIGAR[i],
      sam_records$RNEXT[i],
      sam_records$PNEXT[i],
      sam_records$TLEN[i],
      sam_records$SEQ[i],
      sam_records$QUAL[i],
      paste0("NM:i:", sam_records$NM[i]),
      paste0("AS:i:", sam_records$AS[i]),
      sam_records$PI[i],
      sep = "\t"
    )
    cat(record, "\n")
  }
}

sink()

# ─────────────────────────────────────────────
#  Summary Statistics
# ─────────────────────────────────────────────

if (!opt$header_only) {
  total_alignments <- nrow(sam_records)
  total_bp <- sum(alignments$lenAln)
  avg_identity <- mean(alignments$percentID) * 100
  unique_refs <- length(unique(sam_records$RNAME))
  unique_queries <- length(unique(sam_records$QNAME))
  
  cat(sprintf("\n✓ %s completed successfully!\n", script_name))
  cat(sprintf("Output saved: %s\n", opt$output_filename))
  cat(sprintf("Total alignments: %d\n", total_alignments))
  cat(sprintf("Total aligned bp: %.1f Mb\n", total_bp / 1e6))
  cat(sprintf("Average identity: %.1f%%\n", avg_identity))
  cat(sprintf("Reference sequences: %d\n", unique_refs))
  cat(sprintf("Query sequences: %d\n", unique_queries))
} else {
  cat(sprintf("\n✓ SAM header generated successfully!\n"))
  cat(sprintf("Output saved: %s\n", opt$output_filename))
  cat(sprintf("Reference sequences: %d\n", length(unique_refs)))
}

cat("GitHub: https://github.com/ank-man/genodot\n")
