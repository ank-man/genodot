#!/usr/bin/env Rscript
# BED to PAF Converter - Genomics Format Conversion Suite
# Convert BED format to PAF (Pairwise mApping Format) for genome alignments
# Author: Genomics Conversion Suite
# Version: 1.0.0
# GitHub: https://github.com/ank-man/genodot

suppressPackageStartupMessages(library(optparse))

script_version <- "1.0.0"
script_name <- "BED to PAF Converter"

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

#' Extract sequence length from FASTA file
get_fasta_lengths <- function(fasta_file) {
  if (is.null(fasta_file) || !file.exists(fasta_file)) {
    return(NULL)
  }
  
  con <- file(fasta_file, "r")
  lengths <- list()
  current_seq <- NULL
  
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    
    if (grepl("^>", line)) {
      # New sequence header
      if (!is.null(current_seq)) {
        lengths[[current_seq]] <- seq_length
      }
      current_seq <- sub("^>", "", line)
      current_seq <- strsplit(current_seq, "\\s+")[[1]][1]  # Take first word
      seq_length <- 0
    } else {
      # Sequence line
      if (!is.null(current_seq)) {
        seq_length <- seq_length + nchar(gsub("\\s+", "", line))
      }
    }
  }
  
  # Add last sequence
  if (!is.null(current_seq)) {
    lengths[[current_seq]] <- seq_length
  }
  
  close(con)
  return(lengths)
}

# ─────────────────────────────────────────────
#  CLI Options
# ─────────────────────────────────────────────

option_list <- list(
  make_option(c("-o", "--output"), type = "character",
              help = "Output PAF filename [input.paf]",
              dest = "output_filename"),
  make_option(c("-r", "--reference-fasta"), type = "character", default = NULL,
              help = "Reference FASTA file for sequence lengths",
              dest = "reference_fasta"),
  make_option(c("-q", "--query-fasta"), type = "character", default = NULL,
              help = "Query FASTA file for sequence lengths",
              dest = "query_fasta"),
  make_option(c("-m", "--min-alignment-length"), type = "character", default = "1kb",
              help = "Min alignment length, accepts kb/mb/gb [%default]",
              dest = "min_align"),
  make_option(c("-s", "--default-seq-length"), type = "character", default = "1mb",
              help = "Default sequence length when FASTA not provided [%default]",
              dest = "default_seq_len"),
  make_option(c("-i", "--default-identity"), type = "numeric", default = 0.95,
              help = "Default percent identity (0-1) [%default]",
              dest = "default_identity"),
  make_option(c("-q", "--default-mapq"), type = "numeric", default = 60,
              help = "Default mapping quality [%default]",
              dest = "default_mapq"),
  make_option(c("-t", "--strand"), type = "character", default = "+",
              help = "Default strand for all alignments [+/-] [%default]",
              dest = "default_strand"),
  make_option(c("-v", "--version"), action = "store_true", default = FALSE,
              help = "Show version and exit")
)

parser <- OptionParser(
  usage = "%prog [options] input.bed\n\nBED to PAF Converter - Transform BED regions to PAF alignment format\nFor more information, see https://github.com/ank-man/genodot",
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
  opt$output_filename <- sub("\\.bed$", ".paf", input_file)
}

if (!opt$default_strand %in% c("+", "-")) {
  cat("Error: --strand must be '+' or '-'\n")
  quit(status = 1)
}

# parse bp thresholds
min_align <- tryCatch(parse_bp(opt$min_align), error = function(e) { 
  cat(e$message, "\n"); quit(status=1) 
})
default_seq_len <- tryCatch(parse_bp(opt$default_seq_length), error = function(e) { 
  cat(e$message, "\n"); quit(status=1) 
})

# ─────────────────────────────────────────────
#  Read BED File
# ─────────────────────────────────────────────

cat("Reading BED file...\n")
bed_data <- read.table(input_file, stringsAsFactors = FALSE, row.names = NULL,
                       fill = TRUE, header = FALSE)

if (ncol(bed_data) < 3) {
  cat("Error: BED file must have at least 3 columns (chrom, start, end)\n")
  quit(status = 1)
}

# Assign column names
colnames(bed_data)[1:3] <- c("chrom", "start", "end")
if (ncol(bed_data) >= 4) colnames(bed_data)[4] <- "name"
if (ncol(bed_data) >= 5) colnames(bed_data)[5] <- "score"
if (ncol(bed_data) >= 6) colnames(bed_data)[6] <- "strand"

cat(sprintf("Read %d BED regions\n", nrow(bed_data)))

# ─────────────────────────────────────────────
#  Get Sequence Lengths
# ─────────────────────────────────────────────

cat("Getting sequence lengths...\n")

ref_lengths <- get_fasta_lengths(opt$reference_fasta)
query_lengths <- get_fasta_lengths(opt$query_fasta)

if (!is.null(ref_lengths)) {
  cat(sprintf("Found %d reference sequences in FASTA\n", length(ref_lengths)))
} else {
  cat("Reference FASTA not provided, using default sequence length\n")
}

if (!is.null(query_lengths)) {
  cat(sprintf("Found %d query sequences in FASTA\n", length(query_lengths)))
} else {
  cat("Query FASTA not provided, using default sequence length\n")
}

# ─────────────────────────────────────────────
#  Convert BED to PAF
# ─────────────────────────────────────────────

cat("Converting BED to PAF format...\n")

# Initialize PAF data frame
paf_data <- data.frame(
  queryID = character(),
  queryLen = integer(),
  queryStart = integer(),
  queryEnd = integer(),
  strand = character(),
  refID = character(),
  refLen = integer(),
  refStart = integer(),
  refEnd = integer(),
  numResidueMatches = integer(),
  lenAln = integer(),
  mapQ = integer(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(bed_data)) {
  bed_entry <- bed_data[i, ]
  
  # Calculate alignment length
  lenAln <- bed_entry$end - bed_entry$start
  
  # Apply length filter
  if (lenAln < min_align) next
  
  # Determine strand
  if (ncol(bed_data) >= 6 && !is.na(bed_entry$strand)) {
    strand <- bed_entry$strand
  } else {
    strand <- opt$default_strand
  }
  
  # Determine query name
  if (ncol(bed_data) >= 4 && !is.na(bed_entry$name)) {
    queryID <- bed_entry$name
  } else {
    queryID <- paste0("region_", i)
  }
  
  # Get reference length
  refLen <- if (!is.null(ref_lengths) && bed_entry$chrom %in% names(ref_lengths)) {
    ref_lengths[[bed_entry$chrom]]
  } else {
    default_seq_len
  }
  
  # Get query length
  queryLen <- if (!is.null(query_lengths) && queryID %in% names(query_lengths)) {
    query_lengths[[queryID]]
  } else {
    lenAln  # Use alignment length as estimate
  }
  
  # Calculate coordinates
  refStart <- bed_entry$start  # BED is 0-based
  refEnd <- bed_entry$end
  
  # For query coordinates, assume full query is aligned
  queryStart <- 0
  queryEnd <- queryLen
  
  # Calculate matches based on identity
  numResidueMatches <- round(lenAln * opt$default_identity)
  
  # Add to PAF data
  paf_data <- rbind(paf_data, data.frame(
    queryID = queryID,
    queryLen = queryLen,
    queryStart = queryStart,
    queryEnd = queryEnd,
    strand = strand,
    refID = bed_entry$chrom,
    refLen = refLen,
    refStart = refStart,
    refEnd = refEnd,
    numResidueMatches = numResidueMatches,
    lenAln = lenAln,
    mapQ = opt$default_mapq,
    stringsAsFactors = FALSE
  ))
}

cat(sprintf("Converted %d BED regions to PAF format\n", nrow(paf_data)))

if (nrow(paf_data) == 0) {
  cat("Error: no BED regions converted to PAF format. Adjust filters and retry.\n")
  quit(status = 1)
}

# ─────────────────────────────────────────────
#  Write PAF Output
# ─────────────────────────────────────────────

cat("Writing PAF output...\n")

# Sort PAF by reference and position
paf_data <- paf_data[order(paf_data$refID, paf_data$refStart), ]

# Write PAF file
write.table(paf_data, file = opt$output_filename, sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ─────────────────────────────────────────────
#  Summary Statistics
# ─────────────────────────────────────────────

total_alignments <- nrow(paf_data)
total_bp <- sum(paf_data$lenAln)
avg_identity <- mean(paf_data$numResidueMatches / paf_data$lenAln) * 100
unique_refs <- length(unique(paf_data$refID))
unique_queries <- length(unique(paf_data$queryID))

cat(sprintf("\n✓ %s completed successfully!\n", script_name))
cat(sprintf("Output saved: %s\n", opt$output_filename))
cat(sprintf("Total alignments: %d\n", total_alignments))
cat(sprintf("Total aligned bp: %.1f Mb\n", total_bp / 1e6))
cat(sprintf("Average identity: %.1f%%\n", avg_identity))
cat(sprintf("Reference sequences: %d\n", unique_refs))
cat(sprintf("Query sequences: %d\n", unique_queries))
cat("GitHub: https://github.com/ank-man/genodot\n")
