#!/usr/bin/env Rscript
# SAM to PAF Converter - Genomics Format Conversion Suite
# Convert SAM (Sequence Alignment/Map) format to PAF (Pairwise mApping Format)
# Author: Genomics Conversion Suite
# Version: 1.0.0
# GitHub: https://github.com/ank-man/genodot

suppressPackageStartupMessages(library(optparse))

script_version <- "1.0.0"
script_name <- "SAM to PAF Converter"

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

#' Parse SAM flag to get strand information
flag_to_strand <- function(flag) {
  if (bitwAnd(flag, 16) == 16) return("-")  # reverse strand
  else return("+")  # forward strand
}

#' Parse CIGAR string to get alignment length
parse_cigar_length <- function(cigar) {
  # Extract numbers from CIGAR string and sum for M, =, X operations
  matches <- regmatches(cigar, gregexpr("[0-9]+[M=X]", cigar, perl = TRUE))[[1]]
  if (length(matches) == 0) return(0)
  
  total_length <- 0
  for (match in matches) {
    num <- as.numeric(sub("[M=X]$", "", match))
    op <- substr(match, nchar(match), nchar(match))
    if (op %in% c("M", "=", "X")) {
      total_length <- total_length + num
    }
  }
  return(total_length)
}

#' Parse CIGAR string to get number of matches (simplified)
parse_cigar_matches <- function(cigar) {
  # Simplified approach - assume all M operations are matches
  # In practice, this would need more detailed parsing
  matches <- regmatches(cigar, gregexpr("[0-9]+[M=X]", cigar, perl = TRUE))[[1]]
  if (length(matches) == 0) return(0)
  
  total_matches <- 0
  for (match in matches) {
    num <- as.numeric(sub("[M=X]$", "", match))
    total_matches <- total_matches + num
  }
  return(total_matches)
}

#' Read SAM file and extract header information
read_sam_header <- function(filename) {
  con <- file(filename, "r")
  header_refs <- list()
  
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (substr(line, 1, 1) != "@") {
      # Push back the first alignment line
      seek(con, where = 0, origin = "current")
      break
    }
    
    if (grepl("^@SQ", line)) {
      # Parse @SQ lines for reference sequences
      fields <- strsplit(line, "\t")[[1]]
      sn <- sub("^SN:", "", fields[grepl("^SN:", fields)])
      ln <- as.numeric(sub("^LN:", "", fields[grepl("^LN:", fields)]))
      header_refs[[sn]] <- ln
    }
  }
  
  close(con)
  return(header_refs)
}

# ─────────────────────────────────────────────
#  CLI Options
# ─────────────────────────────────────────────

option_list <- list(
  make_option(c("-o", "--output"), type = "character",
              help = "Output PAF filename [input.paf]",
              dest = "output_filename"),
  make_option(c("-m", "--min-alignment-length"), type = "character", default = "1kb",
              help = "Min alignment length, accepts kb/mb/gb [%default]",
              dest = "min_align"),
  make_option(c("-q", "--min-mapq"), type = "numeric", default = 10,
              help = "Min mapping quality [%default]",
              dest = "min_mapq"),
  make_option(c("-u", "--include-unmapped"), action = "store_true", default = FALSE,
              help = "Include unmapped reads in output [%default]",
              dest = "include_unmapped"),
  make_option(c("-s", "--estimate-query-length"), action = "store_true", default = FALSE,
              help = "Estimate query length from alignment if not available [%default]",
              dest = "estimate_query_len"),
  make_option(c("-v", "--version"), action = "store_true", default = FALSE,
              help = "Show version and exit")
)

parser <- OptionParser(
  usage = "%prog [options] input.sam\n\nSAM to PAF Converter - Transform SAM alignments to PAF format\nFor more information, see https://github.com/ank-man/genodot",
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
  opt$output_filename <- sub("\\.(sam|bam)$", ".paf", input_file)
}

# parse bp threshold
min_align <- tryCatch(parse_bp(opt$min_align), error = function(e) { 
  cat(e$message, "\n"); quit(status=1) 
})

# ─────────────────────────────────────────────
#  Read SAM Header
# ─────────────────────────────────────────────

cat("Reading SAM header...\n")
header_refs <- read_sam_header(input_file)
cat(sprintf("Found %d reference sequences in header\n", length(header_refs)))

# ─────────────────────────────────────────────
#  Read and Process SAM
# ─────────────────────────────────────────────

cat("Reading SAM file...\n")

# Read SAM file (skip header lines starting with @)
sam_data <- read.table(input_file, stringsAsFactors = FALSE, row.names = NULL,
                       fill = TRUE, header = FALSE, comment.char = "@")

if (ncol(sam_data) < 11) {
  cat("Error: SAM file must have at least 11 columns\n")
  quit(status = 1)
}

# Assign column names (standard SAM format)
colnames(sam_data)[1:11] <- c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", 
                              "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL")

# Parse optional tags (if present)
if (ncol(sam_data) > 11) {
  sam_data$TAGS <- apply(sam_data[, 12:ncol(sam_data)], 1, paste, collapse = "\t")
} else {
  sam_data$TAGS <- ""
}

# Filter out unmapped reads if requested
if (!opt$include_unmapped) {
  sam_data <- sam_data[!is.na(sam_data$POS) & sam_data$POS > 0, ]
}

cat(sprintf("Total SAM records: %d\n", nrow(sam_data)))

# ─────────────────────────────────────────────
#  Convert SAM to PAF
# ─────────────────────────────────────────────

cat("Converting to PAF format...\n")

# Initialize data frame for PAF
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

for (i in 1:nrow(sam_data)) {
  record <- sam_data[i, ]
  
  # Skip unmapped reads
  if (is.na(record$POS) || record$POS == 0 || record$RNAME == "*") {
    next
  }
  
  # Parse alignment information
  strand <- flag_to_strand(record$FLAG)
  lenAln <- parse_cigar_length(record$CIGAR)
  numResidueMatches <- parse_cigar_matches(record$CIGAR)
  
  # Apply length filter
  if (lenAln < min_align) next
  
  # Apply MAPQ filter
  if (record$MAPQ < opt$min_mapq) next
  
  # Get reference length from header
  refLen <- if (record$RNAME %in% names(header_refs)) {
    header_refs[[record$RNAME]]
  } else {
    NA
  }
  
  # Get query length
  queryLen <- if (record$SEQ != "*") {
    nchar(record$SEQ)
  } else if (opt$estimate_query_len) {
    lenAln  # Estimate from alignment length
  } else {
    NA
  }
  
  # Calculate coordinates
  refStart <- record$POS - 1  # Convert to 0-based
  refEnd <- refStart + lenAln
  
  # For query coordinates, assume full query is aligned (simplified)
  queryStart <- 0
  queryEnd <- if (!is.na(queryLen)) queryLen else lenAln
  
  # Add to PAF data
  paf_data <- rbind(paf_data, data.frame(
    queryID = record$QNAME,
    queryLen = queryLen,
    queryStart = queryStart,
    queryEnd = queryEnd,
    strand = strand,
    refID = record$RNAME,
    refLen = refLen,
    refStart = refStart,
    refEnd = refEnd,
    numResidueMatches = numResidueMatches,
    lenAln = lenAln,
    mapQ = record$MAPQ,
    stringsAsFactors = FALSE
  ))
}

cat(sprintf("Converted %d alignments to PAF format\n", nrow(paf_data)))

if (nrow(paf_data) == 0) {
  cat("Error: no alignments converted to PAF format. Adjust filters and retry.\n")
  quit(status = 1)
}

# ─────────────────────────────────────────────
#  Write PAF Output
# ─────────────────────────────────────────────

cat("Writing PAF output...\n")

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
