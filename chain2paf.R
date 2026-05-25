#!/usr/bin/env Rscript
# Chain to PAF Converter - Genomics Format Conversion Suite
# Convert UCSC chain format to PAF (Pairwise mApping Format) for genome alignments
# Author: Genomics Conversion Suite
# Version: 1.0.0
# GitHub: https://github.com/ank-man/genodot

suppressPackageStartupMessages(library(optparse))

script_version <- "1.0.0"
script_name <- "Chain to PAF Converter"

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

#' Parse UCSC chain format
parse_chain_file <- function(filename) {
  con <- file(filename, "r")
  chains <- list()
  current_chain <- NULL
  
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    
    if (grepl("^chain", line)) {
      # Save previous chain if exists
      if (!is.null(current_chain)) {
        chains <- append(chains, list(current_chain))
      }
      
      # Parse chain header
      fields <- strsplit(line, "\\s+")[[1]]
      current_chain <- list(
        score = as.numeric(fields[2]),
        tName = fields[3],
        tSize = as.numeric(fields[4]),
        tStrand = fields[5],
        tStart = as.numeric(fields[6]),
        tEnd = as.numeric(fields[7]),
        qName = fields[8],
        qSize = as.numeric(fields[9]),
        qStrand = fields[10],
        qStart = as.numeric(fields[11]),
        qEnd = as.numeric(fields[12]),
        chainID = fields[13],
        alignments = list()
      )
    } else if (grepl("^\\d+", line) && !is.null(current_chain)) {
      # Parse alignment line
      fields <- as.numeric(strsplit(line, "\\s+")[[1]])
      if (length(fields) >= 3) {
        current_chain$alignments <- append(current_chain$alignments, list(
          size = fields[1],
          dt = fields[2],
          dq = fields[3]
        ))
      }
    }
  }
  
  # Save last chain
  if (!is.null(current_chain)) {
    chains <- append(chains, list(current_chain))
  }
  
  close(con)
  return(chains)
}

#' Convert chain alignment to PAF coordinates
chain_to_paf <- function(chain, min_align = 1000) {
  paf_records <- data.frame(
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
  
  if (length(chain$alignments) == 0) return(paf_records)
  
  # Convert strand format
  strand <- if (chain$qStrand == "+") "+" else "-"
  
  # Process each alignment block
  q_pos <- chain$qStart
  t_pos <- chain$tStart
  
  for (align in chain$alignments) {
    if (align$size < min_align) {
      # Skip small alignments
      q_pos <- q_pos + align$size + align$dq
      t_pos <- t_pos + align$size + align$dt
      next
    }
    
    # Calculate coordinates
    if (strand == "+") {
      q_start <- q_pos
      q_end <- q_pos + align$size
    } else {
      # Reverse strand - coordinates are from end
      q_start <- chain$qSize - q_pos - align$size
      q_end <- chain$qSize - q_pos
    }
    
    t_start <- t_pos
    t_end <- t_pos + align$size
    
    # Add to PAF records
    paf_records <- rbind(paf_records, data.frame(
      queryID = chain$qName,
      queryLen = chain$qSize,
      queryStart = q_start,
      queryEnd = q_end,
      strand = strand,
      refID = chain$tName,
      refLen = chain$tSize,
      refStart = t_start,
      refEnd = t_end,
      numResidueMatches = align$size,  # Assume perfect matches
      lenAln = align$size,
      mapQ = 60,  # Default high quality
      stringsAsFactors = FALSE
    ))
    
    # Update positions
    q_pos <- q_pos + align$size + align$dq
    t_pos <- t_pos + align$size + align$dt
  }
  
  return(paf_records)
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
  make_option(c("-s", "--min-score"), type = "numeric", default = 1000,
              help = "Minimum chain score [%default]",
              dest = "min_score"),
  make_option(c("-v", "--version"), action = "store_true", default = FALSE,
              help = "Show version and exit")
)

parser <- OptionParser(
  usage = "%prog [options] input.chain\n\nChain to PAF Converter - Transform UCSC chain format to PAF alignment format\nFor more information, see https://github.com/ank-man/genodot",
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
  opt$output_filename <- sub("\\.chain$", ".paf", input_file)
}

# parse bp threshold
min_align <- tryCatch(parse_bp(opt$min_align), error = function(e) { 
  cat(e$message, "\n"); quit(status=1) 
})

# ─────────────────────────────────────────────
#  Read and Parse Chain File
# ─────────────────────────────────────────────

cat("Reading chain file...\n")
chains <- parse_chain_file(input_file)
cat(sprintf("Found %d chains\n", length(chains)))

# Filter chains by score
chains <- chains[sapply(chains, function(x) x$score >= opt$min_score)]
cat(sprintf("After score filtering: %d chains\n", length(chains)))

if (length(chains) == 0) {
  cat("Error: no chains meet the minimum score criteria.\n")
  quit(status = 1)
}

# ─────────────────────────────────────────────
#  Convert Chains to PAF
# ─────────────────────────────────────────────

cat("Converting chains to PAF format...\n")

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

for (i in 1:length(chains)) {
  chain_paf <- chain_to_paf(chains[[i]], min_align)
  paf_data <- rbind(paf_data, chain_paf)
}

cat(sprintf("Generated %d PAF alignments\n", nrow(paf_data)))

if (nrow(paf_data) == 0) {
  cat("Error: no alignments generated. Adjust filters and retry.\n")
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
