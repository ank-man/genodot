#!/usr/bin/env Rscript
# GenoDot - Advanced Genome Alignment Dot Plot Visualizer
# The ultimate PAF to dotplot converter for genome synteny analysis
# Author: Genomics Visualization Suite
# Version: 3.0.0
# GitHub: https://github.com/ank-man/genodot

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

script_version <- "3.0.0"
script_name <- "GenoDot"

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

#' Format bp axis labels with smart scaling
fmt_bp <- function(x, step) {
  if (step > 1e12)      paste(round(x / 1e12, 1), "Tb")
  else if (step > 1e9)  paste(round(x / 1e9, 1), "Gb")
  else if (step > 1e6)  paste(round(x / 1e6, 1), "Mb")
  else if (step > 1e3)  paste(round(x / 1e3, 1), "kb")
  else                  paste(round(x), "bp")
}

#' Generate n visually distinct colors with multiple palettes
generate_colors <- function(n, palette = "Set3") {
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    stop("Package 'RColorBrewer' is needed. Install with install.packages('RColorBrewer')")
  
  if (n <= 8 && palette == "Set2") {
    return(RColorBrewer::brewer.pal(max(n, 3), "Set2")[seq_len(n)])
  } else if (n <= 12) {
    return(RColorBrewer::brewer.pal(max(n, 3), palette)[seq_len(n)])
  } else {
    return(colorRampPalette(RColorBrewer::brewer.pal(12, palette))(n))
  }
}

#' Compute per-sequence cumulative offsets as a named numeric vector
cumulative_offsets <- function(sizes) {
  offsets <- c(0, cumsum(as.numeric(sizes[-length(sizes)])))
  names(offsets) <- names(sizes)
  offsets
}

#' Add enhanced BED annotation lines with legends
add_bed_lines <- function(gp, bed_file, seq_sizes, axis = c("x", "y"), flip_ids = NULL, 
                         legend_title = "Annotations") {
  axis <- match.arg(axis)
  bed <- tryCatch(
    read.table(bed_file, header = FALSE, stringsAsFactors = FALSE),
    error = function(e) stop(paste("Cannot read BED file:", bed_file, "\n", e$message))
  )
  if (ncol(bed) < 3) stop(paste("BED file must have at least 3 columns:", bed_file))
  colnames(bed)[1:3] <- c("seqID", "start", "end")
  
  # Add name column if available (4th column)
  if (ncol(bed) >= 4) {
    colnames(bed)[4] <- "name"
  } else {
    bed$name <- paste0("Region_", seq_len(nrow(bed)))
  }

  unknown <- setdiff(bed$seqID, names(seq_sizes))
  if (length(unknown) > 0)
    warning(paste("BED seqIDs not found in alignments:", paste(unknown, collapse = ", ")))

  bed <- bed[bed$seqID %in% names(seq_sizes), ]
  if (nrow(bed) == 0) return(gp)

  # handle flipped query sequences
  if (!is.null(flip_ids)) {
    rev_idx <- which(bed$seqID %in% flip_ids)
    if (length(rev_idx) > 0) {
      qmax <- seq_sizes[bed$seqID[rev_idx]]
      bed$start[rev_idx] <- qmax - bed$start[rev_idx] + 1
      bed$end[rev_idx]   <- qmax - bed$end[rev_idx]   + 1
    }
  }

  offsets       <- cumulative_offsets(seq_sizes)
  bed$pos_start <- bed$start + offsets[bed$seqID]
  bed$pos_end   <- bed$end   + offsets[bed$seqID]
  
  # Use unique colors for each unique annotation name
  unique_names <- unique(bed$name)
  colors <- generate_colors(length(unique_names))
  names(colors) <- unique_names

  for (i in seq_len(nrow(bed))) {
    positions <- c(bed$pos_start[i], bed$pos_end[i])
    color <- colors[bed$name[i]]
    if (axis == "x") {
      gp <- gp + geom_vline(xintercept = positions, linetype = "dashed",
                            color = color, linewidth = 0.6, alpha = 0.8)
    } else {
      gp <- gp + geom_hline(yintercept = positions, linetype = "twodash",
                            color = color, linewidth = 0.6, alpha = 0.8)
    }
  }
  gp
}

#' Calculate alignment statistics
calc_alignment_stats <- function(alignments) {
  total_alignments <- nrow(alignments)
  total_bp_aligned <- sum(alignments$lenAln)
  avg_identity <- mean(alignments$percentID) * 100
  median_identity <- median(alignments$percentID) * 100
  
  list(
    total_alignments = total_alignments,
    total_bp_aligned = total_bp_aligned,
    avg_identity = avg_identity,
    median_identity = median_identity
  )
}

# ─────────────────────────────────────────────
#  CLI Options - Enhanced
# ─────────────────────────────────────────────

option_list <- list(
  make_option(c("-o", "--output"), type = "character",
              help = "Output filename prefix [input.paf]",
              dest = "output_filename"),
  make_option(c("-p", "--plot-size"), type = "numeric", default = 15,
              help = "Plot width in inches [%default]",
              dest = "plot_size"),
  make_option(c("-f", "--flip"), action = "store_true", default = FALSE,
              help = "Flip queries where most alignments are reverse-complement [%default]",
              dest = "flip"),
  make_option(c("-b", "--break-point"), action = "store_true", default = FALSE,
              help = "Show alignment break points [%default]",
              dest = "break_point"),
  make_option(c("-c", "--identity-floor"), type = "numeric", default = 0,
              help = "Clamp percent identity below this value to this value [%default]",
              dest = "identity_floor"),
  make_option(c("-a", "--alpha"), type = "numeric", default = 0.7,
              help = "Segment transparency, 0–1 [%default]",
              dest = "alpha"),
  make_option(c("-w", "--line-width"), type = "numeric", default = 0.4,
              help = "Segment line width [%default]",
              dest = "line_width"),
  make_option(c("-s", "--sort-by-refid"), action = "store_true", default = FALSE,
              help = "Sort reference IDs alphabetically (default: by length) [%default]",
              dest = "sortbyID"),
  make_option(c("-q", "--min-query-length"), type = "character", default = "400kb",
              help = "Min total alignment length per query, accepts kb/mb/gb [%default]",
              dest = "min_query_aln"),
  make_option(c("-m", "--min-alignment-length"), type = "character", default = "10kb",
              help = "Min individual alignment length, accepts kb/mb/gb [%default]",
              dest = "min_align"),
  make_option(c("-r", "--min-ref-len"), type = "character", default = "1mb",
              help = "Min reference sequence length, accepts kb/mb/gb [%default]",
              dest = "min_ref_len"),
  make_option(c("-i", "--reference-ids"), type = "character", default = NULL,
              help = "Comma-separated reference IDs to keep and order [%default]",
              dest = "refIDs"),
  make_option(c("-F", "--output-format"), type = "character", default = "both",
              help = "Output format: pdf, png, svg, or both [%default]",
              dest = "output_format"),
  make_option(c("-e", "--ref-bed"), type = "character", default = NULL,
              help = "Reference BED file for vertical marker lines",
              dest = "ref_bed_file"),
  make_option(c("-E", "--query-bed"), type = "character", default = NULL,
              help = "Query BED file for horizontal marker lines",
              dest = "query_bed_file"),
  make_option(c("-t", "--title"), type = "character", default = NULL,
              help = "Custom plot title [%default]",
              dest = "plot_title"),
  make_option(c("-C", "--color-palette"), type = "character", default = "RdYlBu",
              help = "Color palette for identity: RdYlBu, Viridis, Plasma, Heat [%default]",
              dest = "color_palette"),
  make_option(c("-S", "--show-stats"), action = "store_true", default = FALSE,
              help = "Show detailed alignment statistics [%default]",
              dest = "show_stats"),
  make_option(c("-v", "--version"), action = "store_true", default = FALSE,
              help = "Show version and exit")
)

options(error = traceback)
parser <- OptionParser(
  usage = "%prog [options] input.paf\n\nGenoDot - Advanced PAF Dot Plot Visualizer\nFor more information, see https://github.com/ank-man/genodot",
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

if (is.null(opt$output_filename)) opt$output_filename <- input_file

opt$output_format <- tolower(trimws(opt$output_format))
valid_formats <- c("pdf", "png", "svg", "both")
if (!opt$output_format %in% valid_formats) {
  cat(paste("Error: --output-format must be one of:", paste(valid_formats, collapse = ", "), "\n"))
  quit(status = 1)
}

# parse bp thresholds
min_query_aln <- tryCatch(parse_bp(opt$min_query_aln), error = function(e) { cat(e$message, "\n"); quit(status=1) })
min_align     <- tryCatch(parse_bp(opt$min_align),     error = function(e) { cat(e$message, "\n"); quit(status=1) })
min_ref_len   <- tryCatch(parse_bp(opt$min_ref_len),   error = function(e) { cat(e$message, "\n"); quit(status=1) })

# ─────────────────────────────────────────────
#  Read & Prepare Alignments
# ─────────────────────────────────────────────

cat("Reading PAF file...\n")
alignments <- read.table(input_file, stringsAsFactors = FALSE, row.names = NULL,
                         fill = TRUE, header = FALSE)[, 1:12]
alignments[, c(2:4, 7:12)] <- lapply(alignments[, c(2:4, 7:12)], as.numeric)
colnames(alignments)[1:12] <- c("queryID", "queryLen", "queryStart", "queryEnd",
                                "strand", "refID", "refLen", "refStart", "refEnd",
                                "numResidueMatches", "lenAln", "mapQ")

# Convert PAF zero-based half-open coords to 1-based inclusive
alignments$queryStart <- alignments$queryStart + 1
alignments$refStart   <- alignments$refStart   + 1

# Percent identity clamped to [identity_floor, 1]
alignments$percentID <- alignments$numResidueMatches / alignments$lenAln
if (opt$identity_floor > 0)
  alignments$percentID <- pmax(alignments$percentID, opt$identity_floor)

# For reverse-strand alignments, swap query start/end so segments draw correctly
rev_idx <- which(alignments$strand == "-")
tmp <- alignments$queryStart[rev_idx]
alignments$queryStart[rev_idx] <- alignments$queryEnd[rev_idx]
alignments$queryEnd[rev_idx]   <- tmp
rm(tmp, rev_idx)

cat(sprintf("Alignments read:    %d\nQuery sequences:    %d\n",
            nrow(alignments), length(unique(alignments$queryID))))

# ─────────────────────────────────────────────
#  Reference Ordering & Filtering
# ─────────────────────────────────────────────

if (is.null(opt$refIDs)) {
  chromMax <- tapply(alignments$refLen, alignments$refID, max)
  refIDsToKeepOrdered <- if (opt$sortbyID) sort(names(chromMax))
                         else names(sort(chromMax, decreasing = TRUE))
} else {
  refIDsToKeepOrdered <- trimws(unlist(strsplit(opt$refIDs, ",")))
  alignments <- alignments[alignments$refID %in% refIDsToKeepOrdered, ]
}

# Single-pass filtering
queryLenAgg <- tapply(alignments$lenAln, alignments$queryID, sum)
alignments  <- alignments[alignments$queryID %in% names(queryLenAgg)[queryLenAgg > min_query_aln], ]
alignments  <- alignments[alignments$lenAln > min_align, ]
alignments  <- alignments[alignments$refLen > min_ref_len, ]

# Re-check query aggregate after other filters
queryLenAgg <- tapply(alignments$lenAln, alignments$queryID, sum)
alignments  <- alignments[alignments$queryID %in% names(queryLenAgg)[queryLenAgg > min_query_aln], ]

cat(sprintf("After filtering:    %d alignments | %d queries\n\n",
            nrow(alignments), length(unique(alignments$queryID))))

if (nrow(alignments) == 0) {
  cat("Error: no alignments remain after filtering. Adjust thresholds and retry.\n")
  quit(status = 1)
}

refIDsToKeepOrdered <- refIDsToKeepOrdered[refIDsToKeepOrdered %in% alignments$refID]

# ─────────────────────────────────────────────
#  Compute Plot Coordinates
# ─────────────────────────────────────────────

alignments$refID <- factor(alignments$refID, levels = refIDsToKeepOrdered)
alignments <- alignments[order(alignments$refID, alignments$refStart), ]

chromMax  <- tapply(alignments$refLen, alignments$refID, max)
refOffset <- cumulative_offsets(chromMax)

alignments$refStart2 <- alignments$refStart + refOffset[as.character(alignments$refID)]
alignments$refEnd2   <- alignments$refEnd   + refOffset[as.character(alignments$refID)]

# Sort queries: first by position of their longest alignment, then by dominant refID
alignments$queryID <- factor(alignments$queryID,
                             levels = unique(as.character(alignments$queryID)))

queryMaxAlnIdx <- tapply(alignments$lenAln, alignments$queryID, which.max, simplify = FALSE)
alignments$queryID <- factor(alignments$queryID,
  levels = unique(as.character(alignments$queryID))[order(mapply(
    function(x, id) alignments$refStart2[alignments$queryID == id][x],
    queryMaxAlnIdx, names(queryMaxAlnIdx)
  ))]
)

queryLenAggPerRef <- sapply(levels(alignments$queryID), function(qid)
  tapply(alignments$lenAln[alignments$queryID == qid],
         alignments$refID[ alignments$queryID == qid], sum))

queryID_Ref <- if (length(levels(alignments$refID)) > 1)
  apply(queryLenAggPerRef, 2, function(x) rownames(queryLenAggPerRef)[which.max(x)])
else
  sapply(queryLenAggPerRef, function(x) names(queryLenAggPerRef)[which.max(x)])

alignments$queryID <- factor(alignments$queryID,
  levels = levels(alignments$queryID)[order(match(queryID_Ref, levels(alignments$refID)))])

queryMax    <- tapply(alignments$queryLen, alignments$queryID, max)
queryRevComp <- character(0)

if (opt$flip) {
  is_rev <- tapply(alignments$queryEnd - alignments$queryStart, alignments$queryID, sum) < 0
  queryRevComp <- names(is_rev)[is_rev]
  rev_rows <- alignments$queryID %in% queryRevComp
  qmax_matched <- queryMax[as.character(alignments$queryID[rev_rows])]
  alignments$queryStart[rev_rows] <- qmax_matched - alignments$queryStart[rev_rows] + 1
  alignments$queryEnd[rev_rows]   <- qmax_matched - alignments$queryEnd[rev_rows]   + 1
}

queryOffset <- cumulative_offsets(queryMax)
alignments$queryStart2 <- alignments$queryStart + queryOffset[as.character(alignments$queryID)]
alignments$queryEnd2   <- alignments$queryEnd   + queryOffset[as.character(alignments$queryID)]

# ─────────────────────────────────────────────
#  Calculate Statistics
# ─────────────────────────────────────────────

stats <- calc_alignment_stats(alignments)

# ─────────────────────────────────────────────
#  Build Plot
# ─────────────────────────────────────────────

plot_title <- if (!is.null(opt$plot_title)) opt$plot_title else paste(script_name, "- Genome Alignment Dot Plot")

gp <- ggplot(alignments) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(angle = 15, hjust = 1),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    legend.position    = "right",
    plot.caption       = element_text(size = 8, color = "grey50", hjust = 0),
    plot.title         = element_text(size = 14, face = "bold", margin = margin(b = 10))
  )

# Color palette selection
if (opt$color_palette == "RdYlBu") {
  color_scale <- scale_color_distiller(palette = "RdYlBu", direction = 1, limits = c(ifelse(opt$identity_floor > 0, opt$identity_floor, 0), 1), name = "Identity")
} else if (opt$color_palette == "Viridis") {
  color_scale <- scale_color_viridis_c(option = "viridis", limits = c(ifelse(opt$identity_floor > 0, opt$identity_floor, 0), 1), name = "Identity")
} else if (opt$color_palette == "Plasma") {
  color_scale <- scale_color_viridis_c(option = "plasma", limits = c(ifelse(opt$identity_floor > 0, opt$identity_floor, 0), 1), name = "Identity")
} else if (opt$color_palette == "Heat") {
  color_scale <- scale_color_gradient(low = "blue", high = "red", limits = c(ifelse(opt$identity_floor > 0, opt$identity_floor, 0), 1), name = "Identity")
} else {
  color_scale <- scale_color_distiller(palette = "RdYlBu", direction = 1, limits = c(ifelse(opt$identity_floor > 0, opt$identity_floor, 0), 1), name = "Identity")
}

gp <- gp + color_scale

caption_text <- paste0(
  script_name, " v", script_version, "  |  ",
  "alignments: ", nrow(alignments), "  |  ",
  "queries: ", length(unique(alignments$queryID)), "  |  ",
  "coverage: ", round(stats$total_bp_aligned / 1e6, 1), "Mb  |  ",
  "avg identity: ", round(stats$avg_identity, 1), "%"
)

if (opt$show_stats) {
  caption_text <- paste0(caption_text, "\n",
    "median identity: ", round(stats$median_identity, 1), "%  |  ",
    "min aln: ", opt$min_align, "  |  ",
    "min query aln: ", opt$min_query_aln)
}

gp <- gp + labs(
  title   = plot_title,
  caption = caption_text
)

# ── X axis (reference) ──────────────────────
if (length(levels(alignments$refID)) == 1) {
  reflen  <- unique(alignments$refLen)
  xbreaks <- seq(0, reflen, length.out = 11)
  step    <- diff(xbreaks)[1]
  gp <- gp +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, reflen + 1),
                       breaks = xbreaks,
                       labels = fmt_bp(xbreaks, step)) +
    xlab(unique(as.character(alignments$refID)))
} else {
  gp <- gp +
    theme(panel.grid.major.x = element_blank()) +
    geom_vline(xintercept = cumsum(as.numeric(chromMax)),
               color = "#cccccc", linewidth = 0.4) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, sum(as.numeric(chromMax)) + 1),
      breaks = cumsum(as.numeric(chromMax)) - chromMax / 2,
      labels = substr(levels(alignments$refID), 1, 20)
    ) +
    xlab("Reference")
}

# ── Y axis (query) ───────────────────────────
if (length(levels(alignments$queryID)) == 1) {
  queryLen <- unique(alignments$queryLen)
  ybreaks  <- seq(0, queryLen, length.out = 11)
  step     <- diff(ybreaks)[1]
  gp <- gp +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, queryLen + 1),
                       breaks = ybreaks,
                       labels = fmt_bp(ybreaks, step)) +
    ylab(unique(as.character(alignments$queryID)))
} else {
  gp <- gp +
    theme(panel.grid.major.y = element_blank()) +
    geom_hline(yintercept = cumsum(as.numeric(queryMax)),
               color = "#cccccc", linewidth = 0.4) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, sum(as.numeric(queryMax)) + 1),
      breaks = cumsum(as.numeric(queryMax)) - queryMax / 2,
      labels = substr(levels(alignments$queryID), 1, 20)
    ) +
    ylab("Query")
}

# ── Alignment segments ───────────────────────
gp <- gp + geom_segment(
  aes(x = refStart2, xend = refEnd2,
      y = queryStart2, yend = queryEnd2,
      color = percentID),
  linewidth = opt$line_width,
  alpha     = opt$alpha
)

# ── Break points ─────────────────────────────
if (opt$break_point) {
  pt_size <- opt$plot_size / 60
  gp <- gp +
    geom_point(aes(x = refStart2, y = queryStart2, color = percentID),
               size = pt_size, shape = 19) +
    geom_point(aes(x = refEnd2,   y = queryEnd2,   color = percentID),
               size = pt_size, shape = 19)
}

# ── BED annotations ─────────────────────────
if (!is.null(opt$ref_bed_file))
  gp <- add_bed_lines(gp, opt$ref_bed_file, chromMax, axis = "x")

if (!is.null(opt$query_bed_file))
  gp <- add_bed_lines(gp, opt$query_bed_file, queryMax, axis = "y",
                      flip_ids = queryRevComp)

# ─────────────────────────────────────────────
#  Save Output
# ─────────────────────────────────────────────

w <- opt$plot_size
h <- opt$plot_size * 0.8

suppressWarnings({
  if (opt$output_format %in% c("pdf", "both"))
    ggsave(paste0(opt$output_filename, ".pdf"), plot = gp,
           width = w, height = h, units = "in", dpi = 300, limitsize = FALSE)
  if (opt$output_format %in% c("png", "both"))
    ggsave(paste0(opt$output_filename, ".png"), plot = gp,
           width = w, height = h, units = "in", dpi = 300, limitsize = FALSE)
  if (opt$output_format == "svg")
    ggsave(paste0(opt$output_filename, ".svg"), plot = gp,
           width = w, height = h, units = "in", dpi = 300, limitsize = FALSE)
})

output_formats <- if (opt$output_format == "both") "PDF + PNG" 
                  else toupper(opt$output_format)

cat(sprintf("\n✓ %s completed successfully!\n", script_name))
cat(sprintf("Output saved: %s [%s]\n", opt$output_filename, output_formats))
cat(sprintf("Total alignments: %d | Coverage: %.1f Mb | Avg identity: %.1f%%\n", 
            stats$total_alignments, stats$total_bp_aligned / 1e6, stats$avg_identity))
cat("GitHub: https://github.com/ank-man/genodot\n")
