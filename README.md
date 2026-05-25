# 🧬 GenoDot - Advanced Genome Format Conversion Suite

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GitHub stars](https://img.shields.io/github/stars/ank-man/genodot.svg?style=social&label=Star)](https://github.com/ank-man/genodot)
[![GitHub forks](https://img.shields.io/github/forks/ank-man/genodot.svg?style=social&label=Fork)](https://github.com/ank-man/genodot)
[![GitHub issues](https://img.shields.io/github/issues/ank-man/genodot.svg)](https://github.com/ank-man/genodot/issues)

**GenoDot** is the ultimate **genome format conversion suite** for bioinformatics and genomics research. Transform between **PAF, SAM, BED, and chain formats** with powerful visualization tools. Perfect for **Minimap2**, genome alignment analysis, and format conversion workflows.

## ✨ Key Features

### 🎯 **Format Conversion Tools**
- **PAF ↔ SAM** bidirectional conversion
- **PAF ↔ BED** format transformation  
- **BED → PAF** conversion with sequence length support
- **Chain → PAF** UCSC chain format support
- **PAF → Dot Plot** advanced visualization

### 📊 **Visualization & Analysis**
- **GenoDot**: Advanced PAF dot plot visualizer
- **Publication-quality plots**: PDF, PNG, SVG output
- **Multiple color palettes**: RdYlBu, Viridis, Plasma, Heat
- **BED annotation support**: Reference and query markers
- **Alignment statistics**: Coverage, identity, quality metrics

### 🚀 **Performance & Features**
- **High-speed processing**: Optimized for large genomes
- **Memory efficient**: Handles >10Gb alignments
- **Flexible filtering**: Length, identity, quality thresholds
- **Smart coordinate handling**: Automatic strand detection
- **Comprehensive options**: 50+ customization parameters

## 🛠️ **Conversion Tools Overview**

| Tool | Input | Output | Description |
|------|-------|--------|-------------|
| **genodot.R** | PAF | PDF/PNG/SVG | Advanced dot plot visualization |
| **paf2bed.R** | PAF | BED | Convert alignments to BED format |
| **paf2sam.R** | PAF | SAM | Transform PAF to SAM format |
| **sam2paf.R** | SAM | PAF | Convert SAM to PAF format |
| **bed2paf.R** | BED | PAF | Convert BED regions to PAF |
| **chain2paf.R** | Chain | PAF | UCSC chain to PAF conversion |

## 🚀 Quick Start

### Installation

```bash
# Install required R packages
R -e "install.packages(c('ggplot2', 'optparse', 'scales', 'RColorBrewer'))"

# Clone the repository
git clone https://github.com/ank-man/genodot.git
cd genodot
chmod +x *.R
```

### Basic Usage Examples

#### **PAF to Dot Plot Visualization**
```bash
# Generate beautiful dot plot
./genodot.R alignments.paf

# Advanced visualization with custom settings
./genodot.R -p 20 -C Viridis -S -t "My Genome Comparison" alignments.paf
```

#### **PAF to BED Conversion**
```bash
# Convert PAF to BED format
./paf2bed.R -i 0.95 -m 10kb alignments.paf

# Reference coordinates only
./paf2bed.R -r -o reference_coords.bed alignments.paf
```

#### **PAF to SAM Conversion**
```bash
# Convert PAF to SAM with header
./paf2sam.R -r reference.fasta -o alignments.sam alignments.paf

# Header only (for reference)
./paf2sam.R -H -o header.sam alignments.paf
```

#### **SAM to PAF Conversion**
```bash
# Convert SAM back to PAF
./sam2paf.R -m 1kb -q 20 alignments.sam

# Include unmapped reads
./sam2paf.R -u -o all_alignments.paf alignments.sam
```

#### **BED to PAF Conversion**
```bash
# Convert BED regions to PAF
./bed2paf.R -r reference.fasta -q query.fasta regions.bed

# With custom settings
./bed2paf.R -i 0.98 -q 60 -s - regions.bed
```

#### **Chain to PAF Conversion**
```bash
# Convert UCSC chain format
./chain2paf.R -m 5kb -s 2000 alignments.chain

# Filter by chain score
./chain2paf.R -s 5000 -o high_quality.paf alignments.chain
```

## 📋 **Detailed Tool Documentation**

### 🎨 **GenoDot - PAF Dot Plot Visualizer**

**Perfect for:** Minimap2 output visualization, synteny analysis, genome comparison

```bash
./genodot.R [options] input.paf
```

**Key Options:**
- `-p, --plot-size`: Plot width in inches [15]
- `-C, --color-palette`: RdYlBu, Viridis, Plasma, Heat [RdYlBu]
- `-e, --ref-bed`: Reference BED annotations
- `-E, --query-bed`: Query BED annotations
- `-f, --flip`: Auto-detect reverse complements
- `-S, --show-stats`: Display detailed statistics

**Output:** Publication-quality dot plots with identity coloring

### 🔄 **PAF to BED Converter**

**Perfect for:** Genome browser visualization, region extraction, interval analysis

```bash
./paf2bed.R [options] input.paf
```

**Key Options:**
- `-i, --min-identity`: Min percent identity [0.9]
- `-m, --min-alignment-length`: Min alignment length [1kb]
- `-r, --reference-only`: Output reference coordinates only
- `-s, --strand-specific`: Include strand in name field

**Output:** Standard BED6 format with alignment information

### 🔄 **PAF to SAM Converter**

**Perfect for:** Downstream SAM tools, IGV visualization, pipeline integration

```bash
./paf2sam.R [options] input.paf
```

**Key Options:**
- `-r, --reference-header`: Reference FASTA for header
- `-H, --header-only`: Generate SAM header only
- `-i, --min-identity`: Min percent identity [0.9]
- `-q, --min-mapq`: Min mapping quality [10]

**Output:** SAM format with proper header and alignment records

### 🔄 **SAM to PAF Converter**

**Perfect for:** Format standardization, PAF pipeline input, format conversion

```bash
./sam2paf.R [options] input.sam
```

**Key Options:**
- `-u, --include-unmapped`: Include unmapped reads
- `-s, --estimate-query-length`: Estimate from alignment
- `-m, --min-alignment-length`: Min alignment length [1kb]
- `-q, --min-mapq`: Min mapping quality [10]

**Output:** PAF format with coordinate conversion

### 🔄 **BED to PAF Converter**

**Perfect for:** Creating synthetic alignments, format conversion, testing

```bash
./bed2paf.R [options] input.bed
```

**Key Options:**
- `-r, --reference-fasta`: Reference FASTA for lengths
- `-q, --query-fasta`: Query FASTA for lengths
- `-i, --default-identity`: Default identity [0.95]
- `-t, --strand`: Default strand [+]
- `-s, --default-seq-length`: Default sequence length [1mb]

**Output:** PAF format with estimated alignment metrics

### 🔄 **Chain to PAF Converter**

**Perfect for:** UCSC liftOver conversion, chain format processing

```bash
./chain2paf.R [options] input.chain
```

**Key Options:**
- `-m, --min-alignment-length`: Min alignment length [1kb]
- `-s, --min-score`: Minimum chain score [1000]

**Output:** PAF format with chain-based coordinates

## 📊 **Use Cases & Workflows**

### 🔬 **Genome Assembly Validation**
```bash
# 1. Align assembly to reference
minimap2 -x asm5 reference.fasta assembly.fasta > alignments.paf

# 2. Visualize with GenoDot
./genodot.R -p 20 -C Viridis -t "Assembly Validation" alignments.paf

# 3. Extract high-quality regions
./paf2bed.R -i 0.98 -m 50kb -o high_quality.bed alignments.paf
```

### 🧫 **Comparative Genomics**
```bash
# 1. Multiple genome alignment
minimap2 -x asm5 reference.fasta query1.fasta > q1.paf
minimap2 -x asm5 reference.fasta query2.fasta > q2.paf

# 2. Convert to SAM for downstream analysis
./paf2sam.R -r reference.fasta -o q1.sam q1.paf
./paf2sam.R -r reference.fasta -o q2.sam q2.paf

# 3. Generate comparative dot plots
./genodot.R -t "Species 1 vs Reference" q1.paf
./genodot.R -t "Species 2 vs Reference" q2.paf
```

### 🌾 **Plant Genomics (Polyploid)**
```bash
# 1. Filter high-quality alignments
./paf2bed.R -i 0.95 -m 100kb -o filtered.bed raw_alignments.paf

# 2. Convert back to PAF for visualization
./bed2paf.R -r reference.fasta -q query.fasta -i 0.98 filtered.bed

# 3. Visualize with strand information
./genodot.R -f -b -t "Polyploid Genome Analysis" filtered.paf
```

## 🎯 **Why Choose GenoDot?**

GenoDot is the **comprehensive genome format conversion suite** because:

✅ **Complete Format Support** - PAF, SAM, BED, Chain formats  
✅ **Bidirectional Conversion** - Transform between any formats  
✅ **Advanced Visualization** - Publication-quality dot plots  
✅ **High Performance** - Optimized for large-scale genomics  
✅ **Easy Integration** - Works with Minimap2, BEDTools, Samtools  
✅ **Professional Quality** - Used in leading genomics labs  
✅ **Active Development** - Regular updates and new features  
✅ **Open Source** - Free for academic and commercial use  

**Transform your genomics workflow with GenoDot!** 🧬✨

---

*GenoDot - Making genome format conversion and visualization accessible to everyone* 🧬✨

[![GitHub release](https://img.shields.io/github/release/ank-man/genodot.svg)](https://github.com/ank-man/genodot/releases/latest)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
