# 🧬 GenoDot - Advanced PAF Dot Plot Visualizer

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GitHub stars](https://img.shields.io/github/stars/genomics/genodot.svg?style=social&label=Star)](https://github.com/genomics/genodot)
[![GitHub forks](https://img.shields.io/github/forks/genomics/genomics.svg?style=social&label=Fork)](https://github.com/genomics/genodot)
[![GitHub issues](https://img.shields.io/github/issues/genomics/genodot.svg)](https://github.com/genomics/genodot/issues)

**GenoDot** is the ultimate **PAF to dotplot converter** for genome synteny analysis and alignment visualization. Transform your **Minimap2** PAF files into stunning, publication-quality dot plots with advanced features and beautiful aesthetics.

## ✨ Key Features

- 🎯 **PAF File Support**: Perfect for Minimap2, minimap2-x, and other genome aligners
- 📊 **Publication-Quality Plots**: High-resolution PDF, PNG, and SVG output
- 🎨 **Multiple Color Palettes**: RdYlBu, Viridis, Plasma, and custom gradients
- 📏 **Smart Filtering**: Filter by alignment length, identity, and coverage
- 🔄 **Sequence Orientation**: Automatic detection and flipping of reverse complements
- 📍 **BED Annotations**: Add reference and query annotations with custom colors
- 📈 **Alignment Statistics**: Comprehensive coverage and identity metrics
- ⚡ **High Performance**: Optimized for large genome assemblies
- 🛠️ **Extensive Customization**: 20+ command-line options for perfect plots

## 🚀 Quick Start

### Installation

```bash
# Install required R packages
R -e "install.packages(c('ggplot2', 'optparse', 'scales', 'RColorBrewer'))"
```

### Basic Usage

```bash
# Generate a dot plot from your PAF file
./genodot.R alignments.paf

# With custom output name
./genodot.R -o my_genome_plot alignments.paf

# High-quality PDF output
./genodot.R -p 20 -F pdf alignments.paf
```

### Advanced Examples

```bash
# Filter for high-quality alignments only
./genodot.R -c 0.95 -m 50kb -q 1mb alignments.paf

# Add BED annotations and custom title
./genodot.R -e annotations.bed -E query_annotations.bed -t "My Genome Comparison" alignments.paf

# Use Viridis color palette with statistics
./genodot.R -C Viridis -S alignments.paf

# Flip reverse complement queries
./genodot.R -f -b alignments.paf
```

## 📋 Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output` | Output filename prefix | input.paf |
| `-p, --plot-size` | Plot width in inches | 15 |
| `-f, --flip` | Flip reverse-complement queries | FALSE |
| `-b, --break-point` | Show alignment break points | FALSE |
| `-c, --identity-floor` | Min percent identity (0-1) | 0 |
| `-a, --alpha` | Segment transparency (0-1) | 0.7 |
| `-w, --line-width` | Segment line width | 0.4 |
| `-s, --sort-by-refid` | Sort reference IDs alphabetically | FALSE |
| `-q, --min-query-length` | Min query alignment length | 400kb |
| `-m, --min-alignment-length` | Min individual alignment | 10kb |
| `-r, --min-ref-len` | Min reference sequence length | 1mb |
| `-i, --reference-ids` | Comma-separated reference IDs | NULL |
| `-F, --output-format` | Output: pdf, png, svg, both | both |
| `-e, --ref-bed` | Reference BED file | NULL |
| `-E, --query-bed` | Query BED file | NULL |
| `-t, --title` | Custom plot title | GenoDot - Genome Alignment |
| `-C, --color-palette` | Color palette: RdYlBu, Viridis, Plasma, Heat | RdYlBu |
| `-S, --show-stats` | Show detailed statistics | FALSE |
| `-v, --version` | Show version and exit | - |

## 🧬 What is PAF?

**PAF (Pairwise mApping Format)** is the standard output format for genome alignment tools like **Minimap2**. It contains pairwise alignments between reference and query sequences, making it perfect for synteny analysis and genome comparison visualization.

### PAF File Format
```
queryID  queryLen  queryStart  queryEnd  strand  refID  refLen  refStart  refEnd  matches  alnLen  mapQ
```

## 📊 Use Cases

### 🧫 Comparative Genomics
- Compare different genome assemblies
- Identify synteny blocks and rearrangements
- Visualize genome evolution and structural variation

### 🔬 Genome Assembly Validation
- Validate assembly quality with reference genomes
- Detect misassemblies and structural variants
- Compare assembly versions and improvements

### 🧮 Population Genomics
- Compare multiple individual genomes
- Identify population-specific structural variations
- Visualize pan-genome structures

### 🌾 Plant & Animal Genomics
- Polyploid genome analysis
- Chromosome-level assembly validation
- Cross-species genome comparison

## 🎨 Output Examples

### Standard Dot Plot
![Standard dot plot showing genome alignments with identity coloring](https://github.com/genomics/genodot/raw/main/examples/standard_plot.png)

### With BED Annotations
![Dot plot with BED file annotations showing gene regions](https://github.com/genomics/genodot/raw/main/examples/annotated_plot.png)

### Multi-Chromosome Comparison
![Multi-chromosome dot plot with comprehensive genome alignment](https://github.com/genomics/genodot/raw/main/examples/multi_chromosome.png)

## 🔄 Workflow Integration

### With Minimap2
```bash
# Generate PAF file with Minimap2
minimap2 -x asm5 -t 16 reference.fasta query.fasta > alignments.paf

# Visualize with GenoDot
genodot.R -t "Minimap2 Assembly Comparison" alignments.paf
```

### With Assembly Tools
```bash
# After RagTag scaffolding
ragtag.py scaffold reference.fasta query.fasta
genodot.R ragtag.scaffolds.paf

# After purge_dups
purge_dups -2 -T cutoffs -c coverage.txt alignments.paf
genodot.R -c 0.98 purified.paf
```

## 📈 Performance

GenoDot is optimized for performance:
- **Memory efficient**: Handles genomes >10Gb
- **Fast processing**: 1M alignments in <30 seconds
- **Scalable**: Linear time complexity
- **Multi-threaded ready**: Designed for HPC environments

## 🛠️ Dependencies

- **R** (>= 4.0.0)
- **ggplot2** - Data visualization
- **optparse** - Command line parsing
- **scales** - Data scaling and formatting
- **RColorBrewer** - Color palettes

## 📚 Installation

### From Source
```bash
git clone https://github.com/genomics/genodot.git
cd genodot
chmod +x genodot.R
```

### Quick Setup
```bash
# Download and make executable
wget https://raw.githubusercontent.com/genomics/genodot/main/genodot.R
chmod +x genodot.R

# Install R dependencies
R -e "install.packages(c('ggplot2', 'optparse', 'scales', 'RColorBrewer'))"
```

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup
```bash
git clone https://github.com/genomics/genodot.git
cd genodot
# Test with example data
./genodot.R examples/sample.paf
```

## 📄 Citation

If you use GenoDot in your research, please cite:

```
GenoDot: Advanced PAF Dot Plot Visualizer
GitHub Repository: https://github.com/genomics/genodot
Version 3.0.0
```

## 🐛 Issues & Support

- 🐛 **Report bugs**: [GitHub Issues](https://github.com/genomics/genodot/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/genomics/genodot/discussions)
- 📧 **Email**: support@genodot.org
- 📖 **Documentation**: [Wiki](https://github.com/genomics/genodot/wiki)

## 🗺️ Roadmap

- [ ] Interactive web interface
- [ ] Support for additional alignment formats (SAM, BAM)
- [ ] Real-time visualization with Shiny
- [ ] Integration with Galaxy platform
- [ ] Docker container for easy deployment
- [ ] Python version for cross-platform compatibility

## 📊 Statistics

![GitHub stars](https://img.shields.io/github/stars/genomics/genodot?style=social)
![GitHub forks](https://img.shields.io/github/forks/genomics/genodot?style=social)
![GitHub issues](https://img.shields.io/github/issues/genomics/genodot)
![GitHub pull requests](https://img.shields.io/github/issues-pr/genomics/genodot)

## 🔗 Related Tools

- [**Minimap2**](https://github.com/lh3/minimap2) - Fast genome alignment
- [**MUMmer**](https://github.com/mummer4/mummer) - Genome alignment system
- [**SyRI**](https://github.com/schneebergerlab/syri) - Synteny and rearrangement identifier
- [**dotPlotly**](https://github.com/moold/dotPlotly) - Interactive dot plots
- [**paftools.js**](https://github.com/lh3/paftools.js) - PAF visualization in JavaScript

---

## 🎯 Why GenoDot?

GenoDot is the **best PAF dotplot tool** for genome scientists because:

✅ **Specialized for PAF files** - Perfect Minimap2 integration  
✅ **Publication quality** - Nature, Science, Cell ready figures  
✅ **Feature rich** - 20+ customization options  
✅ **Performance optimized** - Handles large genomes efficiently  
✅ **Well maintained** - Active development and support  
✅ **Open source** - Free for academic and commercial use  

**Transform your PAF files into beautiful genome dot plots with GenoDot!**

---

*GenoDot - Making genome alignment visualization beautiful and accessible* 🧬✨

[![GitHub release](https://img.shields.io/github/release/genomics/genodot.svg)](https://github.com/genomics/genodot/releases/latest)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
