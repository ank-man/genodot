# Installation Guide

## System Requirements

- **R** (>= 4.0.0)
- **Unix-like system** (Linux, macOS, Windows with WSL)
- **2GB RAM** minimum for small genomes, **8GB+** recommended for large genomes
- **1GB disk space** for R packages

## Quick Installation

### 1. Install R Dependencies

```bash
# Install required R packages
R -e "install.packages(c('ggplot2', 'optparse', 'scales', 'RColorBrewer'))"
```

### 2. Download GenoDot

```bash
# Option 1: Direct download
wget https://raw.githubusercontent.com/genomics/genodot/main/genodot.R
chmod +x genodot.R

# Option 2: Clone repository
git clone https://github.com/genomics/genodot.git
cd genodot
chmod +x genodot.R
```

### 3. Test Installation

```bash
# Download test data
wget https://raw.githubusercontent.com/genomics/genodot/main/examples/sample.paf

# Run GenoDot
./genodot.R sample.paf

# Should output: sample.paf.pdf and sample.paf.png
```

## Detailed Installation

### R Installation

#### Ubuntu/Debian
```bash
sudo apt update
sudo apt install r-base r-base-dev
```

#### CentOS/RHEL
```bash
sudo yum install R
# or on newer systems
sudo dnf install R
```

#### macOS
```bash
# Using Homebrew
brew install r

# Or download from CRAN
# https://cran.r-project.org/bin/macosx/
```

#### Windows
1. Download R from https://cran.r-project.org/bin/windows/base/
2. Install R with default settings
3. Use WSL or Git Bash for command-line usage

### Package Installation

#### From CRAN
```bash
R -e "install.packages(c('ggplot2', 'optparse', 'scales', 'RColorBrewer'))"
```

#### From Source (if needed)
```bash
R -e "install.packages('ggplot2', type='source')"
R -e "install.packages('optparse', type='source')"
R -e "install.packages('scales', type='source')"
R -e "install.packages('RColorBrewer', type='source')"
```

## Docker Installation

```bash
# Pull GenoDot Docker image
docker pull genodot/genodot:latest

# Run GenoDot
docker run -v $(pwd):/data genodot/genodot genodot.R /data/alignments.paf
```

## Conda Installation

```bash
# Create conda environment
conda create -n genodot r-base r-ggplot2 r-optparse r-scales r-rcolorbrewer
conda activate genodot

# Download GenoDot
wget https://raw.githubusercontent.com/genomics/genodot/main/genodot.R
chmod +x genodot.R
```

## Troubleshooting

### Common Issues

#### 1. "Permission denied"
```bash
chmod +x genodot.R
```

#### 2. "R not found"
Install R following the system-specific instructions above.

#### 3. "Package not found"
```bash
R -e "install.packages('PACKAGE_NAME', repos='https://cran.rstudio.com/')"
```

#### 4. Memory issues with large genomes
```bash
# Increase R memory limit
export R_MAX_VSIZE=16Gb
./genodot.R large_genome.paf
```

### Performance Optimization

#### For Large Genomes (>1Gb)
```bash
# Use conservative filtering
./genodot.R -m 50kb -q 1mb -c 0.9 large_genome.paf

# Reduce plot size for faster rendering
./genodot.R -p 10 -F png large_genome.paf
```

#### For HPC Clusters
```bash
# Submit as batch job
#!/bin/bash
#SBATCH --mem=8G
#SBATCH --time=02:00:00
module load R
./genodot.R -p 20 -F pdf alignments.paf
```

## Verification

Test your installation with the provided example:

```bash
# Download test data
curl -O https://raw.githubusercontent.com/genomics/genodot/main/examples/sample.paf

# Run GenoDot
./genodot.R -o test_output sample.paf

# Check output
ls -la test_output.*
# Should show test_output.pdf and test_output.png
```

## Getting Help

- 📖 [Documentation](https://github.com/genomics/genodot/wiki)
- 🐛 [Issues](https://github.com/genomics/genodot/issues)
- 💬 [Discussions](https://github.com/genomics/genodot/discussions)
- 📧 support@genodot.org
