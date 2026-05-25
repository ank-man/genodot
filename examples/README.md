# GenoDot Examples

This directory contains example files and outputs to demonstrate GenoDot's capabilities.

## Files

### 📊 **sample.paf**
Sample PAF file with simulated genome alignments between human chromosomes.
- **28 alignments** across **24 query sequences**
- **Coverage**: 1.4 Mb
- **Identity**: 100% (simulated perfect alignments)
- **Reference**: Human chromosomes (chr1-chrY)
- **Format**: Standard PAF with 12 columns

### 🎨 **GenoDot_demo.pdf / GenoDot_demo.png**
Demo dot plot generated from sample.paf using GenoDot with the following parameters:
```bash
./genodot.R -q 10kb -m 1kb -r 100kb -p 12 -C Viridis -S -t "GenoDot Demo - Human Genome Alignment" -o GenoDot_demo sample.paf
```

**Features demonstrated:**
- Viridis color palette for identity visualization
- Detailed statistics display
- Professional formatting with title
- Multi-chromosome reference layout
- High-resolution output (300 DPI)

## Usage Examples

### Basic Visualization
```bash
# Generate basic dot plot
./genodot.R examples/sample.paf

# Generate with custom title
./genodot.R -t "My Analysis" examples/sample.paf
```

### Advanced Visualization
```bash
# Reproduce the demo plot
./genodot.R -q 10kb -m 1kb -r 100kb -p 12 -C Viridis -S -t "GenoDot Demo" -o demo examples/sample.paf

# Different color palette
./genodot.R -C Plasma -p 15 examples/sample.paf

# With break points
./genodot.R -b examples/sample.paf
```

### Format Conversion Examples
```bash
# Convert to BED format
./paf2bed.R examples/sample.paf

# Convert to SAM format
./paf2sam.R examples/sample.paf

# Convert back from SAM to PAF
./sam2paf.R examples/sample.sam
```

## Expected Output

When you run GenoDot on the sample data, you should see:
```
Reading PAF file...
Alignments read:    28
Query sequences:    24
After filtering:    28 alignments | 24 queries

✓ GenoDot completed successfully!
Output saved: examples/sample [PDF + PNG]
Total alignments: 28 | Coverage: 1.4 Mb | Avg identity: 100.0%
```

## Integration with Other Tools

The sample data can be used to test integration with:
- **BEDTools**: `bedtools intersect -a examples/sample.bed -b regions.bed`
- **Samtools**: `samtools view -bS examples/sample.sam | samtools sort -o sorted.bam`
- **IGV**: Load `examples/sample.sam` for visualization
- **UCSC Genome Browser**: Load `examples/sample.bed` as custom track

## Performance Testing

Use the sample data to benchmark:
```bash
time ./genodot.R examples/sample.paf
time ./paf2bed.R examples/sample.paf
time ./paf2sam.R examples/sample.paf
```

## Troubleshooting

If you encounter issues with the sample data:

1. **Check file permissions**: `chmod +x examples/sample.paf`
2. **Validate PAF format**: `head -n 5 examples/sample.paf`
3. **Check R packages**: `R -e "library(ggplot2); library(optparse)"`

## Creating Your Own Examples

To create similar test data:
```bash
# Generate synthetic PAF from BED
echo -e "chr1\t10000\t60000\tregion1\t1000\t+\nchr2\t20000\t70000\tregion2\t1000\t+" > test.bed
./bed2paf.R test.bed

# Or use Minimap2 with real data
minimap2 -x asm5 reference.fasta query.fasta > your_data.paf
```

---

*These examples showcase GenoDot's powerful visualization and conversion capabilities for genomics research.*
