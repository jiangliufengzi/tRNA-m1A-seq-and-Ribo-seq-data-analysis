# tRNA-seq Analysis Pipeline

## Project Overview

This project contains a comprehensive bioinformatics pipeline for analyzing tRNA sequencing (tRNA-seq) and ribosome profiling (Ribo-seq) data. The primary focus is investigating the impact of tRNA modifications (particularly m1A modifications) on translation efficiency and changes in codon usage bias under different gene knockout conditions.

## Research Background

This study focuses on:
- **tRNA m1A Modification**: Analyzing changes in N1-methyladenosine (m1A) modification levels in tRNA
- **Translation Efficiency**: Calculating differential translation efficiency of genes through Ribo-seq data
- **Codon Usage Bias**: Studying changes in codon usage frequencies under different conditions
- **Gene Expression Regulation**: Exploring the relationship between tRNA modifications and gene expression regulation

## File Structure and Functions

### Python Scripts

#### Core Analysis Tools
- **`bam_file_cds_codon_usage.py`** - Calculate codon usage frequency in CDS regions from BAM files
- **`bam2cov.py`** - Convert BAM/SAM files to coverage files, supporting both count and coverage modes
- **`longest_cds_transcripts.py`** - Extract longest transcript information and CDS sequences from genes
- **`gene_codon_usage.py`** - Calculate codon usage frequency for specific genes

#### Data Processing Tools
- **`filter_cov_by_pause.py`** - Filter coverage data based on ribosome pause sites
- **`get_ribosome_pause_site.py`** - Identify and analyze ribosome pause sites
- **`reverse_complement.py`** - DNA sequence reverse complement functionality

### R Scripts

#### Visualization Functions
- **`plot_diff_tRNA_exp.R`** - Generate volcano plots for differential tRNA expression
- **`plot_tRNA_mutation.R`** - Visualize tRNA m1A modification level changes
- **`read_m1a.R`** - Read and process m1A modification data
- **`pvalue.R`** - Statistical significance testing tools
- **`reverse_complement.R`** - R version of sequence reverse complement functionality

### Analysis Workflow Documentation

#### Ribo-seq Analysis (`wl-ribo-seq.Rmd`)
1. **Translation Efficiency Calculation**
   - Calculate differential translation efficiency using Xtail algorithm
   - Filter protein-coding genes for analysis
   - Generate cumulative distribution plots (CDF) of translation efficiency
   - Create translation efficiency volcano plots

2. **Codon Usage Analysis**
   - Extract longest transcript information
   - Calculate codon usage frequency for up-regulated/down-regulated genes
   - Analyze relationship between codon usage and m1A modifications
   - Generate codon usage preference plots

3. **Specific Gene Analysis**
   - Analyze codon usage of key genes (e.g., fech, nrf1, ndufb7)
   - Compare codon preference differences between different genes

4. **Functional Enrichment Analysis**
   - GO, KEGG, Reactome, and WikiPathways enrichment analysis
   - GSEA gene set enrichment analysis

#### tRNA-seq Analysis (`wl-tRNA-seq.Rmd`)
1. **Data Preprocessing**
   - Quality control and adapter trimming
   - UMI processing and sequence orientation correction
   - Sequence alignment (Bowtie and BWA)

2. **Expression Analysis**
   - tRNA expression counting and normalization
   - Differential expression analysis
   - Merging of identical codon variants

3. **Modification Analysis**
   - m1A modification site detection
   - Modification level quantification
   - Statistical analysis of modification rate differences

4. **Comprehensive Analysis**
   - Association analysis between tRNA expression and modification levels
   - Comparison across different treatment conditions
   - Codon-specific modification patterns

## Key Features

### Data Processing Capabilities
- **Multi-format Support**: Supports various bioinformatics standard formats including BAM, SAM, FASTA, GTF
- **High-throughput Processing**: Optimized algorithms suitable for large-scale sequencing data analysis
- **Quality Control**: Comprehensive data quality control and filtering pipelines

### Analysis Functions
- **Translation Efficiency Calculation**: Precise Ribo-seq translation efficiency analysis
- **Codon Analysis**: Comprehensive codon usage bias analysis
- **Modification Detection**: Accurate tRNA modification site identification and quantification
- **Statistical Analysis**: Rigorous statistical testing and multiple comparison correction

### Visualization Functions
- **Volcano Plots**: Intuitive display of differential expression and translation efficiency changes
- **Box Plots**: Statistical distribution and group comparisons
- **Heatmaps**: Visualization of modification levels and expression profiles
- **Circular Plots**: Codon usage preference display

## Requirements

### Python Dependencies
```python
pandas
numpy
biopython
pyfaidx
pysam
gffutils
```

### R Dependencies
```r
ggplot2
dplyr
tidyr
readxl
viridis
clusterProfiler
org.Dr.eg.db
xtail
```

### System Requirements
- Python ≥ 3.7
- R ≥ 4.0
- Sufficient memory for processing large genomic datasets
- Multi-threading capable CPU

## Quick Start

### 1. Translation Efficiency Analysis
```python
from longest_cds_transcripts import enhanced_transcript_analysis
from bam_file_cds_codon_usage import output_codon

# Extract longest transcripts
enhanced_transcript_analysis(
    by="cds",
    fasta_file="genome.fasta",
    db_file="annotation.db"
)

# Calculate codon usage
output_codon(
    input_dir="coverage_files/",
    gene_list=None,
    fasta_file="genome.fasta"
)
```

### 2. tRNA Modification Analysis
```r
# Load analysis scripts
source("plot_tRNA_mutation.R")
source("read_m1a.R")

# Plot modification level changes
plot_tRNA_mutation(
    input_file="tRNA_data.xlsx",
    sheet=6,
    treat="treatment_group"
)
```

### 3. Coverage Analysis
```python
from bam2cov import calculate_rna_coverage

# Calculate RNA coverage
calculate_rna_coverage(
    input_file="alignment.bam",
    output_file="coverage.txt",
    type="coverage"
)
```

## Output Results

### Main Output Files
- **Codon Usage Tables**: Excel format codon frequency statistics
- **Translation Efficiency Results**: Differential translation efficiency tables with statistical tests
- **Modification Level Data**: tRNA modification quantification results
- **Visualization Charts**: PDF format analysis plots and charts

### Statistical Reports
- Differential expression gene lists
- Functional enrichment analysis results
- Quality control reports
- Parameter and threshold setting records

## Data Sources

This pipeline is applicable to the following types of sequencing data:
- **Ribo-seq**: Ribosome-bound mRNA fragment sequencing
- **RNA-seq**: Total RNA or mRNA sequencing
- **tRNA-seq**: tRNA-specific sequencing
- **Reference Genome**: Supports genome annotations for various model organisms

## Citation and References

If you use this pipeline in your research, please cite the relevant methodological papers:

- Xtail algorithm for translation efficiency calculation
- Bioinformatics methods for tRNA modification detection
- Codon usage bias analysis methods

## Technical Support

This project focuses on bioinformatics analysis of tRNA and translation regulation, providing a complete analysis pipeline from raw data to biological interpretation. All scripts have been optimized to ensure accuracy and reproducibility of analysis results.

## Update Log

- Initial version: Complete tRNA-seq and Ribo-seq analysis pipeline
- Feature optimization: Improved data processing efficiency and statistical analysis methods
- Visualization enhancement: Added multiple chart types and customization options

---

**Note**: Please ensure all dependency packages are properly installed before use, and adjust corresponding parameters according to your data characteristics.
