
# Gut Microbiota Analysis in EAM Mice

This repository contains an R script (`Part2.R`) for analyzing 16S rRNA sequencing data to study gut microbiota dynamics in a mouse model of experimental autoimmune myocarditis (EAM). It includes processing, filtering, normalization, and visualization of taxonomic data, along with statistical testing and differential abundance analysis.

## 🔬 Main Analyses

- **Data Input & Preprocessing**  
  Loads DADA2 output, metadata, and constructs a phyloseq object with ASV table, taxonomy, metadata, representative sequences, and a phylogenetic tree.

- **Filtering & Normalization**  
  Filters out unwanted taxa (e.g., chloroplasts, mitochondria), and applies relative abundance transformation.

- **Taxonomic Composition Visualization**  
  Generates stacked bar plots by Family, Genus, and Species across microbiome groups and timepoints using `ggplot2`.

- **Alpha Diversity (Shannon, Chao1)**  
  Calculates and compares alpha diversity indices across groups with statistical testing (`wilcox.test`).

- **Beta Diversity (Bray-Curtis PCoA)**  
  Performs ordination using Bray-Curtis distances and visualizes using PCoA plots.

- **PERMANOVA**  
  Performs pairwise PERMANOVA tests for community structure differences across groups and timepoints.

- **Differential Abundance Analysis (Log2 Fold Change)**  
  Uses DESeq2 to identify differentially abundant taxa across microbiome groups and timepoints. Outputs significant ASVs by log2 fold change and adjusted p-value.

- **Focused Genus Analysis**  
  Plots the abundance of specific taxa (e.g., *Oscillibacter*) across treatment groups.

## 📁 Input Files

- `dada2_Lacto2_nochim.fa` — Representative sequences  
- `dada2_Lacto2_nochim.tree` — Phylogenetic tree  
- `metadata_file.tsv` — Sample metadata  
- `EAM_Lacto2_nochim.rda` — Preprocessed data

## 📦 R Packages Used

- `phyloseq`, `ggplot2`, `DESeq2`, `vegan`, `pairwiseAdonis`, `gridExtra`, `ggpubr`, `Biostrings`, `tidyverse`

## 📈 Output

- Publication-ready barplots by Family, Genus, Species
- PCoA plots showing beta diversity across groups
- Boxplots for alpha diversity (Shannon, Chao1)
- Volcano-style plots of DESeq2 log2 fold change
- CSV/TSV tables of significant taxa

## 🧑‍🔬 Author

Phoebe (Xu) Shi  
PhD Candidate, University of Nebraska–Lincoln  
Email: xshi9@huskers.unl.edu  
[LinkedIn](https://www.linkedin.com/in/phoebe-xu-shi/)
