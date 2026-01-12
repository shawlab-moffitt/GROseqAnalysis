# GRO/PRO-seq Global Pausing Index Analysis

This repository provides an R script to compute a **global transcriptional pausing index** from **single-cell GRO-seq / PRO-seq data**, integrate normalized gene expression (illustrated using **Malat1**), and evaluate the relationship between transcriptional pausing and expression at the **per-cell level**.

The workflow is designed around data structured like **GSE242176**, but can be adapted to other GRO/PRO-seq experiments with minimal modification.

---

## Overview

The analysis performs the following steps:

1. Load GRO-seq / PRO-seq count matrices and consolidated read alignments  
2. Construct unique per-cell identifiers  
3. Normalize GRO-seq counts (counts-per-10k)  
4. Define **strand-aware promoter and gene-body regions** (mm10)  
5. Compute a **global pausing index per cell**  
6. Merge Malat1 expression with pausing metrics  
7. Quantify and visualize the correlation between pausing and expression  

---

## Global Pausing Index

For each cell, the global pausing index is defined as:

\[
\text{Global Pausing Index} =
\frac{\text{Promoter read density}}{\text{Gene body read density}}
\]

Where:
- **Promoter density** = reads overlapping all promoters ÷ total promoter bp  
- **Gene body density** = reads overlapping all gene bodies ÷ total gene body bp  

This yields a **cell-level summary of transcriptional pausing**, not a gene-specific pausing index.

---

## Repository Structure

