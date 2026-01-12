############################################################
## Load required libraries
############################################################

library(Seurat)  # Used here mainly for handling single-cell style objects


############################################################
## Load GRO-seq / PRO-seq count matrices and metadata
############################################################

# Load aggregated GRO-seq counts (AGT, collapsed counts)
GSE242176_AGTuC125cv2p8_counts <- readRDS(
  "/Users/4472414/Projects/EricPadron/GROseqExample/GSE242176_AGTuC125cv2p8_mapq3qc_filtered_counts.rds"
)

# Load single-cell GRO-seq counts
GSE242176_scGROv2p8_counts <- readRDS(
  "/Users/4472414/Projects/EricPadron/GROseqExample/GSE242176_scGROv2p8_mapq3qc_filtered_counts.rds"
)

# Load consolidated single-cell GRO-seq reads (GRanges with metadata)
GSE242176_scGROv2p8_consolidated <- readRDS(
  "/Users/4472414/Projects/EricPadron/GROseqExample/GSE242176_scGROv2p8_consolidated.rds"
)

# Construct a unique cell identifier from experimental metadata
GSE242176_scGROv2p8_consolidated$cell_id <- paste(
  GSE242176_scGROv2p8_consolidated$Exp,
  GSE242176_scGROv2p8_consolidated$Plate,
  GSE242176_scGROv2p8_consolidated$Cell,
  sep = "-"
)

# Inspect the consolidated object
GSE242176_scGROv2p8_consolidated


############################################################
## Inspect count matrices
############################################################

# Preview the single-cell GRO-seq count matrix
head(as.matrix(GSE242176_scGROv2p8_counts))

# Load PRO-seq counts (bulk-style)
GSE242176_PROv2p8_counts <- readRDS(
  "/Users/4472414/Projects/EricPadron/GROseqExample/GSE242176_PROv2p8_mapq3qc_filtered_counts.rds"
)

# Load consolidated PRO-seq reads
GSE242176_PROv2p8_consolidated <- readRDS(
  "/Users/4472414/Projects/EricPadron/GROseqExample/GSE242176_PROv2p8_consolidated.rds"
)

# Inspect PRO-seq metadata
head(GSE242176_PROv2p8_consolidated)


############################################################
## Normalize GRO-seq counts to CPM-like units
############################################################

# Convert sparse matrix to base R matrix
GRO_counts <- as.matrix(GSE242176_scGROv2p8_counts)

# Compute library size per cell (column sums)
GRO_lib_size <- colSums(GRO_counts)

# Normalize counts to counts-per-10k (CP10K-style normalization)
GRO_norm <- sweep(GRO_counts, 2, GRO_lib_size, FUN = "/") * 1e4

# Inspect distribution of normalized values
hist(log2(GRO_norm))

# Identify rows corresponding to Malat1
rownames(GRO_norm)[grepl("Malat1", rownames(GRO_norm))]

# Extract Malat1 signal across all cells
GRO_norm["GN-Malat1", ]

# Inspect PRO-seq count matrix
head(as.matrix(GSE242176_PROv2p8_counts))


############################################################
## Calculate the GRO-seq Global Pausing Index
############################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(AnnotationDbi)
})

# Load mouse mm10 gene annotations from UCSC knownGene
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Extract gene ranges (one range per gene, named by Entrez ID)
genes_gr <- genes(txdb)

############################################################
## Define promoter and gene-body parameters
############################################################

prom_up   <- 1000  # upstream of TSS
prom_down <- 150   # downstream of TSS

body_start_offset <- 500  # skip immediate promoter-proximal region
body_end_trim     <- 0    # optional trimming near TES

# Generate strand-aware promoter regions
prom_gr <- promoters(
  genes_gr,
  upstream = prom_up,
  downstream = prom_down
)

############################################################
## Define strand-aware gene-body regions
############################################################

# Determine transcription start site (TSS) and end site (TES)
tss <- ifelse(strand(genes_gr) == "-", end(genes_gr), start(genes_gr))
tes <- ifelse(strand(genes_gr) == "-", start(genes_gr), end(genes_gr))

# Compute gene-body start and end coordinates
body_start <- ifelse(
  strand(genes_gr) == "-",
  tes + body_end_trim,
  tss + body_start_offset
)

body_end <- ifelse(
  strand(genes_gr) == "-",
  tss - body_start_offset,
  tes - body_end_trim
)

# Build gene-body GRanges object
body_gr <- GRanges(
  seqnames = seqnames(genes_gr),
  ranges   = IRanges(
    start = pmin(body_start, body_end),
    end   = pmax(body_start, body_end)
  ),
  strand   = strand(genes_gr),
  gene_id  = names(genes_gr)
)

# Remove invalid or extremely short gene bodies
body_gr <- body_gr[width(body_gr) > 0]

# Annotate promoters with gene IDs
prom_gr$gene_id <- names(prom_gr)


############################################################
## Pausing index calculation settings
############################################################

ignore_strand <- FALSE   # GRO/PRO-seq is strand-specific
pseudo <- 1e-9           # pseudocount to avoid divide-by-zero

# Ensure cell_id metadata exists
stopifnot("cell_id" %in% colnames(mcols(GSE242176_scGROv2p8_consolidated)))

# Drop reads without valid cell IDs
GSE242176_scGROv2p8_consolidated <- 
  GSE242176_scGROv2p8_consolidated[
    !is.na(mcols(GSE242176_scGROv2p8_consolidated)$cell_id)
  ]

cell_id <- as.character(
  mcols(GSE242176_scGROv2p8_consolidated)$cell_id
)


############################################################
## Compute total promoter and gene-body base pairs
############################################################

# Reduce overlapping regions to avoid double-counting
prom_bp <- sum(width(reduce(prom_gr)))
body_bp <- sum(width(reduce(body_gr)))


############################################################
## Count GRO-seq reads overlapping promoters
############################################################

hits_p <- findOverlaps(
  prom_gr,
  GSE242176_scGROv2p8_consolidated,
  ignore.strand = ignore_strand
)

cells_p <- cell_id[subjectHits(hits_p)]
prom_counts_by_cell <- table(cells_p)


############################################################
## Count GRO-seq reads overlapping gene bodies
############################################################

hits_b <- findOverlaps(
  body_gr,
  GSE242176_scGROv2p8_consolidated,
  ignore.strand = ignore_strand
)

cells_b <- cell_id[subjectHits(hits_b)]
body_counts_by_cell <- table(cells_b)


############################################################
## Assemble global pausing index per cell
############################################################

# Union of all observed cells
all_cells <- union(
  names(prom_counts_by_cell),
  names(body_counts_by_cell)
)

# Align promoter and body counts
p_counts <- as.numeric(
  prom_counts_by_cell[match(all_cells, names(prom_counts_by_cell))]
)
b_counts <- as.numeric(
  body_counts_by_cell[match(all_cells, names(body_counts_by_cell))]
)

# Replace missing values with zero
p_counts[is.na(p_counts)] <- 0
b_counts[is.na(b_counts)] <- 0

# Compute read densities
prom_density <- p_counts / (prom_bp + pseudo)
body_density <- b_counts / (body_bp + pseudo)

# Global pausing index = promoter density / gene-body density
global_pause_index <- prom_density / (body_density + pseudo)

# Final per-cell pausing index table
global_pi_df <- data.frame(
  cell_id = all_cells,
  promoter_reads = p_counts,
  body_reads = b_counts,
  promoter_bp = prom_bp,
  body_bp = body_bp,
  promoter_density = prom_density,
  body_density = body_density,
  global_pause_index = global_pause_index
)

# Inspect most paused cells
global_pi_df[order(-global_pi_df$global_pause_index), ][1:10, ]


############################################################
## Merge Malat1 expression with pausing index
############################################################

# Extract normalized Malat1 vector
malat1_vec <- GRO_norm["GN-Malat1", ]
malat1_vec <- as.numeric(malat1_vec)

# Build Malat1 expression data frame
malat1_df <- data.frame(
  cell_id = colnames(GRO_norm),
  Malat1_norm = malat1_vec
)

# Merge with pausing index table
global_pi_df <- merge(
  global_pi_df,
  malat1_df,
  by = "cell_id",
  all.x = TRUE
)

head(global_pi_df)


############################################################
## Correlation analysis and visualization
############################################################

library(dplyr)
library(ggplot2)

# Remove missing or zero Malat1 values
df <- global_pi_df %>%
  filter(!is.na(Malat1_norm), Malat1_norm > 0)

# Log2-transform Malat1 expression
df <- df %>%
  mutate(log2_Malat1_norm = log2(Malat1_norm))

# Pearson correlation
r_pearson <- cor(
  df$log2_Malat1_norm,
  df$global_pause_index,
  method = "pearson"
)

# Statistical test
cor_test <- cor.test(
  df$log2_Malat1_norm,
  df$global_pause_index,
  method = "pearson"
)

# Scatter plot with linear regression
p <- ggplot(df, aes(x = log2_Malat1_norm, y = global_pause_index)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  labs(
    x = "log2(Malat1 counts per 10k reads)",
    y = "Global pausing index",
    title = sprintf(
      "Global pausing index vs Malat1 (r = %.3f, p = %.2e)",
      r_pearson,
      cor_test$p.value
    )
  ) +
  theme_classic()

p
