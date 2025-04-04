# Script to prepare sample test data for the package
# This script should not be included in the built package

library(Seurat)
library(SeuratObject)

# Create a small example Seurat object
set.seed(123)
example_counts <- matrix(rpois(20000, lambda = 0.1), nrow = 100, ncol = 200)
rownames(example_counts) <- paste0("Gene", 1:100)
colnames(example_counts) <- paste0("Cell", 1:200)

# Create Seurat object
example_seurat <- CreateSeuratObject(counts = example_counts)

# Add some metadata
example_seurat$sample <- "example"
example_seurat$condition <- sample(c("A", "B"), ncol(example_seurat), replace = TRUE)

# Run standard preprocessing
example_seurat <- NormalizeData(example_seurat)
example_seurat <- FindVariableFeatures(example_seurat)
example_seurat <- ScaleData(example_seurat)
example_seurat <- RunPCA(example_seurat, npcs = 10)
example_seurat <- RunUMAP(example_seurat, dims = 1:10)
example_seurat <- FindNeighbors(example_seurat, dims = 1:10)
example_seurat <- FindClusters(example_seurat, resolution = 0.5)

# Save to inst/extdata
dir.create("inst/extdata/sample", recursive = TRUE, showWarnings = FALSE)
saveRDS(example_seurat, file = "inst/extdata/sample/example_seurat.rds", compress = TRUE, version = 2)

# Create a small 10X-like counts matrix for testing raw data import
sparse_counts <- Matrix::Matrix(
  rpois(5000, 0.1), 
  nrow = 50, 
  ncol = 100, 
  sparse = TRUE
)
rownames(sparse_counts) <- paste0("Gene", 1:50)
colnames(sparse_counts) <- paste0("Cell", 1:100)

# Create output directory
raw_dir <- "inst/extdata/sample/raw_example"
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# Write files in 10X format
Matrix::writeMM(sparse_counts, file = file.path(raw_dir, "matrix.mtx"))
write.table(
  data.frame(gene = rownames(sparse_counts), name = rownames(sparse_counts)),
  file = file.path(raw_dir, "features.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
  data.frame(cell = colnames(sparse_counts)),
  file = file.path(raw_dir, "barcodes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Indicate that the script has completed
message("Test data preparation complete. Files saved to inst/extdata/sample/")