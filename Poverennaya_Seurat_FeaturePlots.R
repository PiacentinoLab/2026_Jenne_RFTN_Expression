# Load required libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(reticulate)

# Install Python packages (only needed first time)
py_install("anndata", pip = TRUE)
py_install("scanpy", pip = TRUE)

# Import Python modules
ad <- import("anndata")
sc <- import("scanpy")

# Read the h5ad file
adata <- ad$read_h5ad("Poverennaya_2026_mouse_cranial_dataset_E8-E10.NC_and_dorsal_NT_lineages.h5ad")

# Manually extract data from AnnData object
# Get the count matrix
counts <- t(py_to_r(adata$X))

# Get gene names and cell barcodes
gene_names <- py_to_r(adata$var_names$to_list())
cell_names <- py_to_r(adata$obs_names$to_list())

# Set row and column names
rownames(counts) <- gene_names
colnames(counts) <- cell_names

# Get metadata (cell annotations)
metadata <- py_to_r(adata$obs)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, 
                                 meta.data = metadata,
                                 project = "Poverennaya_NC")
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj <- NormalizeData(seurat_obj)


print("Checking for common embeddings:")

# Check for UMAP
tryCatch({
  umap_coords <- py_to_r(adata$obsm["X_umap"])
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  rownames(umap_coords) <- cell_names
  seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = umap_coords, 
                                               key = "UMAP_", 
                                               assay = "RNA")
  print("UMAP coordinates added!")
}, error = function(e) {
  print(paste("No UMAP found:", e$message))
})

# Check for PCA
tryCatch({
  pca_coords <- py_to_r(adata$obsm["X_pca"])
  colnames(pca_coords) <- paste0("PC_", 1:ncol(pca_coords))
  rownames(pca_coords) <- cell_names
  seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = pca_coords, 
                                              key = "PC_", 
                                              assay = "RNA")
  print("PCA coordinates added!")
}, error = function(e) {
  print(paste("No PCA found:", e$message))
})

# Check for tSNE
tryCatch({
  tsne_coords <- py_to_r(adata$obsm["X_tsne"])
  colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")
  rownames(tsne_coords) <- cell_names
  seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = tsne_coords, 
                                               key = "tSNE_", 
                                               assay = "RNA")
  print("tSNE coordinates added!")
}, error = function(e) {
  print(paste("No tSNE found:", e$message))
})

# Explore the object
print("Available reductions in Seurat object:")
print(names(seurat_obj@reductions))
print(seurat_obj)
print(paste("Number of cells:", ncol(seurat_obj)))
print(paste("Number of features:", nrow(seurat_obj)))

# Check available metadata columns
print("Available metadata columns:")
print(colnames(seurat_obj@meta.data))

# Print cluster distribution
print("Cluster distribution:")
print(table(seurat_obj$seurat_clusters))

# If no UMAP exists, create one
if (!"umap" %in% names(seurat_obj@reductions)) {
  print("No UMAP found, computing one...")
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
}

# ====== VISUALIZATION SECTION ======

# 1. UMAP by cluster
DimPlot(seurat_obj, 
        group.by = "celltype",
        pt.size = 0.5) +
  ggtitle("UMAP by Cluster") +
  #NoLegend()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 2. Feature plots for specific genes
FeaturePlot(seurat_obj, 
            features = "Rftn2",
            pt.size = 0.8,
            order = TRUE,
            cols = c("lightgrey", "red")) +  # Add red color scale
  #NoLegend()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
