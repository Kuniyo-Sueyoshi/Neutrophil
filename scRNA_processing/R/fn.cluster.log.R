# R/helper_seurat_log.R
# Single-sample clustering pipeline (LogNormalize + clustering)

fn.cluster.log <- function(
  data,
  vars.to.regress = c("percent.mt", "nCount_RNA"),
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  nfeatures = 3000,
  npcs = 50,
  dims = 1:30,
  resolutions = c(0.2, 0.4)
) {

  data <- Seurat::NormalizeData(
    data,
    normalization.method = normalization.method,
    scale.factor = scale.factor,
    verbose = FALSE
  )

  data <- Seurat::FindVariableFeatures(
    data,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )

  data <- Seurat::ScaleData(
    data,
    features = Seurat::VariableFeatures(data),
    vars.to.regress = vars.to.regress,
    verbose = FALSE
  )

  data <- Seurat::RunPCA(
    data,
    features = Seurat::VariableFeatures(data),
    npcs = npcs,
    verbose = FALSE
  )

  data <- Seurat::RunUMAP(data, reduction = "pca", dims = dims, verbose = FALSE)
  data <- Seurat::FindNeighbors(data, reduction = "pca", dims = dims, verbose = FALSE)
  data <- Seurat::FindClusters(data, resolution = resolutions, verbose = FALSE)

  return(data)
}
