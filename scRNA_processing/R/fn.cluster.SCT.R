
# Sample integration / clustering pipeline (SCTransform + optional Harmony)

fn.cluster.SCT <- function(
  data,
  batch = "batch",                         # sample/condition column name
  do_harmony = TRUE,                       # TRUE for sample integration
  vars.to.regress = c("G2M.Score", "S.Score", "percent.mt"),
  nfeatures = 3000,
  npcs = 50,
  dims = 1:30,
  resolutions = c(0.2, 0.4, 0.8, 1.2, 2, 3)
) {

  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
  if (!requireNamespace("sctransform", quietly = TRUE)) stop("sctransform is required.")

  if (do_harmony) {
    if (!requireNamespace("harmony", quietly = TRUE)) stop("harmony is required when do_harmony=TRUE.")
    if (!batch %in% colnames(data@meta.data)) {
      stop(sprintf('Batch column "%s" not found in meta.data.', batch))
    }
  }

  use_glmGamPoi <- requireNamespace("glmGamPoi", quietly = TRUE)
  sct_method <- if (use_glmGamPoi) "glmGamPoi" else "regular"

  so <- Seurat::SCTransform(
    object = data,
    method = sct_method,
    do.correct.umi = TRUE,
    vars.to.regress = vars.to.regress,
    variable.features.n = nfeatures,
    verbose = FALSE
  )
  Seurat::DefaultAssay(so) <- "SCT"

  so <- Seurat::RunPCA(so, assay = "SCT", npcs = npcs, verbose = FALSE)

  reduction_use <- "pca"
  if (do_harmony) {
    so <- Seurat::RunHarmony(
      object = so,
      group.by.vars = batch,
      assay.use = "SCT",
      dims.use = dims,
      reduction = "pca"
    )
    reduction_use <- "harmony"
  }

  so <- Seurat::RunUMAP(so, reduction = reduction_use, dims = dims, assay = "SCT", verbose = FALSE)
  so <- Seurat::FindNeighbors(so, reduction = reduction_use, dims = dims, verbose = FALSE)
  so <- Seurat::FindClusters(so, resolution = resolutions, verbose = FALSE)
  so <- Seurat::PrepSCTFindMarkers(so, assay = "SCT")

  return(so)
}
  