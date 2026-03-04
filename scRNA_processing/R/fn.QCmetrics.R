# R/fn.QCmetrics.R
# QC metric annotation for Seurat objects
# - percent.mt
# - percent.ribo
# - percent.hb
# - log10GenesPerUMI

add_qc_metrics <- function(sdata, cap = FALSE) {

  mt.pattern <- ifelse(
    cap,
    "^GRCh38-MT-|^mm10---mt-",
    "^MT-|^mt-"
  )

  ribo.pattern <- ifelse(
    cap,
    "^GRCh38-RP[SL]|^mm10---Rp[sl]",
    "^RP[SL]|^Rp[sl]"
  )

  # HB genes only at gene start (exclude HBP naturally)
  hb.pattern <- ifelse(
    cap,
    "^GRCh38-HB[A-Z]|^mm10---Hb[a-z]",
    "^HB[A-Z]|^Hb[a-z]"
  )

  sdata[["percent.mt"]]   <- PercentageFeatureSet(sdata, pattern = mt.pattern)
  sdata[["percent.ribo"]] <- PercentageFeatureSet(sdata, pattern = ribo.pattern)
  sdata[["percent.hb"]]   <- PercentageFeatureSet(sdata, pattern = hb.pattern)

  sdata[["log10GenesPerUMI"]] <-
    log10(sdata$nFeature_RNA) / log10(sdata$nCount_RNA)

  return(sdata)
}

