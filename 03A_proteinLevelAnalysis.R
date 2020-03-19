source("/Users/isabell/Desktop/projects/THP_SEC_SWATH_MS/THP_CCprofiler_scripts/THP_header.R")

#' # Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
design_matrix <- readRDS("design_matrix.rds")
calibration <- readRDS("calibration.rds")

#' Set output directory
outdir = "proteinCentric"
if (! dir.exists(outdir)) {
  dir.create(outdir)
}
setwd(outdir)

#' # Integrate peptide traces across conditions
pepTracesIntegrated_sum <- integrateTraceIntensities(pepTracesList_filtered,
                                                     design_matrix = NULL,
                                                     integrate_within = NULL,
                                                     aggr_fun = "sum")

#' # Perform protein-centric analysis
# Suppress messages for clarity
suppressMessages(
  proteinFeatures <-  findProteinFeatures(traces = pepTracesIntegrated_sum,
                                          corr_cutoff = 0.9,
                                          window_size = 7,
                                          parallelized = T,
                                          n_cores = 3,
                                          collapse_method = "apex_only",
                                          perturb_cutoff = "5%",
                                          rt_height = 1,
                                          smoothing_length = 7)
)
saveRDS(proteinFeatures,"proteinFeatures.rda")

#' Score protein features and filter for a detection FDR of 5%
proteinFeaturesFiltered <- scoreFeatures(proteinFeatures, FDR=0.05, PDF=TRUE, name="qvalueStats_proteinFeatures")
saveRDS(proteinFeaturesFiltered,"proteinFeaturesFiltered.rda")

#' Evaluate protein feature stats
summarizeFeatures(proteinFeaturesFiltered,
                  plot=TRUE,
                  PDF=TRUE,
                  name="feature_summary_proteinFeaturesFiltered")

#' ## Extract feature values
proteinFeaturesFiltered_featureVals <- extractFeatureVals(traces = pepTracesList_filtered,
                                             features = proteinFeaturesFiltered,
                                             design_matrix = design_matrix,
                                             extract = "subunits_detected",
                                             imputeZero = T,
                                             verbose = F)
saveRDS(proteinFeaturesFiltered_featureVals, "proteinFeaturesFiltered_featureVals.rda")

#' ## Fill feature values
proteinFeaturesFiltered_featureValsFilled <- fillFeatureVals(featureVals = proteinFeaturesFiltered_featureVals,
                                                tracesList = pepTracesList_filtered,
                                                design_matrix = design_matrix)
saveRDS(proteinFeaturesFiltered_featureValsFilled, "proteinFeaturesFiltered_featureValsFilled.rda")
