source("/Users/isabell/Desktop/projects/THP_SEC_SWATH_MS/THP_CCprofiler_scripts/THP_header.R")

#' # Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
protein_traces_list <- readRDS("protein_traces_list.rds")
design_matrix <- readRDS("design_matrix.rds")
calibration <- readRDS("calibration.rds")

#' Set output directory
outdir = "complexCentric"
if (! dir.exists(outdir)) {
  dir.create(outdir)
}
setwd(outdir)

#' ## Summarize protein traces across conditions
protein_traces_sum <- integrateTraceIntensities(protein_traces_list,
                                                design_matrix = NULL,
                                                integrate_within = NULL,
                                                aggr_fun = "sum")

#' ## Save protein traces
saveRDS(protein_traces_sum,"protein_traces_sum.rds")

stringTargets <- readRDS("../../StringComplexHypothesisGeneration/stringComplexes_clusterOne_final.rds")
stringTargets[, complex_id := paste0(complex_id,"_string")]
stringTargets[, complex_name := paste0(complex_name,"_string")]
corumTargetes <- corumComplexHypothesesRedundant
corumTargetes[, complex_id := paste0(complex_id,"_corum")]

combiTargets <- rbind(stringTargets,corumTargetes)
saveRDS(combiTargets,"combiTargets.rds")

plotSummarizedMScoverage(hypotheses = combiTargets,
                         protTraces = protein_traces_sum,
                         PDF = TRUE)

binaryCombiTargets <- generateBinaryNetwork(combiTargets)

pathLengthCombiTargets <- calculatePathlength(binaryCombiTargets)

combiTargetsPlusDecoys <- generateComplexDecoys(target_hypotheses = combiTargets,
                                                 dist_info = pathLengthCombiTargets,
                                                 min_distance = 3,
                                                 append = TRUE)

saveRDS(combiTargetsPlusDecoys,"combiTargetsPlusDecoys.rds")

#' ## Complex feature finding
complexFeatures <-  findComplexFeatures(traces = protein_traces_sum,
                                        complex_hypothesis = combiTargetsPlusDecoys,
                                        corr_cutoff = 0.9,
                                        window_size = 7,
                                        parallelized = F,
                                        n_cores = 1,
                                        collapse_method = "apex_network",
                                        perturb_cutoff = "5%",
                                        rt_height = 1,
                                        smoothing_length = 7)

saveRDS(complexFeatures,"combiComplexFeatures.rda")
# complexFeatures <- readRDS("combiComplexFeatures.rda")

complexFeatures <- filterFeatures(complexFeatures,
                                  complex_ids = NULL,
                                  protein_ids = NULL,
                                  min_feature_completeness = NULL,
                                  min_hypothesis_completeness = NULL,
                                  min_subunits = NULL,
                                  min_peak_corr = NULL,
                                  min_monomer_distance_factor = 2)

hypothesis <- "combiTargets"

plotSummarizedMScoverage(hypotheses = combiTargets, protein_traces_sum, PDF = T)

complexFeaturesBest <- getBestFeatures(complexFeatures)
saveRDS(complexFeaturesBest,"complexFeaturesBest.rda")

complexFeaturesBestFiltered <- scoreFeatures(complexFeaturesBest, FDR=0.05, PDF=T, name="qvalueStats_complexFeatures")

scoredDataAll <- appendSecondaryComplexFeatures(scoredPrimaryFeatures = complexFeaturesBestFiltered, allFeatures = complexFeatures, peakCorr_cutoff = 0.7)
saveRDS(scoredDataAll,"scoredDataAll.rds")

plotSummarizedComplexes(scoredDataAll, combiTargets, protein_traces_sum, PDF=T, name="complex_completeness_pie")

summarizeFeatures(scoredDataAll,
                plot=TRUE,
                PDF=T,
                name="feature_summary_scoredDataAll")

summarizeFeatures(complexFeaturesBestFiltered,
                  plot=TRUE,
                  PDF=T,
                  name="feature_summary_complexFeaturesBestFiltered")


complexFeaturesUnique <- getUniqueFeatureGroups(
                            feature_table = scoredDataAll,
                            rt_height = 0,
                            distance_cutoff = 1.25)

complexFeaturesCollapsed <- callapseByUniqueFeatureGroups(
                              feature_table = complexFeaturesUnique,
                              rm_decoys = TRUE)
saveRDS(complexFeaturesCollapsed, "complexFeaturesCollapsed.rda")

#' ## Extract feature values
complex_featureVals <- extractFeatureVals(traces = pepTracesList_filtered,
                                          features = complexFeaturesCollapsed,
                                          design_matrix = design_matrix,
                                          extract = "subunits",
                                          imputeZero = T,
                                          verbose = F)
saveRDS(complex_featureVals, "complex_featureVals.rda")

complex_featureValsFilled <- fillFeatureVals(featureVals = complex_featureVals,
                                             tracesList = pepTracesList_filtered,
                                             design_matrix = design_matrix)
saveRDS(complex_featureValsFilled, "complex_featureValsFilled.rda")
