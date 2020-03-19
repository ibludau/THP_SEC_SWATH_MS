source("/Users/isabell/Desktop/projects/THP_SEC_SWATH_MS/THP_CCprofiler_scripts/THP_header.R")

#' # Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
protein_traces_list <- readRDS("protein_traces_list.rds")
design_matrix <- readRDS("design_matrix.rds")
calibration <- readRDS("calibration.rds")

setwd("complexCentric")
complex_featureValsFilled <- readRDS("complex_featureValsFilled.rda")
scoredDataAll <- readRDS("complexFeaturesCollapsed.rda") 

library(ggrepel)

n_conditions <- length(unique(design_matrix$Condition))

for (i in seq(1,n_conditions,1)) {
  conditions <- unique(design_matrix$Condition)
  conditions <- conditions[-i]
  sub_design_matrix <- subset(design_matrix, Condition %in% conditions)
  sub_complex_featureValsFilled <- subset(complex_featureValsFilled, Condition %in% conditions)
  sub_pepTracesList_filtered <- pepTracesList_filtered[sub_design_matrix$Sample_name]
  class(sub_pepTracesList_filtered) <- "tracesList"
  sub_protein_traces_list <- protein_traces_list[sub_design_matrix$Sample_name]
  class(sub_protein_traces_list) <- "tracesList"

  #' ## Perform peptide-level differential expression testing for all features
  complex_DiffExprPep <- testDifferentialExpression(featureVals = sub_complex_featureValsFilled,
                                                    compare_between = "Condition",
                                                    level = "peptide",
                                                    measuredOnly = FALSE)

  saveRDS(complex_DiffExprPep, paste0("complex_DiffExprPep_",conditions[1],"_",conditions[2],".rda"))

  pdf(paste0("complex_DiffExprPep_stats_",conditions[1],"_",conditions[2],".pdf"), width = 4, height=4)
    hist(complex_DiffExprPep$pVal, breaks = 100)
    hist(complex_DiffExprPep$global_pVal, breaks = 100)
    hist(log(complex_DiffExprPep$qint1), breaks = 100)
    hist(log(complex_DiffExprPep$qint2), breaks = 100)
    hist(log(complex_DiffExprPep$global_int2_imp), breaks = 100)
    hist(log(complex_DiffExprPep$global_int1_imp), breaks = 100)
    hist(complex_DiffExprPep$log2FC, breaks = 50)
    hist(complex_DiffExprPep$global_log2FC, breaks = 50)
  dev.off()

  #' ## Aggregate differential expression results to the protein level
  complex_DiffExprProtein <- aggregatePeptideTests(complex_DiffExprPep)
  saveRDS(complex_DiffExprProtein, paste0("complex_DiffExprProtein_",conditions[1],"_",conditions[2],".rda"))

  pdf(paste0("complex_DiffExprProtein_stats_",conditions[1],"_",conditions[2],".pdf"), width = 4, height=4)
    hist(complex_DiffExprProtein$pVal, breaks = 100)
    hist(complex_DiffExprProtein$global_pVal, breaks = 100)
    hist(complex_DiffExprProtein$medianLog2FC, breaks = 50)
    hist(complex_DiffExprProtein$global_medianLog2FC, breaks = 50)
  dev.off()

  #' ## Aggregate differential expression results to the complex level
  complex_DiffExprComplex <- aggregateProteinTests(complex_DiffExprProtein)
  saveRDS(complex_DiffExprComplex, paste0("complex_DiffExprComplex_",conditions[1],"_",conditions[2],".rda"))

  pdf(paste0("complex_DiffExprComplex_stats_",conditions[1],"_",conditions[2],".pdf"), width = 4, height=4)
    hist(complex_DiffExprComplex$pVal, breaks = 100)
    hist(complex_DiffExprComplex$global_pVal, breaks = 100)
    hist(complex_DiffExprComplex$medianLog2FC, breaks = 50)
    hist(complex_DiffExprComplex$global_medianLog2FC, breaks = 50)
  dev.off()

  #' # Make volcano plots
  #' General volcano plots
  plotVolcano(complex_DiffExprPep, PDF = T, pBHadj_cutoff = 0.05, name = paste0("complex_DiffExprPep_",conditions[1],"_",conditions[2]))
  plotVolcano(complex_DiffExprProtein, PDF = T, pBHadj_cutoff = 0.05, name = paste0("complex_DiffExprProtein_",conditions[1],"_",conditions[2]))
  plotVolcano(complex_DiffExprComplex, PDF = T, pBHadj_cutoff = 0.05, name = paste0("complex_DiffExprComplex_",conditions[1],"_",conditions[2]))

  #' Get differential complexes
  allComplexes <- unique(complex_DiffExprComplex$complex_id)
  diffComplexes <- unique(complex_DiffExprComplex[pBHadj<0.05][abs(medianLog2FC)>1]$complex_id)
  write.table(allComplexes, paste0("allComplexes_",conditions[1],"_",conditions[2],".txt"), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(diffComplexes, paste0("diffComplexes_",conditions[1],"_",conditions[2],".txt"), col.names = F, row.names = F, quote = F, sep = "\t")

  #' Volcanoplts highlighting different complexes
  highlight_complex <- diffComplexes[1]
  plotVolcano(complex_DiffExprComplex, highlight=highlight_complex, PDF = T, pBHadj_cutoff = 0.05, name = paste0("complex_DiffExprComplex_",highlight_complex,"_",conditions[1],"_",conditions[2]))

  #' # Plot some differential protein features
  pdf(paste0("diffComplexes_",conditions[1],"_",conditions[2],".pdf"), width=8, height=7)
  for (id in diffComplexes) {
    plotFeatures(feature_table = scoredDataAll,
                 traces = sub_protein_traces_list,
                 feature_id = id,
                 annotation_label = "Entry_name",
                 design_matrix=sub_design_matrix,
                 calibration=calibration,
                 peak_area = T,
                 legend = F,
                 onlyBest = F,
                 monomer_MW=T)
  }
  dev.off()

  pdf(paste0("diffComplexes_aggregated_",conditions[1],"_",conditions[2],".pdf"), width=8, height=7)
  for (id in diffComplexes) {
    plotFeatures(feature_table = scoredDataAll,
                 traces = sub_protein_traces_list,
                 feature_id = id,
                 annotation_label = "Entry_name",
                 design_matrix=sub_design_matrix,
                 calibration=calibration,
                 peak_area = T,
                 legend = F,
                 onlyBest = F,
                 monomer_MW=T,
                 aggregateReplicates=T)
  }
  dev.off()
}
