source("/Users/isabell/Desktop/projects/THP_SEC_SWATH_MS/THP_CCprofiler_scripts/THP_header.R")

#' # Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
design_matrix <- readRDS("design_matrix.rds")
calibration <- readRDS("calibration.rds")

setwd("proteinCentric")
proteinFeaturesFiltered_featureValsFilled <- readRDS("proteinFeaturesFiltered_featureValsFilled.rda")
proteinFeaturesFiltered <- readRDS("proteinFeaturesFiltered.rda")

n_conditions <- length(unique(design_matrix$Condition))
for (i in seq(1,n_conditions,1)) {
  conditions <- unique(design_matrix$Condition)
  conditions <- conditions[-i]
  sub_design_matrix <- subset(design_matrix, Condition %in% conditions)
  sub_proteinFeaturesFiltered_featureValsFilled <- subset(proteinFeaturesFiltered_featureValsFilled, Condition %in% conditions)
  sub_pepTracesList_filtered <- pepTracesList_filtered[sub_design_matrix$Sample_name]
  class(sub_pepTracesList_filtered) <- "tracesList"

  #' ## Perform peptide-level differential expression testing for all features
  protein_DiffExprPep <- testDifferentialExpression(featureVals = sub_proteinFeaturesFiltered_featureValsFilled,
                                                    compare_between = "Condition",
                                                    level = "peptide",
                                                    measuredOnly = FALSE)

  saveRDS(protein_DiffExprPep, paste0("protein_DiffExprPep_",conditions[1],"_",conditions[2],".rda"))

  pdf(paste0("protein_DiffExprPep_stats_",conditions[1],"_",conditions[2],".pdf"), width = 4, height=4)
    hist(protein_DiffExprPep$pVal, breaks = 100)
    hist(protein_DiffExprPep$global_pVal, breaks = 100)
    hist(log(protein_DiffExprPep$qint1), breaks = 100)
    hist(log(protein_DiffExprPep$qint2), breaks = 100)
    hist(log(protein_DiffExprPep$global_int2_imp), breaks = 100)
    hist(log(protein_DiffExprPep$global_int1_imp), breaks = 100)
    hist(protein_DiffExprPep$log2FC, breaks = 50)
    hist(protein_DiffExprPep$global_log2FC, breaks = 50)
  dev.off()

  #' ## Aggregate differential expression results to the protein level
  protein_DiffExprProtein <- aggregatePeptideTests(protein_DiffExprPep)
  saveRDS(protein_DiffExprProtein, paste0("protein_DiffExprProtein_",conditions[1],"_",conditions[2],".rda"))

  pdf(paste0("protein_DiffExprProtein_stats_",conditions[1],"_",conditions[2],".pdf"), width = 4, height=4)
    hist(protein_DiffExprProtein$pVal, breaks = 100)
    hist(protein_DiffExprProtein$global_pVal, breaks = 100)
    hist(protein_DiffExprProtein$medianLog2FC, breaks = 50)
    hist(protein_DiffExprProtein$global_medianLog2FC, breaks = 50)
  dev.off()

  #' # Make volcano plots
  #' General volcano plots
  plotVolcano(protein_DiffExprPep, PDF = T, pBHadj_cutoff = 0.05, name = paste0("protein_DiffExprPep_",conditions[1],"_",conditions[2]))
  plotVolcano(protein_DiffExprProtein, PDF = T, pBHadj_cutoff = 0.05, name = paste0("protein_DiffExprProtein_",conditions[1],"_",conditions[2]))

  plotVolcano(protein_DiffExprPep, PDF = T, pBHadj_cutoff = 0.05, name = paste0("protein_DiffExprPep_",conditions[1],"_",conditions[2]), level = "global")
  plotVolcano(protein_DiffExprProtein, PDF = T, pBHadj_cutoff = 0.05, name = paste0("protein_DiffExprProtein_",conditions[1],"_",conditions[2]), level = "global")

  #' Get differential proteins
  allProteins <- unique(protein_DiffExprProtein$feature_id)
  diffProteins <- unique(protein_DiffExprProtein[pBHadj<0.05][abs(medianLog2FC)>1]$feature_id)
  write.table(allProteins, paste0("allProteins_",conditions[1],"_",conditions[2],".txt"), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(diffProteins, paste0("diffProteins_",conditions[1],"_",conditions[2],".txt"), col.names = F, row.names = F, quote = F, sep = "\t")
  diffProteins_global <- unique(protein_DiffExprProtein[global_pBHadj<0.05][abs(global_medianLog2FC)>1]$feature_id)
  write.table(diffProteins_global, paste0("diffProteins_global_",conditions[1],"_",conditions[2],".txt"), col.names = F, row.names = F, quote = F, sep = "\t")

  #' Volcanoplts highlighting different Proteins
  highlight_protein <- diffProteins[1]
  plotVolcano(protein_DiffExprProtein, highlight=highlight_protein, PDF = T, pBHadj_cutoff = 0.05, name = paste0("protein_DiffExprProtein_",highlight_protein,"_",conditions[1],"_",conditions[2]))
  plotVolcano(protein_DiffExprProtein, highlight=highlight_protein, PDF = T, pBHadj_cutoff = 0.05, name = paste0("protein_DiffExprProtein_",highlight_protein,"_",conditions[1],"_",conditions[2]), level = "global")

  #' # Plot some differential protein features
  pdf(paste0("diffProteins_",conditions[1],"_",conditions[2],"_top20.pdf"), width=8, height=7)
  for (id in diffProteins[1:20]) {
    plotFeatures(feature_table = proteinFeaturesFiltered,
                 traces = sub_pepTracesList_filtered,
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

  testFilled_norm <- sub_proteinFeaturesFiltered_featureValsFilled[feature_id %in% diffProteins]

  local_vs_global_test <- testLocalVsGlobal(testFilled_norm,name=paste0("local_vs_global_stats_",conditions[1],"_",conditions[2]))
  saveRDS(local_vs_global_test,paste0("local_vs_global_test_",conditions[1],"_",conditions[2],".rds"))

  local_vs_global_test_protein <- aggregateLocalVsGlobal(local_vs_global_test)
  saveRDS(local_vs_global_test_protein,paste0("local_vs_global_test_protein_",conditions[1],"_",conditions[2],".rds"))

  local_vs_global_test_sig <- unique(local_vs_global_test_protein[pBHadj<0.05][abs(medianDiff) > 0.3]$feature_id)
  write.table(local_vs_global_test_sig, paste0("local_vs_global_test_sig_",conditions[1],"_",conditions[2],".txt"), col.names = F, row.names = F, quote = F, sep = "\t")

  pdf(paste0("local_vs_global_test_sig_features_",conditions[1],"_",conditions[2],".pdf"), width=10,height=7)
  for (id in local_vs_global_test_sig) {
    plotFeatures(feature_table = proteinFeaturesFiltered,
                 traces = sub_pepTracesList_filtered,
                 feature_id = id,
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
