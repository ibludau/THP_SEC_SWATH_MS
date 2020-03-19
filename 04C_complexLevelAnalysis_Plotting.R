source("/Users/isabell/Desktop/projects/THP_SEC_SWATH_MS/THP_CCprofiler_scripts/THP_header.R")

#' # Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
protein_traces_list <- readRDS("protein_traces_list.rds")
design_matrix <- readRDS("design_matrix.rds")
calibration <- readRDS("calibration.rds")

setwd("complexCentric")
combiTargetsPlusDecoys <- readRDS("combiTargetsPlusDecoys.rds")
complex_featureValsFilled <- readRDS("complex_featureValsFilled.rda")
scoredDataAll <- readRDS("complexFeaturesCollapsed.rda")

all_complexes <- unique(complex_featureValsFilled$complex_id)
pdf("allComplexes_entryName_combi.pdf", width=8, height=7)
for (id in all_complexes) {
  plotFeatures(feature_table = scoredDataAll,
               traces = protein_traces_list,
               feature_id = id,
               annotation_label = "Entry_name",
               design_matrix=design_matrix,
               calibration=calibration,
               peak_area = F,
               legend = F,
               onlyBest = F,
               monomer_MW=T)
}
dev.off()


#' each pairwise comparisons!
diffComplexes_UD <- fread("diffComplexes_differentiated_undifferentiated.txt", sep="\t", header=F)
diffComplexes_US <- fread("diffComplexes_stimulated_undifferentiated.txt", sep="\t", header=F)
diffComplexes_DS <- fread("diffComplexes_stimulated_differentiated.txt", sep="\t", header=F)

diffComplexes_combi <- unique(c(diffComplexes_UD$V1,diffComplexes_US$V1,diffComplexes_DS$V1))


pdf("diffComplexes_combi.pdf", width=8, height=7)
for (id in diffComplexes_combi) {
  plotFeatures(feature_table = scoredDataAll,
               traces = protein_traces_list,
               feature_id = id,
               design_matrix=design_matrix,
               calibration=calibration,
               peak_area = T,
               legend = F,
               onlyBest = F,
               monomer_MW=T,
              aggregateReplicates=T)
}
dev.off()

pdf("diffComplexes_entryName_combi.pdf", width=8, height=7)
for (id in diffComplexes_combi) {
  plotFeatures(feature_table = scoredDataAll,
               traces = protein_traces_list,
               feature_id = id,
               annotation_label = "Entry_name",
               design_matrix=design_matrix,
               calibration=calibration,
               peak_area = T,
               legend = F,
               onlyBest = F,
               monomer_MW=T,
              aggregateReplicates=T)
}
dev.off()

