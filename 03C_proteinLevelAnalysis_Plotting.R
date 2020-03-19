source("/Users/isabell/Desktop/projects/THP_SEC_SWATH_MS/THP_CCprofiler_scripts/THP_header.R")

#' # Load data
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")
design_matrix <- readRDS("design_matrix.rds")
calibration <- readRDS("calibration.rds")

setwd("proteinCentric")
proteinFeaturesFiltered_featureValsFilled <- readRDS("proteinFeaturesFiltered_featureValsFilled.rda")
proteinFeaturesFiltered <- readRDS("proteinFeaturesFiltered.rda")

diffProteins_UD <- fread("diffProteins_differentiated_undifferentiated.txt", sep="\t", header=F)
diffProteins_US <- fread("diffProteins_stimulated_undifferentiated.txt", sep="\t", header=F)
diffProteins_DS <- fread("diffProteins_stimulated_differentiated.txt", sep="\t", header=F)

diffProteins_combi <- unique(c(diffProteins_UD$V1,diffProteins_US$V1,diffProteins_DS$V1))

#' # Plot some differential protein features
n_plots <- ceiling(length(diffProteins_combi)/100)
for(r in seq(1,n_plots,1)){
  pdf(paste0("diffProteinTraces_combi_part",r,".pdf"), width=10,height=7)
  for (id in diffProteins_combi[((r-1)*100 + 1):min((r*100),length(diffProteins_combi))]) {
    plotFeatures(feature_table = proteinFeaturesFiltered,
                 traces = pepTracesList_filtered,
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
}




# local vs global

local_vs_global_diffProteins_UD <- fread("local_vs_global_test_sig_differentiated_undifferentiated.txt", sep="\t", header=F)
local_vs_global_diffProteins_US <- fread("local_vs_global_test_sig_stimulated_undifferentiated.txt", sep="\t", header=F)
local_vs_global_diffProteins_DS <- fread("local_vs_global_test_sig_stimulated_differentiated.txt", sep="\t", header=F)

local_vs_global_diffProteins_combi <- unique(c(local_vs_global_diffProteins_UD$V1,local_vs_global_diffProteins_US$V1,local_vs_global_diffProteins_DS$V1))

#' # Plot some differential protein features
n_plots <- ceiling(length(local_vs_global_diffProteins_combi)/100)
for(r in seq(1,n_plots,1)){
  pdf(paste0("local_vs_global_diffProteinTraces_combi_part",r,".pdf"), width=10,height=7)
  for (id in local_vs_global_diffProteins_combi[((r-1)*100 + 1):min((r*100),length(local_vs_global_diffProteins_combi))]) {
    plotFeatures(feature_table = proteinFeaturesFiltered,
                 traces = pepTracesList_filtered,
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
}

