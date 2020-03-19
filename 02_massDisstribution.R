source("/Users/isabell/Desktop/projects/THP_SEC_SWATH_MS/THP_CCprofiler_scripts/THP_header.R")

#' Load data
design_matrix <- readRDS("design_matrix.rds")
protTraces <- readRDS("protein_traces_list.rds")
pepTracesList_filtered <- readRDS("pepTracesList_filtered.rds")

#' Set output directory
outdir = "globalAssemblyChanges"
if (! dir.exists(outdir)) {
  dir.create(outdir)
}
setwd(outdir)

# evaluate protein mass distribution
protTraces_assembly <- annotateMassDistribution(protTraces)
saveRDS(protTraces_assembly,"protTraces_assembly.rds")

uniqueProteins <- Reduce(intersect, lapply(protTraces_assembly, function(x){x$trace_annotation$id}))
protTraces_assembly <- subset(protTraces_assembly, trace_subset_ids = uniqueProteins)

lapply(names(protTraces_assembly), function(x){
  summarizeMassDistribution(protTraces_assembly[[x]], 
                            PDF=T,
                            name=paste0("massDistribution_",x))})

n_conditions <- length(unique(design_matrix$Condition))
for (i in seq(1,n_conditions,1)) {
  conditions <- unique(design_matrix$Condition)
  conditions <- conditions[-i]
  sub_design_matrix <- subset(design_matrix, Condition %in% conditions)
  sub_protTraces_assembly <- protTraces_assembly[sub_design_matrix$Sample_name]
  class(sub_protTraces_assembly) <- "tracesList"
  sub_pepTracesList_filtered <- pepTracesList_filtered[sub_design_matrix$Sample_name]
  class(sub_pepTracesList_filtered) <- "tracesList"

  diffAssemblyState <- getMassAssemblyChange(sub_protTraces_assembly, sub_design_matrix, 
                                             plot=T, PDF=T,
                                             name = paste0("beta_pvalue_histogram_",conditions[1],"_",conditions[2]))
  saveRDS(diffAssemblyState, paste0("diffAssemblyState_",conditions[1],"_",conditions[2],".rds"))

  meanDiff_cutoff=0.3

  write.table(diffAssemblyState,paste0("diffAssemblyState_",conditions[1],"_",conditions[2],".txt"),sep="\t",quote=F,row.names = F,col.names=T)

  if (length(unique(design_matrix$Replicate)) > 1) {
    all_proteins_independentOfAssemblyState <- unique(diffAssemblyState$protein_id)
    proteins_noAssemblyChanges <- subset(diffAssemblyState, (abs(meanDiff) <= meanDiff_cutoff) & (betaPval_BHadj < 0.05))$protein_id
    more_assembled_proteins <- subset(diffAssemblyState, (meanDiff < -meanDiff_cutoff) & (betaPval_BHadj < 0.05))$protein_id
    less_assembled_proteins <- subset(diffAssemblyState, (meanDiff > meanDiff_cutoff) & (betaPval_BHadj < 0.05))$protein_id
    proteins_withAssemblyChanges <- subset(diffAssemblyState, (abs(meanDiff) > meanDiff_cutoff) & (betaPval_BHadj < 0.05))$protein_id
  } else {
    all_proteins_independentOfAssemblyState <- unique(diffAssemblyState$protein_id)
    proteins_noAssemblyChanges <- subset(diffAssemblyState, (abs(meanDiff) <= meanDiff_cutoff))$protein_id
    more_assembled_proteins <- subset(diffAssemblyState, (meanDiff < -meanDiff_cutoff))$protein_id
    less_assembled_proteins <- subset(diffAssemblyState, (meanDiff > meanDiff_cutoff))$protein_id
    proteins_withAssemblyChanges <- subset(diffAssemblyState, (abs(meanDiff) > meanDiff_cutoff))$protein_id

  }

  pdf(paste0("diffAssemblyState_AMF_AMF_",conditions[1],"_",conditions[2],".pdf"), width=6, height=4)
    p <- ggplot(diffAssemblyState, aes(x=meanAMF1, y=meanAMF2, colour=-log10(betaPval_BHadj))) +
    geom_point() +
    geom_abline(intercept = meanDiff_cutoff, slope = 1) +
    geom_abline(intercept = -meanDiff_cutoff, slope = 1) +
    theme_classic()
    print(p)
  dev.off()

  write.table(all_proteins_independentOfAssemblyState,paste0("all_proteins_independentOfAssemblyState_",conditions[1],"_",conditions[2],".txt"),sep="\t",quote=F,row.names = F,col.names=F)
  write.table(proteins_noAssemblyChanges,paste0("proteins_noAssemblyChanges_",conditions[1],"_",conditions[2],".txt"),sep="\t",quote=F,row.names = F,col.names=F)
  write.table(more_assembled_proteins,paste0("more_assembled_proteins_",conditions[1],"_",conditions[2],".txt"),sep="\t",quote=F,row.names = F,col.names=F)
  write.table(less_assembled_proteins,paste0("less_assembled_proteins_",conditions[1],"_",conditions[2],".txt"),sep="\t",quote=F,row.names = F,col.names=F)
  write.table(proteins_withAssemblyChanges,paste0("proteins_withAssemblyChanges_",conditions[1],"_",conditions[2],".txt"),sep="\t",quote=F,row.names = F,col.names=F)

  #' # Plot proteins_withAssemblyChanges
  examples <- proteins_withAssemblyChanges
  pdf(paste0("proteins_withAssemblyChanges_pepTraces_",conditions[1],"_",conditions[2],".pdf"), width=6, height=5)
  for(test_proteins in examples){
    pepTest <- subset(sub_pepTracesList_filtered, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
    plot(pepTest,
      design_matrix = sub_design_matrix,
      log = FALSE, monomer_MW=T, legend = F,
      PDF=F, name = paste0("PeptideTraces_",test_proteins),
      aggregateReplicates=T)
  }
  dev.off()
  pdf(paste0("proteins_withAssemblyChanges_protTraces_",conditions[1],"_",conditions[2],".pdf"), width=5, height=4)
  for(test_proteins in examples){
    protTest <- subset(sub_protTraces_assembly, trace_subset_ids = test_proteins)
    plot(protTest,
      design_matrix = sub_design_matrix,
      log = FALSE,
      monomer_MW=T,
      legend = T,
      PDF=F,
      name = paste0("ProteinTraces_",test_proteins),
      aggregateReplicates=T,
      collapse_conditions=T)
  }
  dev.off()
}

conditions <- unique(design_matrix$Condition)
compare1 <- fread(paste0("proteins_withAssemblyChanges_",conditions[1],"_",conditions[2],".txt"), header=F)
compare2 <- fread(paste0("proteins_withAssemblyChanges_",conditions[1],"_",conditions[3],".txt"), header=F)$V1
compare3 <- fread(paste0("proteins_withAssemblyChanges_",conditions[2],"_",conditions[3],".txt"), header=F)$V1

all_diff_proteins <- unique(c(compare3,compare2))
write.table(all_diff_proteins,paste0("proteins_withAssemblyChanges_allConditions.txt"),sep="\t",quote=F,row.names = F,col.names=F)

pdf(paste0("proteins_withAssemblyChanges_protTraces_allConditions.pdf"), width=5, height=4)
for(test_proteins in all_diff_proteins){
  protTest <- subset(protTraces_assembly, trace_subset_ids = test_proteins)
  plot(protTest,
    design_matrix = design_matrix,
    log = FALSE,
    monomer_MW=T,
    legend = T,
    PDF=F,
    name = paste0("ProteinTraces_",test_proteins),
    aggregateReplicates=T,
    collapse_conditions=T)
}
dev.off()

pdf(paste0("proteins_withAssemblyChanges_protTraces_allConditions_sep.pdf"), width=5, height=5)
for(test_proteins in all_diff_proteins){
  protTest <- subset(protTraces_assembly, trace_subset_ids = test_proteins)
  plot(protTest,
    design_matrix = design_matrix,
    log = FALSE,
    monomer_MW=T,
    legend = T,
    PDF=F,
    name = paste0("ProteinTraces_",test_proteins),
    aggregateReplicates=T,
    collapse_conditions=F)
}
dev.off()
