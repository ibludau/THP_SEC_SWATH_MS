source("/Users/isabell/Desktop/projects/THP_SEC_SWATH_MS/THP_CCprofiler_scripts/THP_header.R")

#' # Load data
pepTraces_raw <- readRDS("pepTraces_raw.rds")
design_matrix <- readRDS("design_matrix.rds")

#' # Subset raw traces to fractions 1-49 (missing data after this)
pepTraces_raw <- subset(pepTraces_raw, fraction_ids=seq(1,49,1))
saveRDS(pepTraces_raw, "pepTraces_raw.rds")

#' # QC
#' # Assess how well SEC traces align across all samples
alignTraces(pepTraces_raw, min_lag = -3, max_lag = 3, plot = T, PDF=T)
#' # Assess total intensity as a proxi for extraction efficiency across samples
plotGlobalIntensities(pepTraces_raw, plot = T, PDF=T)

#' Find missing values
#' (defined as having identifications in left and right neigbouring fractions):
pepTracesMV <- findMissingValues(pepTraces_raw,
                                 bound_left = 1,
                                 bound_right = 1,
                                 consider_borders = TRUE)

#' Impute NA values by fitting a spline:
pepTracesImp <- imputeMissingVals(pepTracesMV, method = "spline")

#' Plot imputation summary:
plotImputationSummary(pepTracesMV, pepTracesImp, PDF = T,
                      plot_traces = T, max_n_traces = 2)

saveRDS(pepTracesImp,"pepTracesImp.rds")

#' Normalize intensity values across samples
pepTracesNormalized <- normalizeByCyclicLoess(pepTracesImp, window = 3, step = 1, plot = TRUE, PDF = TRUE, name = "normalizeByCyclicLoess")
saveRDS(pepTracesNormalized,"pepTracesNormalized.rds")

#' # Assess total intensity after normalization
plotGlobalIntensities(pepTracesNormalized, plot = T, PDF=T, name = "IntensitySummary_postNormalization")

#' # Filter by consecutive IDs and sibling peptide correlation
# Combine all traces for filtering:
traces_combined <- combineTracesMutiCond(pepTracesNormalized)
# Filter by consecutive IDs:
pepTracesConsIds <- filterConsecutiveIdStretches(traces_combined,
                                                 min_stretch_length = 3,
                                                 remove_empty = T)

# Filter by maximum correlation for outlier removal:
pepTracesMaxCorr <- filterByMaxCorr(pepTracesConsIds,
                                cutoff = 0.5,
                                plot = T, PDF = T)

# Filter by consecutive SPC:
pepTracesSPC <- filterBySibPepCorr(pepTracesMaxCorr,
                                   absolute_spcCutoff = 0.2,
                                   plot = T, PDF = T)

#' # Subset pepTracesNormalized to valid peptides selected 
#' by consecutive filtering and SPC cutoff
validPeps <- unique(pepTracesSPC$trace_annotation$id)
pepTracesList_filtered <- lapply(pepTracesNormalized, function(x){
  subset(x, trace_subset_ids=validPeps)
})
class(pepTracesList_filtered) <- "tracesList"

#' Update traces with additional metrics for each fraction:
pepTracesList_filtered <- updateTraces(pepTracesList_filtered)

#' Inspect traces list:
summary(pepTracesList_filtered)

saveRDS(pepTracesList_filtered, "pepTracesList_filtered.rds")

#' # Plot some examplary peptide traces
examples <- unique(pepTracesList_filtered[[1]]$trace_annotation$protein_id)[1:20]
pdf("examplePeptideProfiles.pdf")
for(test_proteins in examples){
  pepTest <- subset(pepTracesList_filtered, trace_subset_ids = test_proteins, trace_subset_type = "protein_id")
  plot(pepTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF=F, name = paste0("PeptideTraces_",test_proteins))
}
dev.off()

#' # Global stats
traces_names <- c(
  "pepTraces_raw",
  "pepTracesList_filtered"
)

countDT <- data.table(evidence=character(), sample=character(), count=integer(), traces=character())

for (traces in traces_names) {
  #traces <- readRDS(paste0(traces_name,".rda"))
  traces_name <- traces
  traces <- eval(as.name(traces))

  getIDs <- lapply(traces, function(x){
    peptides <- unique(x$trace_annotation[grep("DECOY",protein_id, invert = T)]$id)
    proteins <- unique(x$trace_annotation[grep("DECOY",protein_id, invert = T)]$protein_id)
    list(peptides=length(peptides),proteins=length(proteins))
  })

  count_dt <- as.data.table(getIDs)
  count_dt[, evidence := c("peptides", "proteins")]
  count_dt <- melt(count_dt, id.vars = "evidence", variable.name = "sample", value.name = "count")
  count_dt[,traces:=traces_name]
  count_dt[, count := unlist(count)]
  countDT <- rbind(countDT,count_dt)
}

traces_names_dt <- data.table(
  traces=c("pepTraces_raw",
           "pepTracesList_filtered"),
  traces_name=c("raw",
                "filtered")
)

countDT <- merge(countDT, traces_names_dt, by="traces")

countDT$traces_name <- factor(countDT$traces_name, levels=traces_names_dt$traces_name)

write.table(countDT, "traces_stat_counts.txt", sep="\t",quote = F, col.names = T, row.names = F)

pdf("traces_stat_counts.pdf",height=9,width=12)
g <- ggplot(countDT,aes(x=traces_name,y=count,fill=evidence, group=sample)) +
  geom_bar(stat="identity") +
  facet_wrap(sample ~ evidence, scales="free",nrow=3) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_classic() +
  theme(legend.position="bottom") + theme(legend.title = element_blank()) +
  labs(y = "count") +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(g)
dev.off()

#' ## Quantify on protein level
protein_traces_list <- proteinQuantification(pepTracesList_filtered ,
                                             quantLevel="protein_id",
                                             topN = 2,
                                             keep_less = TRUE,
                                             rm_decoys = TRUE,
                                             use_sibPepCorr = FALSE,
                                             use_repPepCorr = FALSE,
                                             full_intersect_only = FALSE,
                                             verbose = FALSE)

#' ## Update fraction annotation for protein traces
protein_traces_list <- updateTraces(protein_traces_list)
saveRDS(protein_traces_list,"protein_traces_list.rds")

#' # Plot some examplary protein traces
pdf("exampleProteinProfiles.pdf")
for(test_proteins in examples){
  protTest <- subset(protein_traces_list, trace_subset_ids = test_proteins)
  plot(protTest, design_matrix = design_matrix, log = FALSE, legend = F, PDF=F, name = paste0("PeptideTraces_",test_proteins))
}
dev.off()
