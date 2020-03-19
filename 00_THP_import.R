source("/Users/isabell/Desktop/projects/THP_SEC_SWATH_MS/THP_CCprofiler_scripts/THP_header.R")

#' # Read TRIC output
data <-  fread("../CCprofiler_input/aligned_filtered.csv")
# change column names according to CCprofiler conventions
setnames(data,c("ProteinName","FullPeptideName","filename","Intensity"),c("protein_id","peptide_id","filename","intensity"))
# remove non-proteotypic peptides
data <- data[grep(";",data$protein_id, invert = TRUE)]
# remove peak group rank > 1
data <- subset(data, peak_group_rank == 1)
# reduce data to necessary columns
data_sub <- subset(data,select=c("protein_id","peptide_id","filename","intensity"))
# get uniprot IDs
data_sub[, protein_id := gsub(".*\\|(.*?)\\|.*", "\\1", protein_id)]

#' # Load annotation table
fraction_annotation_table <- unique(subset(data_sub, select="filename"))
fraction_annotation_table[,fraction_number := gsub('^marclaud_.*_.*_\\s*|\\s*\\.mzXML.*$', '', filename)]
fraction_annotation_table[,condition_id := gsub('^marclaud_.*_*_(\\w+)_.*.mzXML$', '\\1', filename)]
fraction_annotation_table[,replicate_id := gsub('^marclaud_.*_(\\w)_.*.mzXML$', '\\1', filename)]
fraction_annotation_table[,replicate_id := str_replace_all(replicate_id, c("A" = "1", "B" = "2", "C" = "3"))]
fraction_annotation_table[,sample := paste(c(condition_id, replicate_id), collapse = "_"), by=c("filename","fraction_number")]

fraction_annotation_table[,condition_id := ifelse(condition_id=="Blue", "stimulated", condition_id)]
fraction_annotation_table[,condition_id := ifelse(condition_id=="red", "differentiated", condition_id)]
fraction_annotation_table[,condition_id := ifelse(condition_id=="white", "undifferentiated", condition_id)]

fraction_annotation_table[,sample := gsub("Blue", "stimulated", sample)]
fraction_annotation_table[,sample := gsub("red", "differentiated", sample)]
fraction_annotation_table[,sample := gsub("white", "undifferentiated", sample)]

fraction_annotation_sub <- subset(fraction_annotation_table, select=c("filename","fraction_number","sample"))

#' # Find missing fractions and create mock fractions:
unique_samples <- unique(fraction_annotation_sub$sample)

get_missing_fractions <- function(x, max_frac){
  measured_fractions <- unique(x$fraction_number)
  all_fractions <- as.character(format(seq(1, max_frac, 1), digits = 2))
  all_fractions <- gsub(" ","0", all_fractions)
  missing_fractions <- all_fractions[which(! all_fractions %in% measured_fractions)]
  mock_fraction_annotation <- data.table(fraction_number=missing_fractions)
  mock_fraction_annotation[,filename:=paste0("mock_",unique(x$sample),"_",fraction_number)]
  mock_fraction_annotation[,sample:=unique(x$sample)]
}

for (s in unique_samples) {
  dt <- subset(fraction_annotation_sub, sample == s)
  max_fraction <- max(fraction_annotation_sub$fraction_number)
  fraction_annotation_sub <- data.table::rbindlist(list(fraction_annotation_sub,get_missing_fractions(dt, max_fraction)), use.names=TRUE)
}

fraction_annotation_sub[, fraction_number := as.integer(fraction_number)]

fraction_annotation_table <- fraction_annotation_sub[,condition_id:=gsub("_.*","",sample)]
fraction_annotation_table <- fraction_annotation_table[,replicate_id:=gsub(".*_","",sample)]
fraction_annotation_table[, fraction_number := as.integer(fraction_number)]

mock_data <- data.table(protein_id="P55011",
                        peptide_id="AAAAAAAAAAAAAAAGAGAGAK",
                        filename=fraction_annotation_table[grep("mock",filename)]$filename,
                        intensity=0)

data_sub <- rbind(data_sub,mock_data)

#' # Molecular weight calibration
molecularWeightCalibration <- fread("../CCprofiler_input/calibrationTable_Yarra.csv")
calibration = calibrateMW(molecularWeightCalibration,
                          PDF=T,
                          plot=TRUE)
saveRDS(calibration,"calibration.rds")

#' # Create design matrix
design_matrix <- unique(subset(fraction_annotation_table, select=c("condition_id","replicate_id","sample")))
setnames(design_matrix, c("condition_id","replicate_id","sample"), c("Condition","Replicate","Sample_name"))
design_matrix$Condition <- factor(design_matrix$Condition, levels=c("undifferentiated","differentiated","stimulated"))
saveRDS(design_matrix,"design_matrix.rds")

#' # Import traces list
samples <- unique(fraction_annotation_table$sample)
# Import data as traces object for each sample
traces_list <- lapply(samples,function(x){
  message(x)
  ann <- subset(fraction_annotation_sub, sample==x)
  ann <- subset(ann, filename %in% data_sub$filename)
  ann <- subset(ann, select = c("filename","fraction_number"))
  setkey(ann,fraction_number)
  data_sub_select <- subset(data_sub,filename %in% ann$filename)
  data_sub_select <- merge(data_sub_select,ann, by=c("filename"))
  setkey(data_sub_select,fraction_number)
  data_sub_select[,fraction_number:=NULL]
  traces <- importPCPdata(input_data=data_sub_select,fraction_annotation=ann)
  return(traces)
})
names(traces_list) = samples
class(traces_list) <- "tracesList"
# remove input data to clean storage
rm(data,data_sub)
gc()

#' # Annotate traces with information from uniprot
pepTraces_raw <- annotateTraces(traces=traces_list,
                            trace_annotation=exampleTraceAnnotation,
                            traces_id_column = "protein_id",
                            trace_annotation_id_column = "Entry",
                            trace_annotation_mass_column = "Mass",
                            uniprot_mass_format = TRUE,
                            replace_whitespace = TRUE)

#' # Annotate traces with molecular weight calibration
pepTraces_raw <- annotateMolecularWeight(pepTraces_raw,
                                     calibration)

summary(pepTraces_raw)
saveRDS(pepTraces_raw, "pepTraces_raw.rds")
