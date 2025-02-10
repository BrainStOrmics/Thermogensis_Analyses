# this file shows an example of the edgeR analysis
library(edgeR)


file <- "tables/GSE201135_rc_3.txt" 
data <- read.table(file, header=TRUE, sep="\t", row.names = 1)
head(data)


### get group information and contrast matrix
group_names <- unique(gsub("_[0-9]+$", "", colnames(data)))
group <- factor(sapply(colnames(data), function(x) {
  if (x != "GeneSymbol") {
    gsub("_[0-9]+$", "", x)
  }
}))
design <- model.matrix(~ 0 + group) 


contrast_list <- list()
for (i in 1:(length(group_names) - 1)) {
  for (j in (i + 1):length(group_names)) {
    contrast_name <- paste(group_names[i], "vs", group_names[j], sep = "_")
    contrast_expr <- paste("group", group_names[i], " - ", "group", group_names[j], sep = "")
    contrast_list[[contrast_name]] <- contrast_expr
  }
}
contrast_exprs <- paste(names(contrast_list), "=", unlist(contrast_list), collapse = ",\n  ")
contrast_matrix_code <- paste("contrast_matrix <- makeContrasts(\n",
                              contrast_exprs, ",\n\n  levels = design\n)", sep = "")
#cat(contrast_matrix_code)
eval(parse(text=contrast_matrix_code))


### DGE analysis

# edge R
data_matrix <- data[, sapply(data, is.numeric)]

dge <- DGEList(counts = data_matrix, group = group)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)
attributes(dge)
fit <- glmFit(dge, design)

### Summary and Output DGE analysis

# dirs
gse <- sub("_.*", "", basename(file))
output_dir <- "data/DEG/dges"
summary_dir <- "data/DEG/summary"

gse_output_dir <- file.path(output_dir, gse)
if (!dir.exists(gse_output_dir)) {
  dir.create(gse_output_dir, recursive = TRUE)
} else {
  cat("Directory for GSE ", gse, " already exists. Skipping creation.\n")
}

# summary table init
summary_table <- data.frame(
  GSE = character(0),
  N_of_Samples_Treatment = integer(0),
  N_of_Samples_Control = integer(0),
  Conditions_Treatment = character(0),
  Conditions_Control = character(0),
  Upregulated_Genes = integer(0),
  Downregulated_Genes = integer(0),
  tm_status = character(0),
  tm_logFC = numeric(0),
  tm_pvalue = numeric(0),
  stringsAsFactors = FALSE
)


contrast_list <- colnames(contrast_matrix)
results_list <- list()

# get DGEs tables
contrast_name = names(contrast_list)[1]
for (contrast_name in contrast_list) {
  cat("Processing contrast: ", contrast_name, "\n")

  fit2 <- glmLRT(fit, contrast = contrast_matrix[, contrast_name])
  #DGE <- topTags(fit2, n = nrow(data_matrix)) # n = Inf
  DGE <- topTags(fit2, n = Inf)
  
  DGE <- as.data.frame(DGE)
  
  conditions <- strsplit(contrast_name, "_vs_")
  conditions_treatment <- conditions[[1]][1]
  conditions_control <- conditions[[1]][2]
  
  N_of_Samples_Treatment <- sum(grepl(conditions_treatment, colnames(dge$counts)))
  N_of_Samples_Control <- sum(grepl(conditions_control, colnames(dge$counts)))
  
  # summary(DGE)
  # ggplot(DGE, aes(x=logFC, y=-log10(FDR))) + geom_point(shape=1) + theme_bw()
  # summary(log10(DGE$PValue))
  # summary(log10(DGE$FDR))
  
  DGE$status <- ifelse(DGE$FDR< 0.05 & DGE$logFC >= 0.2, "Up", 
                       ifelse(DGE$FDR < 0.05 & DGE$logFC <= -0.2, "Down", "None"))
  
  up_genes <- sum(DGE$status == "Up")
  down_genes <- sum(DGE$status == "Down")
  

  cat("Upregulated genes: ", up_genes, "\n")
  cat("Downregulated genes: ", down_genes, "\n")
  

  if ("Tmem41b" %in% rownames(DGE)) {
    tmem41b_status <- DGE["Tmem41b", "status"]
    tmem41b_logFC <- DGE["Tmem41b", "logFC"]
    tmem41b_pvalue <- DGE["Tmem41b", "FDR"]
    cat("Tmem41b status: ", tmem41b_status, "\n")
    cat("Tmem41b logFC: ", tmem41b_logFC, "\n")
    cat("Tmem41b P.Value: ", tmem41b_pvalue, "\n")
  } else {
    tmem41b_status <- NA
    tmem41b_logFC <- NA              
    tmem41b_pvalue <- NA         
  }
  
  results_list[[contrast_name]] <- DGE
  result_file_path <- file.path(gse_output_dir, paste0(contrast_name, "_DGE_results.txt"))
  DGE$Gene_Symbol <- rownames(DGE)
  write.table(DGE, result_file_path, sep = "\t", quote = FALSE, row.names = FALSE)

  summary_table <- rbind(summary_table, data.frame(
    GSE = gse, 
    N_of_Samples_Treatment = N_of_Samples_Treatment,
    N_of_Samples_Control = N_of_Samples_Control,
    Conditions_Treatment = conditions_treatment,
    Conditions_Control = conditions_control,
    Upregulated_Genes = up_genes,
    Downregulated_Genes = down_genes,
    tm_status = tmem41b_status,
    tm_logFC = tmem41b_logFC,
    tm_pvalue = tmem41b_pvalue,
    stringsAsFactors = FALSE
  ))
}

print(summary_table)

summary_file_path <- file.path(summary_dir, paste0(gse, "_summary_table.txt"))
write.table(summary_table, summary_file_path, sep = "\t", quote = FALSE, row.names = FALSE)

#### empty the folder if need
folder_path <- "data/DEG/dges/GSE151905"
unlink(file.path(folder_path, "*"), recursive = TRUE)
list.files(folder_path)
