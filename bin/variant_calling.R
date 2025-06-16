#!/usr/bin/env Rscript

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript variant_calling.R <input_maf_directory> <sample_id> <caller_name>")
}
input_maf_directory <- args[1]
sample_id <- args[2]
caller_name <- args[3]

# Libraries
suppressMessages(library(maftools))

# Define variant classes of interest
interest <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
  "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", 
  "Missense_Mutation", "Silent", "3'Flank", "3'UTR", "5'Flank", "5'UTR", 
  "IGR", "Intron", "RNA", "Splice_Region", "Targeted_Region"
)
# Get all MAF files
maf_files <- list.files(input_maf_directory, pattern = "\\.maf$", full.names = TRUE, recursive = TRUE)
cat("Found", length(maf_files), "MAF files\n")

if (length(maf_files) == 0) {
  stop("No .maf.gz files found in directory: ", input_maf_directory)
}

# Load MAFs
maf_list <- lapply(maf_files, read.maf, vc_nonSyn = FALSE)
names(maf_list) <- basename(dirname(maf_files))

# Merge MAFs
aml <- merge_mafs(maf_list, vc_nonSyn = FALSE)

save_and_show_plot <- function(plot_func, filename, ...) {
  png(filename, width = 20, height = 16, units = "in", pointsize = 30, res = 300)
  plot_func(...)
  dev.off()
  plot_func(...)  # optional if you want to show it in interactive mode
}

# Output plots
out_prefix <- paste0(sample_id, "_", caller_name)
save_and_show_plot(plot_func = oncoplot,
                   filename = paste0(out_prefix, "_waterfallplot.png"),
                   maf = aml, top = 15, legendFontSize = 0.8,
                   fontSize = 0.6, gene_mar = 8, showTumorSampleBarcodes = TRUE)

save_and_show_plot(plot_func = plotmafSummary,
                   filename = paste0(out_prefix, "_mafSummary.png"),
                   maf = aml, rmOutlier = TRUE, addStat = 'mean', dashboard = TRUE, titvRaw = TRUE)

save_and_show_plot(plot_func = tcgaCompare,
                   filename = paste0(out_prefix, "_tcgaCompare.png"),
                   maf = aml, cohortName = "AML")
