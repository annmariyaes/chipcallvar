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
interest <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "IGR", "Intron", "RNA",
              "Silent", "Splice_Region", "Splice_Site", "Targeted_Region",
              "Translation_Start_Site", "Unknown")

# Load and merge MAF files
maf_files <- list.files(input_maf_directory, pattern = ".vep.maf", full.names = TRUE, recursive = TRUE)
print(length(maf_files))
maf_list <- lapply(maf_files, read.maf, vc_nonSyn = interest)
names(maf_list) <- basename(dirname(maf_files))
aml <- merge_mafs(maf_list, vc_nonSyn = interest)


# Output plots
save_and_show_plot <- function(plot_func, filename, ...) {
  png(filename, width = 20, height = 16, units = "in", pointsize = 30, res = 300)
  plot_func(...)
  dev.off()
  plot_func(...)  # optional if you want to show it in interactive mode
}

out_prefix <- paste0(sample_id, ".", caller_name)
save_and_show_plot(plot_func = oncoplot,
                   filename = paste0(out_prefix, ".waterfallplot.png"),
                   maf = aml, top = 15, legendFontSize = 0.8,
                   fontSize = 0.6, gene_mar = 8, showTumorSampleBarcodes = TRUE)

save_and_show_plot(plot_func = plotmafSummary,
                   filename = paste0(out_prefix, ".mafSummary.png"),
                   maf = aml, rmOutlier = TRUE, addStat = 'mean', dashboard = TRUE, titvRaw = TRUE)
