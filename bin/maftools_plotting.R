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

aml = read.maf(maf = input_maf_directory,vc_nonSyn = FALSE)

# Output plots
save_and_show_plot <- function(plot_func, filename, ...) {
  png(filename, width = 20, height = 16, units = "in", pointsize = 30, res = 300)
  plot_func(...)
  dev.off()
  plot_func(...)  # optional if you want to show it in interactive mode
}

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
