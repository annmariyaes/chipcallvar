.libPaths("/home/seb01ann/R/x86_64-pc-linux-gnu-library/4.3/")
library(maftools)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


# Define variant classes of interest
interest <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "IGR", "Intron", "RNA", "Silent", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site", "Unknown")

# Define colors for variant classification
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(interest))
vc_cols <- setNames(my_colors, interest)
titv_cols <- brewer.pal(n = 8, name = 'Accent')
names(titv_cols) <- c('C>G', 'C>T', 'T>C', 'C>A', 'T>A', 'T>G')

# Load and merge MAF files
maf_files <- list.files("/storage/projects/P024_ChIPseq_AE/nf-macs3-merged/mafs_annotated", pattern = ".maf.gz", full.names = TRUE, recursive = TRUE)
print(length(maf_files))
maf_list <- lapply(maf_files, read.maf, vc_nonSyn = interest)
names(maf_list) <- basename(dirname(maf_files))
aml <- merge_mafs(maf_list, vc_nonSyn = interest)

# Define a helper function to save and show plots
save_and_show_plot <- function(plot_func, filename, ...) {
  png(filename, width = 20, height = 16, units = "in", pointsize = 30, res = 300)
  plot_func(...)
  dev.off()
  plot_func(...)
}


save_and_show_plot(plot_func = oncoplot, 
                   filename = "waterfallplot.png", 
                   maf = aml, top = 15, colors = vc_cols, draw_titv = TRUE, titv_col = titv_cols, legendFontSize=0.8, fontSize=0.6, gene_mar=8, showTumorSampleBarcodes=TRUE)

save_and_show_plot(plot_func = plotmafSummary,
                   filename = "mafSummary.png", 
                   maf = aml, rmOutlier = TRUE, addStat = 'mean', dashboard = TRUE, titvRaw = TRUE, color = vc_cols, titvColor = titv_cols)

save_and_show_plot(plot_func = tcgaCompare,
                   filename = "tcgaCompare.png", 
                   maf = aml, cohortName = "AML")

