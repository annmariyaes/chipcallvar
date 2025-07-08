#!/usr/bin/env Rscript

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript variant_calling.R <input_maf_directory> <cohort_id> <caller_name>")
}
input_maf_directory <- args[1]
cohort_id <- args[2]
caller_name <- args[3]

.libPaths("/storage/home/seb01ann/R/x86_64-pc-linux-gnu-library/4.3/")

# Libraries
suppressMessages(library(maftools))
suppressMessages(library(GenomicRanges))
suppressMessages(library(IRanges))
suppressMessages(library(dplyr))

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

# Load enhancer and promoter BED files and convert to GRanges
enhancers <- read.table("/storage/projects/P024_ChIPseq_AE/beds/enhancer_H3K27Ac_peaks.bed")
promoters <- read.table("/storage/projects/P024_ChIPseq_AE/beds/promoter_H3K27Ac_peaks.bed")
he <- read.table("/storage/projects/P024_ChIPseq_AE/beds/patient_enhancers-concat-sort-merge.bed")

colnames(enhancers) <- colnames(he) <- colnames(promoters) <- c("chr", "start", "end")

enhancer_gr <- GRanges(seqnames = enhancers$chr, ranges = IRanges(start = enhancers$start, end = enhancers$end))
promoter_gr <- GRanges(seqnames = promoters$chr, ranges = IRanges(start = promoters$start, end = promoters$end))
he_gr <- GRanges(seqnames = he$chr, ranges = IRanges(start = he$start, end = he$end))

# Create GRanges for variants
aml_gr <- GRanges(seqnames = aml@data$Chromosome,
                  ranges = IRanges(start = aml@data$Start_Position, 
                                   end = aml@data$End_Position))

# Identify overlaps
enhancer_hits <- findOverlaps(aml_gr, enhancer_gr)
promoter_hits <- findOverlaps(aml_gr, promoter_gr)
he_hits <- findOverlaps(aml_gr, he_gr)

# Extract overlapping variants
enhancer_vars <- aml@data[queryHits(enhancer_hits), ]
promoter_vars <- aml@data[queryHits(promoter_hits), ]
he_vars <- aml@data[queryHits(he_hits), ]

enhancer_vars$Region_Type <- "Enhancer"
promoter_vars$Region_Type <- "Promoter"
he_vars$Region_Type <- "Hyperactive Enhancer"

# Combine variants and remove unknown genes
combined_vars <- rbind(promoter_vars, enhancer_vars, he_vars) %>% 
  filter(Hugo_Symbol != "Unknown")

# Apply additional filters
cat("Before filtering:", nrow(combined_vars), "\n")

filtered_variants <- combined_vars %>%
  filter(!grepl("pseudogene", BIOTYPE, ignore.case = TRUE)) %>%   # Remove pseudogenes
  
  # include VAF in the data frame
  mutate(t_vaf = t_alt_count/(t_ref_count+t_alt_count), .after = t_alt_count) %>%
  
  # Remove variants with AF>=0.02, gnomADe_AF>=0.0001 DP<10, #ALT<5 and VAF<10%
  filter(
    (AF <= 0.02 | is.na(AF)) &                        # Keep low AF variants (likely somatic) or missing values
      (gnomADe_AF <= 0.0001 | is.na(gnomADe_AF)) &      # Rare or novel in gnomAD (filter out germline)
      t_depth >= 10 &                                    # Sufficient tumor read depth
      t_vaf >= 0.10 &                                     # Minimizes strand bias
      t_alt_count >= 5 & t_ref_count >= 5 &  
      FILTER != 'common_variant'                        # Remove known common variants
  ) %>% 
  
  arrange(Chromosome, Start_Position) 

cat("After filtering:", nrow(filtered_variants), "\n")

# save filtered variants
filename <- paste0("/storage/projects/P024_ChIPseq_AE/rstudio/", caller_name, "_filtered_variants.csv")
write.csv(filtered_variants, file = filename, row.names = FALSE)

# Output plots
save_and_show_plot <- function(plot_func, filename, ...) {
  png(filename, width = 20, height = 16, units = "in", pointsize = 30, res = 300)
  plot_func(...)
  dev.off()
  plot_func(...)  # optional if you want to show it in interactive mode
}

out_prefix <- paste0(cohort_id, ".", caller_name)
save_and_show_plot(plot_func = oncoplot,
                   filename = paste0(out_prefix, ".waterfallplot.png"),
                   maf = aml, top = 15, legendFontSize = 0.8,
                   fontSize = 0.6, gene_mar = 8, showTumorSampleBarcodes = TRUE)

save_and_show_plot(plot_func = plotmafSummary,
                   filename = paste0(out_prefix, ".mafSummary.png"),
                   maf = aml, rmOutlier = TRUE, addStat = 'mean', dashboard = TRUE, titvRaw = TRUE)

# Find recurrent mutations across samples
recurrent_mutations <- filtered_variants %>%
  # group_by(across(-Tumor_Sample_Barcode)) %>%
  group_by(Hugo_Symbol, Chromosome, Start_Position, End_Position, Strand, Reference_Allele, Allele, Variant_Classification, Variant_Type, Consequence, dbSNP_RS, AF, AFR_AF, AMR_AF, ASN_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF, gnomADe_AF, gnomADe_AFR_AF, gnomADe_AMR_AF, gnomADe_ASJ_AF, gnomADe_EAS_AF, gnomADe_FIN_AF, gnomADe_NFE_AF, gnomADe_OTH_AF, gnomADe_SAS_AF, FILTER, SOMATIC, Region_Type) %>%
  
  summarize(
    Num_Samples = n_distinct(Tumor_Sample_Barcode),
    Sample_Names = paste(unique(Tumor_Sample_Barcode), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(Num_Samples > 1) %>%
  arrange(desc(Num_Samples), desc(Chromosome), desc(Start_Position))

filename2 <- paste0("/storage/projects/P024_ChIPseq_AE/rstudio/", caller_name, "_recurrent_mutations.csv")
write.csv(recurrent_mutations, filename2, row.names = FALSE)
