.libPaths("/home/seb01ann/R/x86_64-pc-linux-gnu-library/4.3/")

library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(ggrepel)
library(coin)


# Metadata
mdata <- read.csv("/storage/projects/P024_ChIPseq_AE/csv/metadata.csv", header = TRUE)
metadata <- data.frame(Samples = mdata$Sample_name_x, Sex = mdata$Gender, Age = mdata$Age, Disease=mdata$Disease.Status, Type = mdata$cell_type_x)




tpm <- read.table("/storage/projects/P034_RNAseq_AE/nf-rnaseq/star_rsem/rsem.merged.gene_tpm.tsv", header=TRUE, row.names=1, sep="\t")
tpm <- tpm[, -1]
subset_tpm <- tpm[ , metadata$Samples]
head(subset_tpm)




variant_callers = read.csv("/storage/projects/P024_ChIPseq_AE/rstudio/mutect2_macs3_freebayes.csv", header=TRUE, sep=",")

# Filter variants as before
filtered_variants <- variant_callers %>%
  group_by(Hugo_Symbol) %>%
  filter(n() > 1) %>%
  arrange(Hugo_Symbol)



# T-test
# Initialize results data frame
variant_info <- data.frame(
  Hugo_Symbol = character(), 
  Chromosome = character(),
  Start_Position = numeric(),
  End_Position = numeric(),
  P_val = numeric(), 
  P_adj = numeric(), 
  log2FC = numeric(), 
  Perm_p_val = numeric(), 
  Perm_p_adj = numeric(), 
  Mutant = character(), 
  Non_mutant = character()
)

snp_colors <- brewer.pal(n = 6, name = "Accent")

# Skip problematic cases
valid_expr_data <- function(mut_expr, non_mut_expr) {
  if (length(mut_expr) < 3 || length(non_mut_expr) < 3) return(FALSE) # Sufficient size
  if (all(is.na(mut_expr)) || all(is.na(non_mut_expr))) return(FALSE) # Not missing
  if (sd(mut_expr, na.rm = TRUE) == 0 || sd(non_mut_expr, na.rm = TRUE) == 0) return(FALSE) # Has variability
  TRUE
}

pval_to_asterisks <- function(p) {
  if (is.na(p)) return("n.s.")             
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("n.s.")
}




# Wilcoxon rank-sum test (Wilcox test) 
### TPM that are often skewed, even after log transformation.
### A non-parametric test used to compare the distribution of two independent groups. 
# Permutation test
### A non-parametric method that assesses the significance of an observed effect by randomly shuffling labels to generate a null distribution of the test statistic.

# Function to perform differential expression analysis
expression_analysis <- function(dge_df, n_permutations = 10000) {
  mutant_expr <- dge_df$Log2TPM[dge_df$Group == "Mutant"]
  non_mutant_expr <- dge_df$Log2TPM[dge_df$Group == "Non-Mutant"]
  
  if (!valid_expr_data(mutant_expr, non_mutant_expr)) {
    return(list(p_val = NA, log2FC = NA, Perm_p_val = NA))
  }
  
  # Observed Wilcoxon test
  obs_test <- wilcox.test(mutant_expr, non_mutant_expr, paired = FALSE, exact = FALSE)
  p_val <- obs_test$p.value
  
  # Log2 fold change
  mean_log2_mutant <- mean(mutant_expr, na.rm = TRUE)
  mean_log2_non_mutant <- mean(non_mutant_expr, na.rm = TRUE)
  log2FC <- mean_log2_mutant - mean_log2_non_mutant
  
  # Permutation test using coin package
  dge_df$Group <- as.factor(dge_df$Group)
  perm_test <- wilcox_test(Log2TPM ~ Group, data = dge_df, distribution = approximate(nresample = n_permutations))
  perm_p_val <- pvalue(perm_test)
  
  return(list(p_val = p_val, log2FC = log2FC, Perm_p_val = perm_p_val))
}




# Store all results for proper multiple testing correction
all_results <- list()

# Main analysis loop
for (gene_of_interest in unique(filtered_variants$Hugo_Symbol)) {
  # Check if gene exists in expression data
  if (!gene_of_interest %in% rownames(subset_tpm)) {
    # warning(paste("Gene", gene_of_interest, "not found in expression data"))
    next
  }
  
  # Get the sample list of the gene of interest
  filtered_rows <- filtered_variants[filtered_variants$Hugo_Symbol == gene_of_interest, ]
  
  # Merge sample names from different callers
  filtered_rows$Sample_Names_merged <- apply(filtered_rows[, 
                                                           c("Sample_Names_mutect2", "Sample_Names_macs3", "Sample_Names_freebayes")], 1, function(x) {
                                                             valid_samples <- x[x != "" & !is.na(x)]
                                                             paste(unique(valid_samples), collapse = ", ")
                                                           })
  
  sample_list <- lapply(filtered_rows$Sample_Names_merged, function(x) {
    unlist(strsplit(x, ",\\s*"))
  })
  
  for (i in seq_along(sample_list)) {
    mut_samples <- unique(sample_list[[i]])
    all_samples <- colnames(subset_tpm)
    non_mut_samples <- setdiff(all_samples, mut_samples)
    
    # Get gene expression for this gene
    gene_expr <- suppressWarnings(as.numeric(subset_tpm[gene_of_interest, ]))
    names(gene_expr) <- colnames(subset_tpm)
    
    # Create analysis data frame
    plot_data <- data.frame(
      Gene = gene_of_interest,
      Sample = names(gene_expr),
      Expression = gene_expr,
      Group = ifelse(names(gene_expr) %in% mut_samples, "Mutant", "Non-Mutant"),
      stringsAsFactors = FALSE
    )
    
    plot_data$Log2TPM <- log2(plot_data$Expression + 1)
    plot_data <- plot_data[!is.na(plot_data$Log2TPM), ]
    
    # Perform differential expression analysis
    dge <- expression_analysis(plot_data, n_permutations = 10000)
    
    if (!is.na(dge$p_val) && !is.na(dge$Perm_p_val)) {
      # Store result for multiple testing correction
      result_id <- paste(gene_of_interest, i, sep = "_")
      all_results[[result_id]] <- list(
        gene = gene_of_interest,
        variant_idx = i,
        chr = filtered_rows[i, "Chromosome"],
        start = filtered_rows[i, "Start_Position"],
        end = filtered_rows[i, "End_Position"],
        p_val = dge$p_val,
        log2FC = dge$log2FC,
        perm_p_val = dge$Perm_p_val,
        mut_samples = mut_samples,
        non_mut_samples = non_mut_samples,
        plot_data = plot_data
      )
    }
  }
}



# False Discovery Rate (FDR) 
### Benjamini-Hochberg/FDR adjustment controls for multiple testing by limiting the expected proportion of false positives among significant results.
# Apply multiple testing correction to all results at once
if (length(all_results) > 0) {
  
  # Extract all p-values
  all_p_vals <- sapply(all_results, function(x) x$p_val)
  all_perm_p_vals <- sapply(all_results, function(x) x$perm_p_val)
  
  # Apply Benjamini-Hochberg correction to ALL p-values together
  p_adj <- p.adjust(all_p_vals, method = "BH")
  perm_p_adj <- p.adjust(all_perm_p_vals, method = "BH")
  
  # Store results and create plots for significant findings
  for (i in seq_along(all_results)) {
    result <- all_results[[i]]
    result$p_adj <- p_adj[i]
    result$perm_p_adj <- perm_p_adj[i]
    
    # Create plots for significant results
    if (result$p_val < 0.05 && result$perm_p_val < 0.05) {
      # Store in variant_info data frame
      variant_info <- rbind(variant_info, data.frame(
        Hugo_Symbol = result$gene,
        Chromosome = result$chr,
        Start_Position = result$start, 
        End_Position = result$end,
        P_val = result$p_val,
        P_adj = result$p_adj,         
        log2FC = result$log2FC,
        Perm_p_val = result$perm_p_val,
        Perm_p_adj = result$perm_p_adj,  
        Mutant = paste(result$mut_samples, collapse = ","),
        Non_mutant = paste(result$non_mut_samples, collapse = ","),
        stringsAsFactors = FALSE
      ))
      # Save results statistically significant genes to CSV
      write.csv(variant_info, "/storage/projects/P024_ChIPseq_AE/rstudio/variants_info.csv", row.names = FALSE)
      
      plot_data <- result$plot_data
      
      # Main expression boxplot
      plot_colors <- c("Mutant" = snp_colors[result$variant_idx %% length(snp_colors) + 1],  "Non-Mutant" = "gray")
      
      p1 <- ggplot(plot_data, aes(x = Group, y = Log2TPM, fill = Group)) +
        geom_violin(trim = FALSE, alpha = 0.7) +  # Violin plot
        geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +  # Boxplot inside violin
        geom_jitter(width = 0.2, alpha = 0.6) +    # Individual points
        labs(title = paste("Expression level -", result$gene), y = "Log2(TPM + 1)", x = "") +
        scale_fill_manual(values = plot_colors) +
        theme_minimal() +
        guides(fill = "none") +
        coord_cartesian(ylim = c(0, max(plot_data$Log2TPM, na.rm = TRUE) * 1.5)) +
        annotate("text", x = 1.5, y = max(plot_data$Log2TPM, na.rm = TRUE) * 1.2, size = 4, hjust = 0.5,
                 label = paste0("Wilcox test ",  pval_to_asterisks(result$p_val), 
                                "\nPermutation test ",  pval_to_asterisks(result$p_val)))
      
      print(p1)
      
    }
  }
}


merged_variants <- merge(variant_callers, variant_info, by=c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position"))

merged_variants <- merged_variants %>%
  mutate(Coordinates = paste0(Chromosome, ":", Start_Position, "-", End_Position), .after = Hugo_Symbol) %>%
  select(-Chromosome, -Start_Position, -End_Position, -Mutant, -Non_mutant, -Sample_Names_merged)

write.csv(merged_variants, "/storage/projects/P024_ChIPseq_AE/rstudio/merged_variants.csv", row.names = FALSE)
head(merged_variants)

