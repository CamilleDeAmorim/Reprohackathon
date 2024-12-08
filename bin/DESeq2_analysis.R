#!/usr/bin/env Rscript
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

# arguments by default 
default_args = c("results/06_Counting_results/counts.txt", # counting file
               "assets/translation_genes.txt", # file containing translation genes
               "assets/tRNA_synthetases_genes.txt", # file containing tRNA synthetases genes
               "assets/GeneSpecificInformation_NCTC8325.txt", # Correspondance file code - symbol of gene
               "results/07_Final_results/DESeq2_Results/") # results directory

# test if there are arguments ; if not by default arguments are used
if (length(args)<length(default_args)) {  
     args[length(args)+1:length(default_args)] <- default_args[length(args)+1:length(default_args)]
} 

counts_file = args[1]
gene_trans_f = args[2]
gene_tRNA_f = args[3]
gene_corres_f = args[4]
results_directory = args[length(args)] 

volcano_file = paste(results_directory,"Volcano_plot.png", sep = "")
MAplot_file = paste(results_directory,"MA_plot.png", sep = "")
MAplot_translation_file1 = paste(results_directory,"MA_plot_translation1.png", sep = "")
MAplot_translation_file2 = paste(results_directory,"MA_plot_translation2.png", sep = "")
if (!dir.exists(results_directory)) {
  dir.create(results_directory, recursive = TRUE) 
  cat("Directory created:", results_directory, "\n")
}

#### installation of data and required libraries ####
library(DESeq2)
library(ggplot2)
counts_data <- read.table(counts_file, skip = 1, header = T, sep = '\t')
counts_data$Geneid <- gsub("gene-","", counts_data$Geneid)
gene_translation <- read.table(gene_trans_f, header = F)[,1]
gene_tRNA <- read.table(gene_tRNA_f, header = F)[,1]
gene_corres <- read.table(gene_corres_f, skip = 1, blank.lines.skip=TRUE, nrows=4116)
colnames(gene_corres) <- c("Geneid", "GeneName")

a_afficher <- c("frr", "infA", "infB", "infC", "tsf", "pth")

## Data formatting ##
# Selecting necessary data
counts_matrix <- counts_data[, 7:ncol(counts_data)] 
# Use 'Geneid" for row names
counts_data$Geneid <- gsub("gene-","", counts_data$Geneid)
rownames(counts_matrix) <- counts_data$Geneid

#  Table of conditions
# The first 3 samples are persisters, the others are control samples 
conditions <- factor(c(rep("2", 3), rep("1", 3)))
coldata <- data.frame(row.names = colnames(counts_matrix), condition = conditions)

# Creation of DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = coldata,
                              design = ~ condition)

#### Analysis #### 
# Normalization and differential expression
dds <- DESeq(dds)
results <- results(dds)
results <- as.data.frame(results) 
results$Geneid <- row.names(results) 
results <- na.omit(results) # removing NA values

# Total number of genes differentially expressed, upregulated and downregulated
n_significant <- nrow(results[results$padj < 0.05 , ])
n_upregulated <- nrow(results[results$padj < 0.05 & results$log2FoldChange >0, ])
n_downregulated <- nrow(results[results$padj < 0.05 & results$log2FoldChange < 0, ])

cat("Number of genes differentially expressed :", n_significant, "\n")
cat("Number of upregulated genes :", n_upregulated, "\n")
cat("Number of downregulated genes :", n_downregulated, "\n")


#### Volcano plot ####

# Groupes based on significance and the expression 
results$group <- "Non-significant"
results$group[results$padj < 0.05 & results$log2FoldChange > 0] <- "Upregulated"
results$group[results$padj < 0.05 & results$log2FoldChange < 0] <- "Downregulated"

# Volcano plot creation
volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Non-significant" = "grey", 
                                "Upregulated" = "red", 
                                "Downregulated" = "blue")) +
  labs(title = "p-value ajustée en fonction du log2FC", 
       x = "Log2 Fold Change", 
       y = "-Log10(padj)", 
       color = "Gene Status") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, size = 1))

#save figure
ggsave(volcano_file, plot = volcano_plot, dpi = 300, width = 8, height = 6)


### MA plots ###
results$significant <- results$padj < 0.05  # identify significant genes

# Adjustment of logFC to constrain values within the range of -4.1 to 4.1
results$adjusted_log2FC <- pmin(pmax(results$log2FoldChange, -4.1), 4.1)
results$shape <- ifelse(results$log2FoldChange > 4.1, "triangle_up",
                                 ifelse(results$log2FoldChange < -4.1, "triangle_down", "circle"))

ggplot(results, aes(x = baseMean, y = adjusted_log2FC)) +
  geom_point(aes(color = significant, fill = significant, shape = shape), 
             size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.8) +
  scale_color_manual(values = c("black", "red"), 
                     labels = c("Non-significant", "Significant"), 
                     name = "Gene Status") + 
  scale_fill_manual(values = c("black", "red"), guide = "none") + #guide = none : to remove legend
  scale_shape_manual(values = c("circle" = 16, "triangle_up" = 24, "triangle_down" = 25), 
                     guide = "none") + 
  scale_x_log10() +
  coord_cartesian(ylim = c(-4.2, 4.2)) +
  labs(title = "Log2 Fold Change vs Mean Counts",
       x = "Mean Counts (log scale)",
       y = "Log2 Fold Change") + 
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, size = 1))
ggsave(MAplot_file, dpi = 300, width = 10, height = 6)


## MA plot for translation genes##

## Version 1 : DESeq2 analysis done on the new set of genes (translation only)  ##
counts_matrix_translation <- counts_matrix[rownames(counts_matrix) %in% gene_translation,] 
# Analysis
dds_trans <- DESeqDataSetFromMatrix(countData = counts_matrix_translation,
                                    colData = coldata,
                                    design = ~ condition)
dds_trans <- DESeq(dds_trans)
results_translation <- results(dds_trans)
results_translation <- as.data.frame(results_translation)
results_translation$Geneid <- row.names(results_translation)
results_translation <- na.omit(results_translation) # removing NA values


results_translation$significant <- results_translation$padj < 0.05  # identifying significant genes
results_translation <- merge(results_translation, gene_corres, by = "Geneid", all.x=TRUE, all.y=F) # Adding gene names

# Plot
ggplot(results_translation, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 2, alpha = 0.9) +
  # Adding Black Outlines for tRNA Genes
  geom_point(data = results_translation[results_translation$Geneid %in% gene_tRNA, ],
             aes(x = log2(baseMean), y = log2FoldChange),
             shape = 21, size = 2, color = "black", fill = NA, stroke = 1.2) +
  #Adding Segments Between Genes and Their Names
  geom_segment(data = results_translation[results_translation$GeneName %in% a_afficher, ],
               aes(x = log2(baseMean), y = log2FoldChange, 
                   xend = log2(baseMean) + 0.2, yend = log2FoldChange + 0.2),
               color = "black", size = 0.7) +
  # Adding gene names
  geom_text(data = results_translation[results_translation$GeneName %in% a_afficher, ],
            aes(label = GeneName),
            size = 4, vjust = -0.5, hjust = 0, nudge_x = 0.2, nudge_y = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey20", size = 0.8) +
  scale_color_manual(values = c("grey35", "red"),
                     labels = c("Non-significant", "Significant"), 
                     name = "Gene Status") +
  labs(title = "Log2 Fold Change vs Mean Counts\n DESeq2 analysis done on translation genes",
       x = "log2(Mean Counts)",
       y = "Log2 Fold Change",
       color = "Significant") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_rect(color = "grey20", fill = NA, size = 1))

ggsave(MAplot_translation_file1, dpi = 300, width = 7.3, height = 6)


## version 2 : DESeq2 analysis done on all dataset, we retrieve results only for translation genes
results_translation2 <- results[results$Geneid %in% gene_translation,]
results_translation2 <- merge(results_translation2, gene_corres, by = "Geneid", all.x=TRUE, all.y=F) #Adding gene names

# plot
ggplot(results_translation2, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 2, alpha = 0.9) +
  # Adding Black Outlines for tRNA Genes
  geom_point(data = results_translation2[results_translation2$Geneid %in% gene_tRNA, ],
             aes(x = log2(baseMean), y = log2FoldChange),
             shape = 21, size = 2, color = "black", fill = NA, stroke = 1.2) +
  # Adding Segments Between Genes and Their Names
  geom_segment(data = results_translation2[results_translation2$GeneName %in% a_afficher, ],
               aes(x = log2(baseMean), y = log2FoldChange, 
                   xend = log2(baseMean) + 0.2, yend = log2FoldChange + 0.2),
               color = "black", size = 0.7) +
  # Adding gene names
  geom_text(data = results_translation2[results_translation2$GeneName %in% a_afficher, ],
            aes(label = GeneName),
            size = 4, vjust = -0.5, hjust = 0, nudge_x = 0.2, nudge_y = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey20", size = 0.8) +
  scale_color_manual(values = c("grey35", "red"),
                     labels = c("Non-significant", "Significant"), 
                     name = "Gene Status") +
  labs(title = "Log2 Fold Change vs Mean Counts\n DESeq2 analysis done on all genes",
       x = "log2(Mean Counts)",
       y = "Log2 Fold Change",
       color = "Significant") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_rect(color = "grey20", fill = NA, size = 1))

ggsave(MAplot_translation_file2, dpi = 300, width = 7.3, height = 6)


####################################
# Generating Histograms for p-values
####################################


# p-values histogram for all the genes

pvalue_histogram <- ggplot(results, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of p-values for all the dataset",
       x = "p-value",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "none")
# save the figure
histogram_file <- paste(results_directory, "/pvalues_histogram_all_genes.png", sep = "")
ggsave(histogram_file, plot = pvalue_histogram, dpi = 300, width = 8, height = 6)

# p-values histogram for translation genes

pvalue_histogram_2 <- ggplot(results_translation, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of p-values for translation genes",
       x = "p-value",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "none")
# save the figure
histogram_file_2 <- paste(results_directory, "/pvalues_histogram_translation_genes.png", sep = "")
ggsave(histogram_file_2, plot = pvalue_histogram_2, dpi = 300, width = 8, height = 6)


