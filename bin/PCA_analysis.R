#!/usr/bin/env Rscript
# Retrieve command line arguments or use default values
args <- commandArgs(trailingOnly = TRUE)
default_args <- c("results/06_Counting_results/counts.txt", 
                  "results/07_Final_results/Supplementary_results",
                  "results/07_Final_results/DESeq2_Results/upregulated_translation_genes.txt",
                  "results/07_Final_results/DESeq2_Results/downregulated_translation_genes.txt")

if (length(args) < 4) {  
    args[length(args) + 1:4] <- default_args[length(args) + 1:4]
} 

counts_file <- args[1]
results_directory <- args[2]
upregulated_genes=args[3]
downregulated_genes=args[4]
if (!dir.exists(results_directory)) {
  dir.create(results_directory, recursive = TRUE) 
  cat("Directory created:", results_directory, "\n")
}
# Define the output PDF path
output_pdf <- paste0(results_directory, "/Supplementary_analysis.pdf")

# Load required libraries
suppressMessages(library(factoextra))
suppressMessages(library(FactoMineR))
suppressMessages(library(ggplot2))

# 1. Load the data
if (!file.exists(counts_file)) {
    stop(paste("File not found:", counts_file))
}
countData <- read.table(counts_file, header = TRUE, row.names = 1, skip = 1)

# Extract only numeric columns (raw counts)
countDataNumeric <- countData[, 6:ncol(countData)]

# Remove constant columns (zero variance)
constant_cols <- apply(countDataNumeric, 2, var) == 0
if (any(constant_cols)) {
    cat("Removing constant columns:", colnames(countDataNumeric)[constant_cols], "\n")
    countDataNumeric <- countDataNumeric[, !constant_cols]
}

# Remove constant rows (zero variance across all columns)
constant_rows <- apply(countDataNumeric, 1, var) == 0
if (any(constant_rows)) {
    countDataNumeric <- countDataNumeric[!constant_rows, ]
}

# Adjust sample names
sample_names <- gsub(".*Mapping_results\\.|\\.bam", "", colnames(countDataNumeric))
colnames(countDataNumeric) <- sample_names

# Define sample groups
groups <- c(rep("Persisters", 3), rep("Controls", 3))  # Assign groups

pdf(output_pdf, width = 10, height = 8)
###################################
#PCA Plot with ggplot2 for samples
###################################
transposed_matrix <- t(countDataNumeric)  # Transpose the data (samples as rows)

pca <- PCA(transposed_matrix, graph = FALSE)  # Use PCA from FactoMineR
pca_data <- as.data.frame(pca$ind$coord)  # Extract PCA coordinates
pca_data$Group <- groups                  # Add group information
pca_data$Sample <- rownames(pca_data)     # Add sample names
# Ensure Group is a clean factor with correct levels
pca_data$Group <- factor(pca_data$Group, levels = c("Persisters", "Controls"))

# Create the PCA plot using ggplot2
print(ggplot(pca_data, aes(x = Dim.1, y = Dim.2, color = Group)) +
        geom_point(size = 4) +
        geom_text(aes(label = Sample), vjust = -1, hjust = 0.5, size = 3) +
        scale_color_manual(
          values = c("Persisters" = "red", "Controls" = "blue"),
          labels = c("Persisters", "Controls")  # Explicit labels
        ) +
        labs(
          title = "PCA Plot for samples",
          x = paste0("Principal Component 1 (", round(pca$eig[1, 2], 1), "%)"),
          y = paste0("Principal Component 2 (", round(pca$eig[2, 2], 1), "%)"),
          color = NULL  # Suppress legend title
        ) +
        theme_minimal() +
        theme(
          legend.position = "top",  
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          plot.title = element_text(size = 14, hjust = 0.5)
        )) 
# Extract PCA coordinates for the samples
sample_coords <- pca$ind$coord

# Calculate contributions for PC1 and PC2
contrib_pc1 <- (sample_coords[, "Dim.1"]^2) / sum(sample_coords[, "Dim.1"]^2) * 100
contrib_pc2 <- (sample_coords[, "Dim.2"]^2) / sum(sample_coords[, "Dim.2"]^2) * 100

# Create a data frame for contributions
sample_contrib <- data.frame(
  Sample = rownames(sample_coords),
  Contribution_PC1 = contrib_pc1,
  Contribution_PC2 = contrib_pc2,
  Group = pca_data$Group  # Add group information (Persisters/Controls)
)
# Plot contributions for PC1 with percentages above bars
ggplot(sample_contrib, aes(x = reorder(Sample, -Contribution_PC1), y = Contribution_PC1, fill = Group)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = round(Contribution_PC1, 1)),  # Add text with rounded percentages
            vjust = -0.5, size = 3) +  # Position and size of text
  scale_fill_manual(values = c("Persisters" = "red", "Controls" = "black")) +
  labs(
    title = "Sample Contributions to PC1",
    x = "Samples",
    y = "Contribution (%)",
    fill = "Sample Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )
# Plot contributions for PC2 with percentages above bars
ggplot(sample_contrib, aes(x = reorder(Sample, -Contribution_PC2), y = Contribution_PC2, fill = Group)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = round(Contribution_PC2, 1)),  # Add text with rounded percentages
            vjust = -0.5, size = 3) +  # Position and size of text
  scale_fill_manual(values = c("Persisters" = "red", "Controls" = "black")) +
  labs(
    title = "Sample Contributions to PC2",
    x = "Samples",
    y = "Contribution (%)",
    fill = "Sample Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )



#########################
#Scree Plot with ggplot2
#########################
eig <- pca$eig[, 2]  # Percentage of explained variance (eigenvalues)
eig_data <- data.frame(
    Axes = seq_along(eig),  # Principal component indices
    Inertia = eig           # Percentage of explained variance
)

print(ggplot(eig_data, aes(x = Axes, y = Inertia)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.6) +  # Add bars
    geom_text(aes(label = round(Inertia, 1)), vjust = -0.5, size = 3.5) +  # Add labels above bars
    ylim(0,80 ) +  # Adjust Y-axis limits
    labs(
        title = "Scree Plot",
        x = "Principal Components",
        y = "Percentage of Variance"
    ) +
    theme_minimal() +
    theme(
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
    ))
#########################################################
#PCA Plot for upregulated/downregulated translation genes
###########################################################

# Load the lists of genes
upregulated_genes <- read.table(upregulated_genes, stringsAsFactors = FALSE)$V1
downregulated_genes <- read.table(downregulated_genes, stringsAsFactors = FALSE)$V1

# Annotate genes with "Upregulated", "Downregulated", or "Neutral"
countDataNumeric$group <- "Neutral"  # Default to "Neutral"
rownames(countDataNumeric) <- gsub("^gene-", "", rownames(countDataNumeric))

countDataNumeric$group[rownames(countDataNumeric) %in% upregulated_genes] <- "Upregulated"
countDataNumeric$group[rownames(countDataNumeric) %in% downregulated_genes] <- "Downregulated"

# Filter only "Upregulated" and "Downregulated" genes
filtered_data <- countDataNumeric[countDataNumeric$group %in% c("Upregulated", "Downregulated"), ]

# Remove the "group" column for PCA
data_for_pca <- filtered_data[, !colnames(filtered_data) %in% "group"]

# Perform PCA
pca_genes <- prcomp(scale(data_for_pca), center = TRUE)
dim(data_for_pca)  # Dimensions de la matrice
head(rownames(data_for_pca))  # Noms des lignes (gènes ou échantillons ?)
# Calculate the percentage of explained variance for PC1 and PC2
explained_variance <- round(summary(pca_genes)$importance[2, 1:2] * 100, 1)  # Variance for PC1 and PC2
x_label <- paste0("Principal Component 1 (", explained_variance[1], "%)")
y_label <- paste0("Principal Component 2 (", explained_variance[2], "%)")


# Extract PCA coordinates
pca_coords <- as.data.frame(pca_genes$x)
pca_coords$Gene <- rownames(filtered_data)  # Assign gene names
pca_coords$group <- filtered_data$group     # Add "group" column

# Create the PCA plot
ggplot(pca_coords, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "grey")) +
  labs(
    title = "PCA of Up- and Down-regulated Genes",
    x = x_label,
    y = y_label,
    color = "Gene Status"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )


# extract scores for PC1
scores <- pca_genes$x[, "PC1"]

# contrib to genes for PC1
contributions <- (scores^2) / sum(scores^2) * 100  # Contributions en pourcentage

# create df for contributions
contrib_df <- data.frame(
  Gene = rownames(filtered_data),  # Les gènes sont dans les lignes de filtered_data
  Contribution = contributions,
  group = filtered_data$group      # Ajouter l'annotation de groupe
)


contrib_df <- contrib_df[order(-contrib_df$Contribution), ]

# Top 20 contributors
top_contrib <- contrib_df[1:20, ]

# Plot contributions to PC1
ggplot(top_contrib, aes(x = reorder(Gene, -Contribution), y = Contribution, fill = group)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8,width=0.6) +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "black")) +
  labs(
    title = "Top 20 Contributors to PC1",
    x = "Genes",
    y = "Contribution (%)",
    fill = "Gene Status"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )
# Extract scores for PC2
scores_pc2 <- pca_genes$x[, "PC2"]

# Calculate contributions of genes to PC2
contributions_pc2 <- (scores_pc2^2) / sum(scores_pc2^2) * 100  # Contributions in percentage

# Create a data frame for PC2 contributions
contrib_df_pc2 <- data.frame(
  Gene = rownames(filtered_data),  # Genes are in the rows of filtered_data
  Contribution = contributions_pc2,
  group = filtered_data$group      # Add group annotation (Upregulated/Downregulated)
)

# Sort contributions in descending order for PC2
contrib_df_pc2 <- contrib_df_pc2[order(-contrib_df_pc2$Contribution), ]

# Optional: Select the top 20 contributors for PC2
top_contrib_pc2 <- contrib_df_pc2[1:20, ]



ggplot(top_contrib_pc2, aes(x = reorder(Gene, -Contribution), y = Contribution, fill = group)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.6) +  # Create bar plot
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "gray")) +  # Custom colors
  labs(
    title = "Top 20 Contributors to PC2",  # Title of the plot
    x = "Genes",  # X-axis label
    y = "Contribution (%)",  # Y-axis label
    fill = "Gene Status"  # Legend title
  ) +
  theme_minimal() +  # Minimal theme for clean aesthetics
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5, size = 14),  # Center-align the title
    axis.title = element_text(size = 12),  # Set axis title size
    axis.text = element_text(size = 10),  # Set axis text size
    legend.position = "top"  # Move legend to the top
  )

# Close the PDF device
dev.off()

