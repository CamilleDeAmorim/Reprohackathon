#!/usr/bin/env Rscript
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

#arguments par défaut 
default_args = c("results/06_Counting_results/counts.txt", #fichier de comptage
               "assets/translation_genes.txt", #fichier avec gènes liés traduction
               "assets/tRNA_synthetases_genes.txt", # fichier avec gènes lié tRNA
               "assets/GeneSpecificInformation_NCTC8325.txt", #fichier de correspondance des noms de gènes
               "results/07_Final_results/") # dossier pour les résultats

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

#### Chargement données et librairies ####
library(DESeq2)
library(ggplot2)
counts_data <- read.table(counts_file, skip = 1, header = T, sep = '\t')
counts_data$Geneid <- gsub("gene-","", counts_data$Geneid)
gene_translation <- read.table(gene_trans_f, header = F)[,1]
gene_tRNA <- read.table(gene_tRNA_f, header = F)[,1]
gene_corres <- read.table(gene_corres_f, skip = 1, blank.lines.skip=TRUE, nrows=4116)
colnames(gene_corres) <- c("Geneid", "GeneName")

a_afficher <- c("frr", "infA", "infB", "infC", "tsf", "pth")

## Mise en forme données ##
# Sélection des colonnes nécessaires
counts_matrix <- counts_data[, 7:ncol(counts_data)] 
# Utiliser 'Geneid' comme noms de lignes
counts_data$Geneid <- gsub("gene-","", counts_data$Geneid)
rownames(counts_matrix) <- counts_data$Geneid

#  Tableau des conditions
# les trois premiers échantillons sont la condition test, les trois suivants le contrôle
conditions <- factor(c(rep("2", 3), rep("1", 3)))
coldata <- data.frame(row.names = colnames(counts_matrix), condition = conditions)

#création d'un objet DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = coldata,
                              design = ~ condition)

#### Analyse #### 
# Normalisation et analyse différentielle
dds <- DESeq(dds)
results <- results(dds)
results <- as.data.frame(results) 
results$Geneid <- row.names(results) 
results <- na.omit(results) #retirer les valeurs NA

# Nombre total de gènes significatifs, sur et sous exprimés
n_significant <- nrow(results[results$padj < 0.05 , ])
n_upregulated <- nrow(results[results$padj < 0.05 & results$log2FoldChange >0, ])
n_downregulated <- nrow(results[results$padj < 0.05 & results$log2FoldChange < 0, ])

cat("Nombre de gènes avec expression significativement différente :", n_significant, "\n")
cat("Nombre de gènes sur-exprimés :", n_upregulated, "\n")
cat("Nombre de gènes sous-exprimés :", n_downregulated, "\n")


#### Volcano plot ####

# Groupes basés sur la significativité et la direction de l'expression
results$group <- "Non-significant"
results$group[results$padj < 0.05 & results$log2FoldChange > 0] <- "Upregulated"
results$group[results$padj < 0.05 & results$log2FoldChange < 0] <- "Downregulated"

# Création du volcano plot
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

#sauvgarde données
ggsave(volcano_file, plot = volcano_plot, dpi = 300, width = 8, height = 6)


### MA plots ###
results$significant <- results$padj < 0.05  # Identifier les gènes significatifs

# Ajustement des logFC pour ne pas avoir de valeurs au dessus de 4,1 ou en dessous de -4,1
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
  scale_fill_manual(values = c("black", "red"), guide = "none") + #guide = none : pour enlever de la légende
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


## MA plot traduction ##

## version 1 : en refaisant toute l'analyse juste pour ce set de gènes ##
counts_matrix_translation <- counts_matrix[rownames(counts_matrix) %in% gene_translation,] 
# Analyse
dds_trans <- DESeqDataSetFromMatrix(countData = counts_matrix_translation,
                                    colData = coldata,
                                    design = ~ condition)
dds_trans <- DESeq(dds_trans)
results_translation <- results(dds_trans)
results_translation <- as.data.frame(results_translation)
results_translation$Geneid <- row.names(results_translation)
results_translation <- na.omit(results_translation) #retirer les valeurs NA


results_translation$significant <- results_translation$padj < 0.05  # Identifier les gènes significatifs
results_translation <- merge(results_translation, gene_corres, by = "Geneid", all.x=TRUE, all.y=F) #ajouter les noms des gènes

# Plot
ggplot(results_translation, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 2, alpha = 0.9) +
  # Ajouter les contours noirs pour les gènes tRNA
  geom_point(data = results_translation[results_translation$Geneid %in% gene_tRNA, ],
             aes(x = log2(baseMean), y = log2FoldChange),
             shape = 21, size = 2, color = "black", fill = NA, stroke = 1.2) +
  # Ajouter les segments entre les gènes et les noms
  geom_segment(data = results_translation[results_translation$GeneName %in% a_afficher, ],
               aes(x = log2(baseMean), y = log2FoldChange, 
                   xend = log2(baseMean) + 0.2, yend = log2FoldChange + 0.2),
               color = "black", size = 0.7) +
  # Ajouter les noms des gènes
  geom_text(data = results_translation[results_translation$GeneName %in% a_afficher, ],
            aes(label = GeneName),
            size = 4, vjust = -0.5, hjust = 0, nudge_x = 0.2, nudge_y = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey20", size = 0.8) +
  scale_color_manual(values = c("grey35", "red"),
                     labels = c("Non-significant", "Significant"), 
                     name = "Gene Status") +
  labs(title = "Log2 Fold Change vs Mean Counts\n analyse done on all genes",
       x = "log2(Mean Counts)",
       y = "Log2 Fold Change",
       color = "Significant") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_rect(color = "grey20", fill = NA, size = 1))

ggsave(MAplot_translation_file1, dpi = 300, width = 7.3, height = 6)


## version 2 : juste en récupérant les résultats de l'analyse globale
results_translation2 <- results[results$Geneid %in% gene_translation,]
results_translation2 <- merge(results_translation2, gene_corres, by = "Geneid", all.x=TRUE, all.y=F) #ajouter les noms des gènes

# plot
ggplot(results_translation2, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 2, alpha = 0.9) +
  # Ajouter les contours noirs pour les gènes tRNA
  geom_point(data = results_translation2[results_translation2$Geneid %in% gene_tRNA, ],
             aes(x = log2(baseMean), y = log2FoldChange),
             shape = 21, size = 2, color = "black", fill = NA, stroke = 1.2) +
  # Ajouter les segments entre les gènes et les noms
  geom_segment(data = results_translation2[results_translation2$GeneName %in% a_afficher, ],
               aes(x = log2(baseMean), y = log2FoldChange, 
                   xend = log2(baseMean) + 0.2, yend = log2FoldChange + 0.2),
               color = "black", size = 0.7) +
  # Ajouter les noms des gènes
  geom_text(data = results_translation2[results_translation2$GeneName %in% a_afficher, ],
            aes(label = GeneName),
            size = 4, vjust = -0.5, hjust = 0, nudge_x = 0.2, nudge_y = 0.2, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey20", size = 0.8) +
  scale_color_manual(values = c("grey35", "red"),
                     labels = c("Non-significant", "Significant"), 
                     name = "Gene Status") +
  labs(title = "Log2 Fold Change vs Mean Counts\n analyse done on all genes",
       x = "log2(Mean Counts)",
       y = "Log2 Fold Change",
       color = "Significant") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_rect(color = "grey20", fill = NA, size = 1))

ggsave(MAplot_translation_file2, dpi = 300, width = 7.3, height = 6)


#### Histogramme des p-valeurs ####

# Générer l'histogramme des p-valeurs
pvalue_histogram <- ggplot(results, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution des p-valeurs",
       x = "p-valeur",
       y = "Fréquence") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "none")

# Sauvegarder l'histogramme en PNG
histogram_file <- paste(results_directory, "/Pvalue_histogram.png", sep = "")
ggsave(histogram_file, plot = pvalue_histogram, dpi = 300, width = 8, height = 6)


## Récupérer les gènes sur/sous exprimés pour l'analyse ACP
############################################################
# Filtrer les gènes up-regulated
upregulated_genes <- subset(results_translation2, group == "Upregulated")

# Filtrer les gènes down-regulated
downregulated_genes <- subset(results_translation2, group == "Downregulated")

# Définir les chemins des fichiers de sortie
output_upregulated <- "results/07_Final_results/DESeq2_Results/upregulated_translation_genes.txt"
output_downregulated <- "results/07_Final_results/DESeq2_Results/downregulated_translation_genes.txt"

# Sauvegarder les listes dans des fichiers texte
write.table(upregulated_genes$geneID, file = output_upregulated, 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(downregulated_genes$geneID, file = output_downregulated, 
            row.names = FALSE, col.names = FALSE, quote = FALSE)



