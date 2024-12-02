#### Chargement données et librairies ####
library(DESeq2)
library(ggplot2)
counts_data <- read.table("counts.txt", header = T)

## Mise en forme données ##
# Sélection des colonnes nécessaires
counts_matrix <- counts_data[, 7:ncol(counts_data)] 
# Utiliser 'Geneid' comme noms de lignes
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
results$geneID <- rownames(results)
results <- as.data.frame(results) #retirer les valeurs NA
results <- na.omit(results)

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
        panel.border = element_rect(color = "black", fill = NA, size = 1)) # Utilisation de 'size' pour la bordure


# Sauvegarder le volcano plot en PDF
ggsave("Volcano_plot.pdf", plot = volcano_plot, dpi = 300, width = 8, height = 6)


## MA plot traduction ##
library(KEGGREST)

#Les ids impliqués dans la traduction : 
# "RNA polymerase - Staphylococcus aureus subsp. aureus NCTC8325" 
#sao03010 
#"Ribosome - Staphylococcus aureus subsp. aureus NCTC8325" 
#sao00970 
#"Aminoacyl-tRNA biosynthesis - Staphylococcus aureus subsp. aureus NCTC8325" 
#sao03060 

#sélection des gènes impliqués dans la traduction
genesList03010 = gsub("sao:","gene-", keggLink("sao03010")[-1,2])
genesList00970 = gsub("sao:","gene-", keggLink("sao00970")[-1,2])
genesList03060 = gsub("sao:","gene-", keggLink("sao03060")[-1,2])

gene_translation = c(genesList03010,genesList00970, genesList03060)
gene_tRNA = genesList00970[-grep("T", genesList00970)]

## version 1 : en refaisant toute l'analyse juste pour ce set de gènes ##
counts_matrix_translation <- counts_matrix[rownames(counts_matrix) %in% gene_translation,] 

# Analyse
dds_trans <- DESeqDataSetFromMatrix(countData = counts_matrix_translation,
                                    colData = coldata,
                                    design = ~ condition)
dds_trans <- DESeq(dds_trans)
results_translation <- results(dds_trans)
results_translation$geneID <- rownames(results_translation)
results_translation <- na.omit(results_translation) #retirer les valeurs NA

results_translation$significant <- results_translation$padj < 0.05  # Identifier les gènes significatifs

# Plot
ggplot(results_translation, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 2, alpha = 0.9) +
  geom_point(data = results_translation[rownames(results_translation) %in% gene_tRNA, ],
             aes(x = log2(baseMean), y = log2FoldChange),
             shape = 21, size = 2, color = "black", fill = NA, stroke = 1.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey20", linewidth = 0.8) +
  scale_color_manual(values = c("grey35", "red"),
                     labels = c("Non-significant", "Significant"), 
                     name = "Gene Status") +
  labs(title = "Log2 Fold Change vs Mean Counts",
       x = "log2(Mean Counts)",
       y = "Log2 Fold Change",
       color = "Significant") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_rect(color = "grey20", fill = NA, linewidth = 1))

ggsave("MA_plot_translation.pdf", dpi = 300, width = 7.3, height = 6)


## version 2 : juste en récupérant les résultats de l'analyse globale
results_translation2 <- results[results$geneID %in% gene_translation,]

# plot
ggplot(results_translation2, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 2, alpha = 0.9) +
  geom_point(data = results_translation2[rownames(results_translation2) %in% gene_tRNA, ],
             aes(x = log2(baseMean), y = log2FoldChange),
             shape = 21, size = 2, color = "black", fill = NA, stroke = 1.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey20", linewidth = 0.8) +
  scale_color_manual(values = c("grey35", "red"),
                     labels = c("Non-significant", "Significant"), 
                     name = "Gene Status") +
  labs(title = "Log2 Fold Change vs Mean Counts",
       x = "log2(Mean Counts)",
       y = "Log2 Fold Change",
       color = "Significant") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.border = element_rect(color = "grey20", fill = NA, linewidth = 1))
ggsave("MA_plot_translation2.pdf", dpi = 300, width = 7.3, height = 6)
