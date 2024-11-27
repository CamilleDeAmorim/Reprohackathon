
library("KEGGREST")
########
# Notes 
########

# The strain ID used in the article is "sao," as referenced in their bibliography (reference 78).
# Retrieving the pathway identifiers involved in translation.
# Pathway IDs related to translation:  
# "RNA polymerase - Staphylococcus aureus subsp. aureus NCTC8325"  : sao03010  
# "Ribosome - Staphylococcus aureus subsp. aureus NCTC8325" : sao00970  
# "Aminoacyl-tRNA biosynthesis - Staphylococcus aureus subsp. aureus NCTC8325"  : sao03060 


# Suppress messages during KEGGREST calls
suppressMessages({
  genesList03010 <- keggLink("sao03010")
  genesList03060 <- keggLink("sao03060")
  genesList00970 <- keggLink("sao00970")
})

# Concatenate gene lists and remove the first row (organism identifier)
globListGenes <- rbind(
  genesList03010[-1, ],
  genesList00970[-1, ],
  genesList03060[-1, ]
)

# Reformat gene identifiers by replacing "sao:" with "gene-"
gene_ids <- gsub("sao:", "", globListGenes[, 2])

# Save the list of all translation-related genes
write.table(gene_ids, file = "assets/translation_genes.txt", row.names = F, col.names = F, quote = F)

# Extract aminoacyl-related genes and reformat identifiers
amino_gene_ids <- gsub("sao:", "", (genesList03060[-1, ])[, 2])

# Save the aminoacyl gene list
write.table(amino_gene_ids, file = "assets/AminoAcyl_tRNA_synthetases_genes.txt", row.names = F, col.names = F, quote = F)

######################
#   two output files:
######################
# translation_genes.txt : Contains all genes involved in translation,
# and Aminoacyl_tRNA_synthetases_genes.txt : Contains only the genes specifically
#involved in aminoacyl-tRNA biosynthesis.

