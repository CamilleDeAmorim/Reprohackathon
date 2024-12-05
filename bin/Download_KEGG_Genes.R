library("KEGGREST")
########
# Notes 
########

# The strain ID used in the article is "sao," as referenced in their bibliography (reference 78).
# Retrieving the pathway identifiers involved in translation.
# Pathway IDs related to translation:   
# "Ribosome - Staphylococcus aureus subsp. aureus NCTC8325" : sao03010   
# "Aminoacyl-tRNA biosynthesis - Staphylococcus aureus subsp. aureus NCTC8325"  : sao00970
# "Translation factors - Staphylococcus aureus subsp. aureus NCTC8325" : sao03012 

#Download genes 
# + remove the first row (organism identifier) 
# + Reformat gene identifiers by replacing "sao:" with "gene-"
suppressMessages({
  genesList03010 = gsub("sao:","", keggLink("sao03010")[-1,2])
  genesList00970 = gsub("sao:","", keggLink("sao00970")[-1,2])
})

#impossible de récupérer sao03012 avec keggLink : besoin de passer par une autre méthode
brite_raw <- keggGet("br:sao03012")[[1]]
brite_data <- strsplit(brite_raw, "\n")

gene_lines <- grep("SAOUHSC_", unlist(brite_data), value = TRUE)
gene_ids <- regmatches(gene_lines, gregexpr("SAOUHSC_\\d+", gene_lines))
gene_ids <- unlist(gene_ids)
genesList03012 <- unique(c(GeneID = gene_ids))

# Concatenate gene lists 
translation_ids <- c(genesList03010,genesList00970, genesList03012)

# Save the list of all translation-related genes
write.table(translation_ids, file = "assets/translation_genes.txt", row.names = F, col.names = F, quote = F)

# Extract aminoacyl-related genes and reformat identifiers

tRNA_ids <- genesList00970[-grep("T", genesList00970)]

# Save the aminoacyl gene list
write.table(tRNA_ids, file = "assets/tRNA_synthetases_genes.txt", row.names = F, col.names = F, quote = F)

######################
#   two output files:
######################
# translation_genes.txt : Contains all genes involved in translation,
# and Aminoacyl_tRNA_synthetases_genes.txt : Contains only the genes specifically
#involved in aminoacyl-tRNA biosynthesis.

