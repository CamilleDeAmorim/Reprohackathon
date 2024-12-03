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

#Download genes 
# + remove the first row (organism identifier) 
# + Reformat gene identifiers by replacing "sao:" with "gene-"
suppressMessages({
  genesList03010 = gsub("sao:","gene-", keggLink("sao03010")[-1,2])
  genesList00970 = gsub("sao:","gene-", keggLink("sao00970")[-1,2])
  genesList03060 = gsub("sao:","gene-", keggLink("sao03060")[-1,2])
})

# Concatenate gene lists 
translation_ids <- c(genesList03010,genesList00970, genesList03060)

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

