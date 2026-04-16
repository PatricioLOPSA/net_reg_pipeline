library(biomaRt)
library(dplyr)
library(data.table)

#Script for obtaining protein coding gene stable ids and gene symbols from biomart. Used for filtering count data to only include protein coding genes.

# Set up the connection to Ensembl BioMart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve the list of protein-coding genes
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                              filters = "biotype",
                              values = "protein_coding",
                              mart = ensembl)

# Save the list of protein-coding genes to a file
fwrite(protein_coding_genes, file = "protein_coding_genes.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
