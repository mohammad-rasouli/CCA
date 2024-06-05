library(biomaRt)
setwd ("C:/Users/mohammad/Desktop/Thesis/datasets")
transcript <- read.csv("RNA_filtered.csv", header = TRUE)
ensg <- colnames(transcript) # taking ensg ids and save them as a list

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") #creating biomart object
conversion <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                    filters = 'ensembl_gene_id',
                    values = ensg,
                    mart = ensembl)

#change and save the converted col names 
colnames(transcript) <- conversion$hgnc_symbol[match(colnames(transcript), conversion$ensembl_gene_id)]
write.csv(transcript, "transcriptome_genesymbol.csv", row.names = FALSE)
