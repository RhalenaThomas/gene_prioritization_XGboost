library(data.table)
library(biomaRt)

ci <- readLines("ci_gwas_genes.txt")
locus <- fread("../locus/gwas_locus_hg19_genes.tab")
locus <- locus[,c("symbol","locus")]

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "GRCh37")
gene <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = c('ensembl_gene_id', 'biotype'), values = list(ci, biotype="protein_coding"), mart = ensembl)

colnames(gene)[2] <- "symbol"


write.table(merge(subset(gene, symbol > 0),locus), "ci_gwas_genes_hgnc.txt", sep = "\t", quote = F, row.names = F)


