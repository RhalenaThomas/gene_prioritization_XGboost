library(data.table)
library(dplyr)
library(biomaRt)

locus <- fread("../locus/gwas_locus_hg19_genes.tab")
locus <- locus[,c("symbol","locus")]

data <- fread("D30.DA.qtl_results_all.sorted_gwas_snp_after_fdr.txt")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "GRCh37")
gene <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = c('ensembl_gene_id', 'biotype'), values = list(data$feature_id, biotype="protein_coding"), mart = ensembl)

colnames(gene) <- c("feature_id", "symbol")

final <- merge(data, gene)

result <- final %>% group_by(symbol) %>% summarize(Cuomo_MEAN_EFFECT = mean(beta), Cuomo_MAX_EFFECT = max(beta))

final22 <- merge(result, locus)

write.table(final22, "D30.DA.qtl_results_all.sorted_gwas_snp_summarized.txt", quote = F, row.names = F, sep = "\t")

final2_subset.list <- lapply(sort(unique(final22$locus)), function(i){
                final2 <- subset(final22, locus == i)
                colnames(final2) <- paste0(colnames(final2),"_NBH")
                if(nrow(final2) == 1){
                        final3 <- final2[,2:(ncol(final2)-1)]
                        final3[final3 > 0] <- 1
                        return(data.frame(symbol = final2$symbol,final3, locus = final2$locus))
                }
                final3 <- apply(final2[,2:(ncol(final2)-1)], 2, function(i){
                        abs(i)/max(abs(i))
                })
                return(data.frame(symbol = final2$symbol,final3, locus = final2$locus))
})

is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))


final3 <- Reduce(rbind, final2_subset.list)
#final3[is.nan(final3)] <- 0

write.table(final3, "D30.DA.qtl_results_all.sorted_gwas_snp_summarized_nbh.txt", quote = F, row.names = F, sep = "\t")
