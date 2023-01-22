library(data.table)
library(dplyr)

locus <- fread("../locus/gwas_locus_hg19_genes.tab")
locus <- locus[,c("symbol","locus")]

cis_eqtl <- fread("BrainMeta_cis_eqtl_summary/gwas_snp_concensus_BrainMeta_cis_eQTL.txt")
cis_sqtl <- fread("BrainMeta_cis_sqtl_summary/gwas_snp_concensus_BrainMeta_cis_sQTL.txt")
cage <- fread("cage_eqtl_data/gwas_snp_concensus_cage_eqtl_data.txt")

cis_eqtl_result <- cis_eqtl %>% group_by(Gene) %>% summarize(BrainMeta_cis_eqtl_MEAN_EFFECT = mean(b), BrainMeta_cis_eqtl_MAX_EFFECT = max(b))
cis_sqtl_result <- cis_sqtl %>% group_by(Gene) %>% summarize(BrainMeta_cis_sqtl_MEAN_EFFECT = mean(b), BrainMeta_cis_sqtl_MAX_EFFECT = max(b))
cage_result <- cage %>% group_by(Gene) %>% summarize(cage_eqtl_MEAN_EFFECT = mean(b), cage_eqtl_MAX_EFFECT = max(b))


result <- Reduce(function(x, y) merge(x, y, by = "Gene", all = T), list(cis_eqtl_result,cis_sqtl_result,cage_result ))
#result[is.na(result)] <- 0
colnames(result)[1] <- "symbol"

final22 <- merge(result, locus)
write.table(final22, "smr_gwas_snp_summarized.txt", quote = F, row.names = F, sep = "\t")

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
write.table(final3, "smr_gwas_snp_summarized_nbh.txt", quote = F, row.names = F, sep = "\t")
