library(data.table)
library(dplyr)
library(reshape)

data <- fread("eqtl_gwas_snp.txt")
locus <- fread("../locus/gwas_locus_hg19_genes.tab")
locus <- locus[,c("symbol","locus")]

final <- data %>% group_by(symbol, tissue) %>% summarize(MEAN_EFFECT = mean(abs(signed_stats)), MAX_EFFECT = max(abs(signed_stats)))
all_data <- data %>% group_by(symbol) %>% summarize(MEAN_EFFECT = mean(abs(signed_stats)), MAX_EFFECT = max(abs(signed_stats)))
colnames(all_data) <- c("symbol", "EQTL_ALL_MEAN_EFFECT", "EQTL_ALL_MAX_EFFECT")


is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))


result <- cast(final[,c("symbol","tissue","MEAN_EFFECT")], symbol~tissue, mean)
#result[is.nan(result)] <- 0

colnames(result)[-1] <- paste0(colnames(result)[-1],"_MEAN_EFFECT")
result <- merge(result, all_data[,c("symbol", "EQTL_ALL_MEAN_EFFECT")], all.x = TRUE)
write.table(merge(result, locus), "eqtl_gwas_snp_summarized_mean.txt", quote = F, row.names = F, sep = "\t")

final22 <- merge(result, locus)
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

final3 <- Reduce(rbind, final2_subset.list)
#final3[is.nan(final3)] <- 0
write.table(final3, "eqtl_gwas_snp_summarized_mean_nbh.txt", quote = F, row.names = F, sep = "\t")


result <- cast(final[,c("symbol","tissue","MAX_EFFECT")], symbol~tissue, mean)
#result[is.nan(result)] <- 0

colnames(result)[-1] <- paste0(colnames(result)[-1],"_MAX_EFFECT")
result <- merge(result, all_data[,c("symbol", "EQTL_ALL_MAX_EFFECT")], all.x = TRUE)
write.table(merge(result, locus) , "eqtl_gwas_snp_summarized_max.txt", quote = F, row.names = F, sep = "\t")

final22 <- merge(result, locus)
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

final3 <- Reduce(rbind, final2_subset.list)
#final3[is.nan(final3)] <- 0
write.table(final3, "eqtl_gwas_snp_summarized_max_nbh.txt", quote = F, row.names = F, sep = "\t")
