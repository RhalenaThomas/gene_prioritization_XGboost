library(data.table)

locus <- fread("../locus/gwas_locus_hg19_genes.tab")
locus <- locus[,c("symbol","locus")]

SOX6_AGTR1 <- fread("SOX6_AGTR1.txt")

result <- Reduce(function(x, y) merge(x, y, by = "symbol", all.x = T), list(locus,SOX6_AGTR1))


result[is.na(result)] <- 0

write.table(result, "SOX6_AGTR1_da_mean_expression.txt", quote = F, row.names = F, sep = "\t")


final22 <- result
final2_subset.list <- lapply(sort(unique(final22$locus)), function(i){
                final2 <- subset(final22, locus == i)
                colnames(final2) <- paste0(colnames(final2),"_NBH")
                if(nrow(final2) == 1){
                        final3 <- final2[,3]
                        final3[final3 > 0] <- 1
                        return(data.frame(symbol = final2$symbol,final3, locus = final2$locus))
                }
                final3 <- apply(final2[,3], 2, function(i){
                        abs(i)/max(abs(i))
                })
                return(data.frame(symbol = final2$symbol,final3, locus = final2$locus))
})

is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))


final3 <- Reduce(rbind, final2_subset.list)
final3[is.nan(final3)] <- 0


write.table(final3, "SOX6_AGTR1_da_mean_expression_nbh.txt", quote = F, row.names = F, sep = "\t")
