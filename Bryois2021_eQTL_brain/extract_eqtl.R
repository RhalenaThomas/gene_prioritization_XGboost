library(data.table)
library(dplyr)

astro <- fread("gwas_snp_Astrocytes.txt")
endo <- fread("gwas_snp_Endothelial.txt")
exc <- fread("gwas_snp_Excitatory.txt")
inhib <- fread("gwas_snp_Inhibitory.txt")
micro <- fread("gwas_snp_Microglia.txt")
oli <- fread("gwas_snp_Oligodendrocytes.txt")
opc <- fread("gwas_snp_OPCs.txt")
peri <- fread("gwas_snp_Pericytes.txt")

locus <- fread("../locus/gwas_locus_hg19_genes.tab")
locus <- locus[,c("symbol","locus")]


#Extract effect size and calcuate mean and max values
result <- lapply(1:8, function(i){
        data <- list(astro,endo,exc,inhib,micro,oli,opc,peri)
        names(data) <- c("astrocytes","endothelial","excitatory","inhibitory","microglia","oligodendrocytes","opc","pericyte")
        colnames(data[[i]]) <- c("symbol","ensembl","SNP","dist","p","b")
        final <- data[[i]] %>% group_by(symbol) %>% summarize(MEAN_EFFECT = mean(abs(b)), MAX_EFFECT = max(abs(b)))
        colnames(final)[-1] <- paste0(names(data)[i],"_",colnames(final)[-1])
        return(final)
})

final <- Reduce(function(x, y) merge(x, y, by = "symbol", all = T), result)
#final[is.na(final)] <- 0
final$BRAIN_EQTL_MEAN_EFFECT <- rowMeans(final[,-1])
final$BRAIN_EQTL_MAX_EFFECT <- apply(final[,-1], 1, max)
final22 <- merge(final, locus)
write.table(final22, "brain_eqtl_gwas_snp_summarized.txt", quote = F, row.names = F, sep = "\t")

#Create neighborhood scores
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
#final3[is.na(final3)] <- 0
write.table(final3, "brain_eqtl_gwas_snp_summarized_nbh.txt", quote = F, row.names = F, sep = "\t")
