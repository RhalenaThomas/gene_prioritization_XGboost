library(data.table)
library(dplyr)

feature <- fread("xgb_fea_imp_v9_imbalanced_v2.csv")

label <- fread("training_labels.txt")
dist <- fread("../dist/dist_data.txt")
tss <- fread("../dist/tss_data.txt")
vep <- fread("../VEP/vep_data.txt")
eqtl_mean <- fread("../FUMA/eqtl_gwas_snp_summarized_mean.txt")
eqtl_max <- fread("../FUMA/eqtl_gwas_snp_summarized_max.txt")
eqtl_mean_nbh <- fread("../FUMA/eqtl_gwas_snp_summarized_mean_nbh.txt")
eqtl_max_nbh <- fread("../FUMA/eqtl_gwas_snp_summarized_max_nbh.txt")
cuomo <- fread("../Cuomo_eQTL_scRNA_DA_ipsc/D30.DA.qtl_results_all.sorted_gwas_snp_summarized.txt")
cuomo_nbh <- fread("../Cuomo_eQTL_scRNA_DA_ipsc/D30.DA.qtl_results_all.sorted_gwas_snp_summarized_nbh.txt")
ci <- fread("../FUMA/ci_gwas_genes_hgnc.txt")
brain_eqtl <- fread("../Bryois2021_eQTL_brain/brain_eqtl_gwas_snp_summarized.txt")
brain_eqtl_nbh <- fread("../Bryois2021_eQTL_brain/brain_eqtl_gwas_snp_summarized_nbh.txt")
smr_qtl <- fread("../SMR_QTL/smr_gwas_snp_summarized.txt")
smr_qtl_nbh <- fread("../SMR_QTL/smr_gwas_snp_summarized_nbh.txt")
AGTR1_da <- fread("../AGTR1_SNpc/AGTR1_da_mean_expression.txt")
AGTR1_da_nbh <- fread("../AGTR1_SNpc/AGTR1_da_mean_expression_nbh.txt")

dist_data <- dist[,c("symbol","DIST_MEAN","DIST_MIN","DIST_MEAN_NBH","DIST_MIN_NBH","locus")]
tss_data <- tss[,c("symbol","TSS_MEAN","TSS_MIN","TSS_MEAN_NBH","TSS_MIN_NBH","locus")]
vep_data <- vep[,c("symbol","VEP_IMPACT_MAX","POLY_MAX","VEP_IMPACT_MAX_NBH","POLY_MAX_NBH","locus")]
ci$ENHANCER_PROMOTER <- 1
ci_data <- ci[,c("symbol","ENHANCER_PROMOTER","locus")]

t1 <- merge(dist_data,tss_data, by = c("symbol","locus"), all = TRUE)
t1[is.na(t1)] <- -10

eqtl_mean <- eqtl_mean[ , EQTL_ALL_MEAN_EFFECT:=NULL]
eqtl_max <- eqtl_max[ , EQTL_ALL_MAX_EFFECT:=NULL]
eqtl_mean_nbh <- eqtl_mean_nbh[ , EQTL_ALL_MEAN_EFFECT_NBH:=NULL]
eqtl_max_nbh <- eqtl_max_nbh[ , EQTL_ALL_MAX_EFFECT_NBH:=NULL]


all_eqtl <- Reduce(function(x, y) merge(x, y, by = c("symbol","locus"), all.x = T), list(eqtl_mean,eqtl_max,cuomo,brain_eqtl,smr_qtl))

mean_eqtl <- all_eqtl[,which(grepl("*MEAN*", colnames(all_eqtl))), with=FALSE]
max_eqtl <- all_eqtl[,which(grepl("*MAX*", colnames(all_eqtl))), with=FALSE]

all_eqtl$EQTL_ALL_MEAN_EFFECT <- rowMeans(mean_eqtl, na.rm = T)
all_eqtl$EQTL_ALL_MAX_EFFECT <- rowMeans(mean_eqtl, na.rm = T)

mean_eqtl_nbh.list <- lapply(sort(unique(all_eqtl$locus)), function(i){
                final2 <- subset(all_eqtl, locus == i)
                if(nrow(final2) == 1){
                        EQTL_ALL_MEAN_EFFECT_NBH <- final2$EQTL_ALL_MEAN_EFFECT
                        EQTL_ALL_MEAN_EFFECT_NBH[EQTL_ALL_MEAN_EFFECT_NBH > 0] <- 1
                        return(data.frame(symbol = final2$symbol,EQTL_ALL_MEAN_EFFECT_NBH, locus = final2$locus))
                }
                EQTL_ALL_MEAN_EFFECT_NBH <- abs(final2$EQTL_ALL_MEAN_EFFECT)/max(abs(final2$EQTL_ALL_MEAN_EFFECT))

                return(data.frame(symbol = final2$symbol,EQTL_ALL_MEAN_EFFECT_NBH, locus = final2$locus))
})

max_eqtl_nbh.list <- lapply(sort(unique(all_eqtl$locus)), function(i){
                final2 <- subset(all_eqtl, locus == i)
                if(nrow(final2) == 1){
                        EQTL_ALL_MAX_EFFECT_NBH <- final2$EQTL_ALL_MAX_EFFECT
                        EQTL_ALL_MAX_EFFECT_NBH[EQTL_ALL_MAX_EFFECT_NBH > 0] <- 1
                        return(data.frame(symbol = final2$symbol,EQTL_ALL_MAX_EFFECT_NBH, locus = final2$locus))
                }
                EQTL_ALL_MAX_EFFECT_NBH <- abs(final2$EQTL_ALL_MAX_EFFECT)/max(abs(final2$EQTL_ALL_MAX_EFFECT))

                return(data.frame(symbol = final2$symbol,EQTL_ALL_MAX_EFFECT_NBH, locus = final2$locus))
})

eqtl <- Reduce(function(x, y) merge(x, y, by = c("symbol","locus"), all.x = T), list(all_eqtl[,c("symbol","locus","EQTL_ALL_MEAN_EFFECT","EQTL_ALL_MAX_EFFECT"), with = F],
                                                                                Reduce(rbind, mean_eqtl_nbh.list),
                                                                                Reduce(rbind, max_eqtl_nbh.list)))


#final <- Reduce(function(x, y) merge(x, y, by = c("symbol","locus"), all.x = T), list(t1,vep_data,eqtl_mean,eqtl_max,eqtl_mean_nbh,eqtl_max_nbh,cuomo,cuomo_nbh,ci_data,brain_eqtl,brain_eqtl_nbh,smr_qtl,smr_qtl_nbh,AGTR1_da,AGTR1_da_nbh,eqtl))

final <- Reduce(function(x, y) merge(x, y, by = c("symbol","locus"), all.x = T), list(t1,vep_data,eqtl_mean,eqtl_max,eqtl_mean_nbh,eqtl_max_nbh,cuomo,cuomo_nbh,ci_data,brain_eqtl,brain_eqtl_nbh,smr_qtl,smr_qtl_nbh,AGTR1_da,AGTR1_da_nbh))

final_feature <- c(1,2, which(colnames(final) %in% feature$feature))

final2 <- merge(final[, ..final_feature],label, by = "symbol", all.y = T)


#final2[is.na(final2)] <- 0

is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))

final2[is.nan(final2)] <- NA

write.csv(final2, "training_data_v10_imbalanced_v2.csv", quote = F, row.names = F)
