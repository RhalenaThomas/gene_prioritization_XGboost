library(data.table)

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


final <- Reduce(function(x, y) merge(x, y, by = c("symbol","locus"), all.x = T), list(t1,vep_data,eqtl_mean,eqtl_max,eqtl_mean_nbh,eqtl_max_nbh,cuomo,cuomo_nbh,ci_data,brain_eqtl,brain_eqtl_nbh,smr_qtl,smr_qtl_nbh,AGTR1_da,AGTR1_da_nbh))


final[is.na(final)] <- 0

write.csv(final, "testing_data_v9.csv", quote = F, row.names = F)
