library(data.table)
library(biomaRt)


locus <- fread("../locus/locus", header = FALSE)
colnames(locus) <-  c("symbol","locus")


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'), filters = 'hgnc_symbol', values = locus$symbol, mart = ensembl)

colnames(gene)[4] <- "symbol"
final <- merge(locus, gene)

nrow(subset(final, chromosome_name %in% 1:22))

library(data.table)
library(dplyr)
library(biomaRt)

data <- fread("pre_parse_vep_output.txt")
locus <- fread("../locus/gwas_locus_hg19_genes.tab")


imp <- strsplit(data$Extra, ";")
imp_high_index <- which(sapply(imp, function(e) 'IMPACT=HIGH' %in% e ))
imp_moderate_index <- which(sapply(imp, function(e) 'IMPACT=MODERATE' %in% e ))
imp_low_index <- which(sapply(imp, function(e) 'IMPACT=LOW' %in% e ))

cons <- strsplit(data$Consequence, ",")

cons_index <- which(sapply(cons, function(e) 'intron_variant' %in% e | "downstream_gene_variant" %in% e | "upstream_gene_variant" %in% e | "NMD_transcript_variant" %in% e))

data$VEP_IMPACT <- 0
data[imp_high_index,]$VEP_IMPACT <- 1
data[imp_moderate_index,]$VEP_IMPACT <- 0.66
data[imp_low_index,]$VEP_IMPACT <- 0.33
data[cons_index,]$VEP_IMPACT <- 0.1

data$POLY <- 0
poly <- strsplit(data$Extra, ";")
poly_index <- which(sapply(sapply(poly, function(e) grepl('Poly',e) ), function(ee) is.element(TRUE, ee)))
poly_score_list <- sapply(poly[poly_index], function(e) grep('Poly', e, value = TRUE))
poly_score <- as.numeric(sapply(strsplit(poly_score_list, "="), function(e) e[2]))
data[poly_index,]$POLY <- poly_score


final <- data %>% group_by(Gene) %>% summarize(VEP_IMPACT_MAX = max(VEP_IMPACT), POLY_MAX = max(POLY))

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "GRCh37")
gene <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = final$Gene, mart = ensembl)
colnames(gene) <- c("Gene", "symbol")
vep <- merge(final, gene)

vep$VEP_IMPACT_MAX_NBH <- 0
vep$POLY_MAX_NBH <- 0

for(i in unique(locus$locus)){
        gene <- subset(locus, locus == i)
        vep_sub <- which(vep$symbol %in% gene$symbol)
        max_impact <- max(vep[vep_sub,]$VEP_IMPACT_MAX)
        max_poly <- max(vep[vep_sub,]$POLY_MAX)
        vep[vep_sub,]$VEP_IMPACT_MAX_NBH <- vep[vep_sub,]$VEP_IMPACT_MAX/max_impact
        vep[vep_sub,]$POLY_MAX_NBH <- vep[vep_sub,]$POLY_MAX/max_poly
}

is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))

vep[is.nan(vep)] <- 0
locus_genes <- locus[,c("symbol","locus")]
write.table(merge(subset(vep,symbol > 0), locus_genes), "vep_data.txt", sep = "\t", quote = F, row.names = F)
