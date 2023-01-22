library(data.table)
library(biomaRt)
library(GenomicRanges)

snp <- fread("../FUMA/FUMA_input.txt")
colnames(snp) <-  c("SNP","chr","start","A1","A2","b","se","p")
locus_gene <- fread("../locus/gwas_locus_hg19_genes.tab")
locus <- fread("../locus/gwas_locus_hg19.tab")
locus$CHR <- paste0("chr",locus$CHR)

tss <- fread("~/runs/go_lab/gencode/gencode.v40lift37.annotation.tss.bed")
colnames(tss)[4] <- "ensembl_gene_id"
snp$chr <- paste0("chr",snp$chr)
tss$ensembl_gene_id <- sapply(strsplit(tss$ensembl_gene_id,"[.]"), `[`, 1)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "GRCh37")
gene <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = c('ensembl_gene_id', 'biotype'), values = list(tss$ensembl_gene_id, biotype="protein_coding"), mart = ensembl)
colnames(gene)[2] <- "symbol"

tss_gene <- merge(tss, merge(gene, locus_gene[,c("symbol","locus")]))


tss_data.list <- lapply(sort(unique(locus$Locus)), function(i){
                locus_subset <- subset(locus, Locus == i)
                snp_subset <- subset(snp, chr %in% locus_subset$CHR & start > locus_subset$START & start < locus_subset$END)
                tss_gene_subset <- subset(tss_gene, chr %in% locus_subset$CHR & locus == i)
                tss_gene_subset$TSS_MEAN <- sapply(tss_gene_subset$start, function(x) round(mean(abs(snp_subset$start - x))))
                tss_gene_subset$TSS_MIN <- sapply(tss_gene_subset$start, function(x) round(min(abs(snp_subset$start - x))))
                tss_gene_subset$TSS_MEAN_NBH <- -log((tss_gene_subset$TSS_MEAN + 1)/min((tss_gene_subset$TSS_MEAN + 1)))
                tss_gene_subset$TSS_MIN_NBH <- -log((tss_gene_subset$TSS_MIN + 1)/min((tss_gene_subset$TSS_MIN + 1)))
                return(tss_gene_subset)
})

dist_data.list <- lapply(sort(unique(locus$Locus)), function(i){
                locus_subset <- subset(locus, Locus == i)
                snp_subset <- subset(snp, chr %in% locus_subset$CHR & start > locus_subset$START & start < locus_subset$END)
                locus_gene_subset <- subset(locus_gene, chr %in% locus_subset$CHR & locus == i)
                locus_gene_subset$DIST_MEAN <- sapply(1:nrow(locus_gene_subset), function(x) {
                        dist <- sapply(1:nrow(snp_subset) , function(y) min(abs(locus_gene_subset[x,]$start - snp_subset[y,]$start), abs(locus_gene_subset[x,]$end - snp_subset[y]$start)))
                        dist[which(locus_gene_subset[x,]$start < snp_subset$start & snp_subset$start < locus_gene_subset[x,]$end)] <- 0
                        round(mean(dist))

                })
                locus_gene_subset$DIST_MIN <- sapply(1:nrow(locus_gene_subset), function(x) {
                        dist <- sapply(1:nrow(snp_subset) , function(y) min(abs(locus_gene_subset[x,]$start - snp_subset[y,]$start), abs(locus_gene_subset[x,]$end - snp_subset[y]$start)))
                        dist[which(locus_gene_subset[x,]$start < snp_subset$start & snp_subset$start < locus_gene_subset[x,]$end)] <- 0
                        round(min(dist))

                })

                locus_gene_subset$DIST_MEAN_NBH <- -log((locus_gene_subset$DIST_MEAN + 1)/min((locus_gene_subset$DIST_MEAN + 1)))
                locus_gene_subset$DIST_MIN_NBH <- -log((locus_gene_subset$DIST_MIN + 1)/min((locus_gene_subset$DIST_MIN + 1)))

                return(locus_gene_subset)
})

tss_data <- Reduce(rbind, tss_data.list)
dist_data <- Reduce(rbind, dist_data.list)


write.table(dist_data, "dist_data.txt", sep = "\t", quote = F, row.names = F)
write.table(tss_data, "tss_data.txt", sep = "\t", quote = F, row.names = F)

