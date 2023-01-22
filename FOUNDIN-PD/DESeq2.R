.libPaths(c("~/runs/eyu8/library/DESeq2/","/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))

library(DESeq2)
library(data.table)
library(plyr)
library(ggplot2)

cts <- as.matrix(fread("~/runs/go_lab/FOUNDIN-PD/RNAB/aggregated_expression/countTable.tsv"),rownames=1)
set <- as.data.frame(fread("~/runs/eyu8/data/gene_prio/FOUNDIN-PD/gene_ensembl_conversion.txt"))
coldata <- data.frame(fread("~/runs/eyu8/data/gene_prio/FOUNDIN-PD/coldata.txt"),row.names=1)

rownames(cts) <- substr(rownames(cts),1,15)
cts <- cts[rownames(cts) %in% set$ensembl_gene_id, rownames(coldata)]
conv <- join(data.frame(ensembl_gene_id=rownames(cts)),set)
rownames(cts) <- conv$hgnc_symbol
coldata[which(coldata$Status == "GENPD" & coldata$Date == "da65"),]$Status <- "PDwVarda65"
#coldata[which(coldata$Status == "GENUN" & coldata$Date == "da65"),]$Status <- "PDwVarda65"
coldata[which(coldata$Status == "PD" & coldata$Date == "da65"),]$Status <- "PDwVarda65"
coldata[which(coldata$Status == "HC" & coldata$Date == "da65"),]$Status <- "HCda65"
coldata$GENDER <- factor(coldata$GENDER)
coldata$BATCH <- factor(coldata$BATCH)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Status + GENDER + BATCH + ageBin)

dds <- DESeq(dds)
res <- results(dds, contrast=c("Status","PDwVarda65","HCda65"))
res

write.csv(res,"dge_result_PDwVar_vs_HC_da65.csv",quote=F)

coldata <- data.frame(fread("~/runs/eyu8/data/gene_prio/FOUNDIN-PD/coldata.txt"),row.names=1)
coldata[which(coldata$Status == "PD" & coldata$Date == "da65"),]$Status <- "PDda65"
coldata[which(coldata$Status == "HC" & coldata$Date == "da65"),]$Status <- "HCda65"
coldata$GENDER <- factor(coldata$GENDER)
coldata$BATCH <- factor(coldata$BATCH)


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Status + GENDER + BATCH + ageBin)

dds <- DESeq(dds)
res <- results(dds, contrast=c("Status","PDda65","HCda65"))
res

write.csv(res,"dge_result_PD_vs_HC_da65.csv",quote=F)


coldata <- data.frame(fread("~/runs/eyu8/data/gene_prio/FOUNDIN-PD/coldata.txt"),row.names=1)
coldata[which(coldata$Status == "PD" & coldata$Date == "da65"),]$Status <- "PDwPROda65"
coldata[which(coldata$Status == "PRODROMA" & coldata$Date == "da65"),]$Status <- "PDwPROda65"
coldata[which(coldata$Status == "GENPD" & coldata$Date == "da65"),]$Status <- "PDwPROda65"
coldata[which(coldata$Status == "GENUN" & coldata$Date == "da65"),]$Status <- "PDwPROda65"
coldata[which(coldata$Status == "HC" & coldata$Date == "da65"),]$Status <- "HCda65"
coldata$GENDER <- factor(coldata$GENDER)


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Status + GENDER + BATCH + ageBin + Extra_status)

dds <- DESeq(dds)
res <- results(dds, contrast=c("Status","PDwPROda65","HCda65"))
res

write.csv(res,"dge_result_PDwPROda65_vs_HC_da65.csv",quote=F)



coldata <- data.frame(fread("~/runs/eyu8/data/gene_prio/FOUNDIN-PD/coldata.txt"),row.names=1)
coldata[which(coldata$Status == "PD"),]$Status <- "PDwPRO"
coldata[which(coldata$Status == "PRODROMA"),]$Status <- "PDwPRO"
coldata[which(coldata$Status == "GENPD"),]$Status <- "PDwPRO"
coldata[which(coldata$Status == "GENUN"),]$Status <- "PDwPRO"

coldata$GENDER <- factor(coldata$GENDER)


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Status + GENDER + BATCH + ageBin + Date + Extra_status)

dds <- DESeq(dds)
res <- results(dds, contrast=c("Status","PDwPRO","HC"))
res

write.csv(res,"dge_result_PDwPRO_vs_HC.csv",quote=F)


coldata <- data.frame(fread("~/runs/eyu8/data/gene_prio/FOUNDIN-PD/coldata.txt"),row.names=1)
coldata$GENDER <- factor(coldata$GENDER)


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Status + GENDER + BATCH + ageBin + Date)

dds <- DESeq(dds)
res <- results(dds, contrast=c("Status","PD","HC"))
res

write.csv(res,"dge_result_PD_vs_HC.csv",quote=F)
