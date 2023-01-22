.libPaths(c("~/runs/eyu8/library/Seurat_v4/","/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))

library(Seurat)

data <- readRDS("~/runs/go_lab/FOUNDIN-PD/SCRN/iPSCsDopaALL_integratedAfterBroadCellType.RDS")

genes_set1 <- unique(readLines("~/runs/eyu8/data/gene_prio/AGTR1_SNpc/gene_set1.txt"))
genes_set2 <- unique(readLines("~/runs/eyu8/data/gene_prio/AGTR1_SNpc/gene_set2.txt"))

genes_set1 <- genes_set1[which(genes_set1 %in% rownames(data))]
genes_set2 <- genes_set2[which(genes_set2 %in% rownames(data))]


da_data <- subset(x = data, idents = "Dopaminergic Neurons")

for(clust in c("2","1")){
        cells.use <- WhichCells(da_data, expression = pheno == clust)
        Idents(object = da_data, cells = cells.use) <- clust
}

ida_data <- subset(x = data, idents = "Immature Dopaminergic Neurons")

for(clust in c("2","1")){
        cells.use <- WhichCells(ida_data, expression = pheno == clust)
        Idents(object = ida_data, cells = cells.use) <- clust
}


set1 <- FindMarkers(data, ident.1 = "Dopaminergic Neurons", features = genes_set1, assay = "SCT", latent.vars = c("percent.mt","pheno","genetic_sex","BATCH"), test.use = "MAST")
write.csv(set1, "set1_Dopaminergic_Neurons_vs_all_sex_status_percentMT_batch.csv", quote = F)

set2 <- FindMarkers(data, ident.1 = "Dopaminergic Neurons", features = genes_set2, assay = "SCT", latent.vars = c("percent.mt","pheno","genetic_sex","BATCH"), test.use = "MAST")
write.csv(set2, "set2_Dopaminergic_Neurons_vs_all_sex_status_percentMT_batch.csv", quote = F)

set1 <- FindMarkers(da_data, ident.1 = "2", ident.2 = "1", features = genes_set1, assay = "SCT", latent.vars = c("percent.mt","genetic_sex","BATCH"), test.use = "MAST", logfc.threshold = 0)
write.csv(set1, "set1_PD_vs_Ctrl_da_sex_status_percentMT_batch.csv", quote = F)

set2 <- FindMarkers(da_data, ident.1 = "2", ident.2 = "1", features = genes_set2, assay = "SCT", latent.vars = c("percent.mt","genetic_sex","BATCH"), test.use = "MAST", logfc.threshold = 0)
write.csv(set2, "set2_PD_vs_Ctrl_da_sex_status_percentMT_batch.csv", quote = F)

set1 <- FindMarkers(ida_data, ident.1 = "2", ident.2 = "1", features = genes_set1, assay = "SCT", latent.vars = c("percent.mt","genetic_sex","BATCH"), test.use = "MAST", logfc.threshold = 0)
write.csv(set1, "set1_PD_vs_Ctrl_ida_sex_status_percentMT_batch.csv", quote = F)

set2 <- FindMarkers(ida_data, ident.1 = "2", ident.2 = "1", features = genes_set2, assay = "SCT", latent.vars = c("percent.mt","genetic_sex","BATCH"), test.use = "MAST", logfc.threshold = 0)
write.csv(set2, "set2_PD_vs_Ctrl_ida_sex_status_percentMT_batch.csv", quote = F)

