library(Matrix)
library(data.table)

full_matrix <- readMM("~/runs/eyu8/data/AGTR1_SNpc/GSE178265_Homo_matrix.mtx")

nrow(full_matrix)
ncol(full_matrix)

genes <- as.data.frame(fread("~/runs/eyu8/data/AGTR1_SNpc/GSE178265_Homo_features.tsv", header = FALSE))
cells <- readLines("~/runs/eyu8/data/AGTR1_SNpc/GSE178265_Homo_bcd.tsv")
colnames(full_matrix) <- cells
rownames(full_matrix) <- genes$V1

da <- fread("da_UMAP.tsv")

for(i in unique(da$Cell_Type)){
   data <- subset(da, Cell_Type == i)
   agtr1_col <- which(cells %in% data$NAME)
   agtr1_cluster <- full_matrix[,agtr1_col]
   write.table(rowMeans(agtr1_cluster),paste0(i,".txt"),quote = F, col.names = F, sep = "\t")
   }

