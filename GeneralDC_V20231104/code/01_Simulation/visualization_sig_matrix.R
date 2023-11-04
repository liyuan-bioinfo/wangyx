#' @desc test novel method of de-convolution
#' @time 20231104
#' @author Li-Y


library(dplyr)
library(pheamap)

setwd("/aaa/zihanwu/yyyli2/projectx_General_Deconv")

# pheatmap for sig.matrix
sig_matrix_df = read.delim("01_datasets/data/Lung_celltype9_sig_matrix.txt",sep="\t",row.names=1,header=T)

head(sig_matrix_df)

pdf("01_datasets/data/Pheatmap_Lung_celltype9_sig_matrix.pdf",width=6,height = 8)
p2=pheatmap::pheatmap(sig_matrix_df,show_rownames = F,cluster_cols = F,cluster_rows = T,scale="row")
dev.off()