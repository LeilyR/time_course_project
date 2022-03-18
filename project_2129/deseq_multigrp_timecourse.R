
# time course analysis with deseq2

.libPaths("/localenv/rabbani/anaconda/miniconda3/envs/snakepipes_RNAseq_environment_0.1/lib/R/library")
suppressPackageStartupMessages({
library("DESeq2")
library("dplyr")
library("tibble")
library("ggplot2")
library(SummarizedExperiment)
library("DEGreport")
#library("ggrepel")
library("pheatmap")
library(splines)
})
# In this case, using the likelihood ratio test with a reduced model which does not contain the interaction terms will test
# whether the condition induces a change in gene expression at any time point after the reference level time point (time 0)
out_path = "/data/akhtar/group2/rabbani/rna_project2129/lrt/"

# Read count matrix
# count exons and summarize on gene level
countdata <- read.table("/data/manke/group/rabbani/rna_project2129/brb_counts.tsv", header=TRUE, check.names = FALSE)
# remove rows with mean <= 5
countdata = countdata[rowMeans(countdata)>5,]

# # Read smaplesheet , detects conditions, generate the formula
sampleInfo <- read.table("/data/manke/group/rabbani/rna_project2129/samplesheet.tsv", header=TRUE, check.names = FALSE)
sampleInfo$condition <- factor(sampleInfo$condition, levels = c('WT', 'vector', 'K197R', 'K197Q'))
sampleInfo$treatment <- factor(sampleInfo$time, level = c('0', '3', '12'))
sampleInfo$time <- as.integer(sampleInfo$time)
# print(ns(sampleInfo$time))
sampleInfo$time <-ns(sampleInfo$time)[1:36,]
print(sampleInfo)
print("----------")
print(sampleInfo$name)
# build the matrix

# from https://support.bioconductor.org/p/62684/#67807

d<-as.formula(~condition + time + condition:time)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata[,sampleInfo$name], colData = sampleInfo, design =d)
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds, test="LRT", reduced = ~condition + time)
# pca plot
rld <- rlog(dds)
pcaplt = plotPCA(rld, intgroup=c("condition","treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaplt, "percentVar"))
pcaplt$name<-make.names(pcaplt$treatment,unique=TRUE)
ggplot(pcaplt, aes(PC1, PC2, color = condition, shape = treatment)) +
        geom_hline(aes(yintercept = 0), colour = "grey") +
        geom_vline(aes(xintercept = 0), colour = "grey") +
        geom_point(size = 5) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle("PCA\n") +
        scale_shape_manual(values = c(0:18,33:17))
ggsave(filename=paste0(out_path, "pcaplot.png", sep = ""), device='png')

# resultsNames(dds)
res_lrt <- DESeq2::results(dds)
all_genes <- res_lrt %>%
               data.frame() %>%
               rownames_to_column(var="gene") %>%
               as_tibble()
write.table(all_genes, file=paste0(out_path,'time_course_res_lrt_vs_wt.tsv', sep =""), quote=FALSE, sep='\t', row.names=FALSE)

sig_res_lrt <- res_lrt %>%
               data.frame() %>%
               rownames_to_column(var="gene") %>%
               as_tibble() %>%
               filter(padj < 0.05)

write.table(sig_res_lrt, file=paste0(out_path,'time_course_sig_res_lrt_vs_wt.tsv', sep =""), quote=FALSE, sep='\t', row.names=FALSE)

# pheatmap of norm counts on sig genes
normalized_counts <- counts(dds, normalized=T)
write.table(normalized_counts, file=paste0(out_path,'normalized_counts_vs_wt.tsv',sep =""), quote=FALSE, sep='\t', row.names=TRUE)
res_sig <- subset(res_lrt, padj<0.05)
res_sig_sorted = res_sig[order(res_sig$padj), ]
norm_sig <- subset(normalized_counts, row.names(normalized_counts) %in% row.names(res_sig_sorted))
png(file = paste0(out_path,'norm_counts_of_sig_genes_vs_wt.png',sep=""))
pheatmap(norm_sig,
         cluster_rows = T,
         show_rownames = F,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         cluster_cols = F,
				 # cutree_rows = 10,
         # fontsize_row = 10,
         height = 20)
dev.off()
#
#
# clustering_sig_genes <- sig_res_lrt %>%
#                         arrange(padj)
#
# # Obtain rlog values for those significant genes
# rld <- rlog(dds, blind=T)
# rld_mat <- assay(rld)
# cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
# row.names(sampleInfo) <-sampleInfo$name
# clusters <- degPatterns(cluster_rlog, metadata = sampleInfo, time="treatment", col="condition", plot = FALSE, minc = 10)
# for (ii in unique(clusters[["normalized"]]$cluster)) {
# 	cluster_id <- paste(ii ,sep='')
# 	cluster_id <- substr(cluster_id,1,nchar(cluster_id))
# 	num_genes <- clusters[["normalized"]] %>% filter(cluster == ii)
# 	num_genes <- unique(num_genes[["genes"]])
# 	num_genes <- length(num_genes)
#   df <- clusters$df %>% filter(cluster == ii)
#   df$gene_id <- sapply(strsplit(as.character(df$genes),'.', fixed = TRUE), "[", 1)
#   cluster_file_to_Save <- paste("/data/processing1/leily/deseq_pairwise/numeric_time/genes_in_cluster",cluster_id,".tsv", sep ="")
#   write.table(df, file=cluster_file_to_Save, quote=FALSE, sep='\t', row.names=FALSE)
# }
#
#
# ggplt<- degPlotCluster(clusters$normalized, time="treatment", col="condition", lines = FALSE, points = FALSE)
# ggplt + theme(text = element_text(size=6))
# ggsave(filename = "/data/manke/group/rabbani/deseq_multi_grp/cluster_figures/cluster_rlt_sig_genes_numeric_all.pdf", plot = ggplt, width = 32, height = 32, units= "cm")
# #
# # # Extract the Group 1 genes
# # #cluster_groups <- clusters$df
# # #group1 <- clusters$df %>%
# # #          filter(cluster == 1)
# # # pdf("/data/manke/group/rabbani/deseq_multi_grp/heatmap_rlt_sig_genes_consensusCluster.pdf")
# # # pheatmap(sig_res_lrt,
# # #          cluster_rows = T,
# # #          kmeans_k = 10,
# # #          show_rownames = F,
# # #          border_color = NA,
# # #          fontsize = 10,
# # #          scale = "row",
# # #          fontsize_row = 10,
# # #          height = 20)
