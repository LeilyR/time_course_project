
# time course analysis with deseq2


.libPaths("/data/processing3/leily/miniconda3/envs/snakepipes_RNAseq_environment_0.1/lib/R/library")
suppressPackageStartupMessages({
library("DESeq2")
library("dplyr")
library("tibble")
library("ggplot2")
library(SummarizedExperiment)
library("DEGreport")
library("pheatmap")
library(splines)
})
# In this case, using the likelihood ratio test with a reduced model which does not contain the interaction terms will test
# whether the condition induces a change in gene expression at any time point after the reference level time point (time 0)
out_path = "/data/akhtar/group2/rabbani/rna_project1904/lrt_numeric_wo_prdx/"

# Read count matrix
# count exons and summarize on gene level

countdata <- read.table("/data/manke/group/rabbani/rna_project1904/brb_counts.tsv", header=TRUE, check.names = TRUE, row.names = 1)

# Read smaplesheet , detects conditions, generate the formula
sampleInfo <- read.table("/data/manke/group/rabbani/rna_project1904/samplesheet_numeric_wo_prdx1.tsv", header=TRUE, check.names = TRUE,
                         stringsAsFactor = F)
sampleInfo$condition <- factor(sampleInfo$condition, levels = c('shCTRL', 'shMOF'))
sampleInfo$time <- as.integer(sampleInfo$time)

sampleInfo$time <-ns(sampleInfo$time)[1:18,]
sampleInfo$treatment <- factor(sampleInfo$treatment, levels = c("LPS0", "LPS3", "LPS12"))

# build the matrix

# from https://support.bioconductor.org/p/62684/#67807

d<-as.formula(~condition + time + condition:time)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata[,sampleInfo$name], colData = sampleInfo, design =d)
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds, test="LRT", reduced = ~condition+ time)
# resultsNames(dds)
res_lrt <- DESeq2::results(dds)
all_genes <- res_lrt %>%
               data.frame() %>%
               rownames_to_column(var="gene") %>%
               as_tibble()
write.table(all_genes, file=paste0(out_path,'time_course_allGenes_lrt_numeric.tsv', sep = ""), quote=FALSE, sep='\t',
                                   row.names=FALSE)

sig_res_lrt <- res_lrt %>%
               data.frame() %>%
               rownames_to_column(var="gene") %>%
               as_tibble() %>%
               filter(padj < 0.05)

write.table(sig_res_lrt, file=paste0(out_path,'time_course_sigGenes_lrt_numeric.tsv', sep = ""), quote=FALSE, sep='\t',
                                     row.names=FALSE)
#
# # Get sig gene lists
# sig_lrt_genes <- sig_res_lrt %>%
#                  pull(gene)
# pheatmap
normalized_counts <- counts(dds, normalized=T)
res_sig <- subset(res_lrt, padj<0.05)
res_sig_sorted = res_sig[order(res_sig$padj), ]
norm_sig <- subset(normalized_counts, row.names(normalized_counts) %in% row.names(res_sig_sorted))
norm_sig <- norm_sig[,c("shCTRL_LPS0_Rep1", "shCTRL_LPS0_Rep2", "shCTRL_LPS0_Rep3", "shMOF_LPS0_Rep1", "shMOF_LPS0_Rep2", "shMOF_LPS0_Rep3",
                        "shCTRL_LPS3_Rep1", "shCTRL_LPS3_Rep2", "shCTRL_LPS3_Rep3", "shMOF_LPS3_Rep1", "shMOF_LPS3_Rep2", "shMOF_LPS3_Rep3",
                        "shCTRL_LPS12_Rep1", "shCTRL_LPS12_Rep2", "shCTRL_LPS12_Rep3", "shMOF_LPS12_Rep1", "shMOF_LPS12_Rep2", "shMOF_LPS12_Rep3")]
write.table(normalized_counts, file=paste0(out_path,'time_course_normCount_lrt_numeric.tsv', sep = ""), quote=FALSE, sep='\t', row.names=TRUE)
png(file = paste0(out_path,'time_course_normCountOnSigGenes_lrt_numeric.png', sep = ""))
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


clustering_sig_genes <- sig_res_lrt %>%
                        arrange(padj)

# Obtain rlog values for those significant genes
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
row.names(sampleInfo) <-sampleInfo$name
clusters <- degPatterns(cluster_rlog, metadata = sampleInfo, time="treatment", col="condition", plot = FALSE, minc = 10)
for (ii in unique(clusters[["normalized"]]$cluster)) {
	cluster_id <- paste(ii ,sep='')
	cluster_id <- substr(cluster_id,1,nchar(cluster_id))
	num_genes <- clusters[["normalized"]] %>% filter(cluster == ii)
	num_genes <- unique(num_genes[["genes"]])
	num_genes <- length(num_genes)
  df <- clusters$df %>% filter(cluster == ii)
  df$gene_id <- sapply(strsplit(as.character(df$genes),'.', fixed = TRUE), "[", 1)
  cluster_file_to_Save <- paste(out_path,"genes_in_clusters/",cluster_id,".tsv", sep ="")
  write.table(df, file=cluster_file_to_Save, quote=FALSE, sep='\t', row.names=FALSE)
}


ggplt<- degPlotCluster(clusters$normalized, time="treatment", col="condition", lines = FALSE, points = FALSE)
ggplt + theme(text = element_text(size=6))
ggsave(filename = paste0(out_path,"cluster_figures/cluster_rlt_sig_genes_numeric_all.pdf", sep =""),
                         plot = ggplt, width = 32, height = 32, units= "cm")
