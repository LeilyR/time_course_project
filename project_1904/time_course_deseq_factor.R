
# time course analysis with deseq2


.libPaths("/data/processing3/leily/miniconda3/envs/snakepipes_RNAseq_environment_0.1/lib/R/library")
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

out_path = "/data/akhtar/group/rabbani/rna_project1904/lrt_factor/"


# Read count matrix
# count exons and summarize on gene level
# filtering count data!! TODO
countdata <- read.table("/data/manke/group/rabbani/rna_project1904/counts.tsv", header=TRUE, check.names = TRUE)

# Read smaplesheet , detects conditions, generate the formula
sampleInfo <- read.table("/data/manke/group/rabbani/rna_project1904/samplesheet.tsv",
												 header=TRUE, check.names = TRUE, stringsAsFactor = F)
sampleInfo$condition <- factor(sampleInfo$condition, levels = c('shCTRL', 'shMOF', 'shPRDX1'))
sampleInfo$treatment <- factor(sampleInfo$treatment, levels = c("LPS0", "LPS3", "LPS12"))
treatment_nested <- unique(sampleInfo$treatment.nested)
sampleInfo$fullname <- paste(sampleInfo$condition,sampleInfo$treatment, sampleInfo$replicates, sep='_')

# build the matrix

# from https://support.bioconductor.org/p/62684/#67807

d<-as.formula(~condition + treatment + condition:treatment)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata[,sampleInfo$fullname], colData = sampleInfo, design =d)
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds, test="LRT", reduced = ~condition+ treatment)
#resultsNames(dds)
res_lrt <- DESeq2::results(dds)
sig_res_lrt <- res_lrt %>%
               data.frame() %>%
               rownames_to_column(var="gene") %>%
               as_tibble() %>%
               filter(padj < 0.05)

write.table(sig_res_lrt, file = paste0(out_path,'time_course_sigGenes_lrt_factor.tsv', sep = ""),
            quote=FALSE, sep='\t', row.names=FALSE)

# Get sig gene lists
sig_lrt_genes <- sig_res_lrt %>%
                 pull(gene)
# pheatmap
normalized_counts <- counts(dds, normalized=T)
res_sig <- subset(res_lrt, padj<0.05)
res_sig_sorted = res_sig[order(res_sig$padj), ]
norm_sig <- subset(normalized_counts, row.names(normalized_counts) %in% row.names(res_sig_sorted))
norm_sig <- norm_sig[,c("shCTRL_LPS0_Rep1", "shCTRL_LPS0_Rep2", "shCTRL_LPS0_Rep3",
												"shMOF_LPS0_Rep1", "shMOF_LPS0_Rep2", "shMOF_LPS0_Rep3",
												"shPRDX1_LPS0_Rep1", "shPRDX1_LPS0_Rep2", "shPRDX1_LPS0_Rep3",
                        "shCTRL_LPS3_Rep1", "shCTRL_LPS3_Rep2", "shCTRL_LPS3_Rep3",
												"shMOF_LPS3_Rep1", "shMOF_LPS3_Rep2", "shMOF_LPS3_Rep3",
												"shPRDX1_LPS3_Rep1", "shPRDX1_LPS3_Rep2", "shPRDX1_LPS3_Rep3",
                        "shCTRL_LPS12_Rep1", "shCTRL_LPS12_Rep2", "shCTRL_LPS12_Rep3",
												"shMOF_LPS12_Rep1", "shMOF_LPS12_Rep2", "shMOF_LPS12_Rep3",
												"shPRDX1_LPS12_Rep1", "shPRDX1_LPS12_Rep2", "shPRDX1_LPS12_Rep3")]
write.table(normalized_counts, file=paste0(out_path,'time_course_normCount_lrt_factor.tsv', sep = ""), quote=FALSE, sep='\t', row.names=TRUE)
png(file = paste0(out_path,'time_course_normCountOnSigGenes_lrt_factor.png', sep = ""))
pheatmap(norm_sig,
         cluster_rows = T,
         show_rownames = F,
         border_color = NA,
         fontsize = 10,
         scale = "row",
				 # cutree_rows = 10,
         # fontsize_row = 10,
         height = 20)
dev.off()

# Subset the result (optional)
# 3000
clustering_sig_genes <- sig_res_lrt %>%
                        arrange(padj)
              # %>% head(n=100)
#head(clustering_sig_genes)
# Obtain rlog values for those significant genes
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
row.names(sampleInfo) <- paste(sampleInfo$condition, sampleInfo$treatment, sampleInfo$replicates, sep = "_")
clusters <- degPatterns(cluster_rlog, metadata = sampleInfo, time="treatment", col="condition", plot = FALSE, minc = 5)
for (ii in unique(clusters[["normalized"]]$cluster)) {
	cluster_id <- paste(ii ,sep='')
	cluster_id <- substr(cluster_id,1,nchar(cluster_id))
	num_genes <- clusters[["normalized"]] %>% filter(cluster == ii)
	num_genes <- unique(num_genes[["genes"]])
	num_genes <- length(num_genes)
	filename_to_save <- paste(out_path,"cluster_figures/",cluster_id,".png", sep ="")
	png(file = filename_to_save)
	plt<- ggplot(clusters[["normalized"]] %>% filter(cluster == ii),
	 		  aes(treatment, value, fill = condition)) +
	  		geom_boxplot() +
	 			# geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
	      # change the method to make it smoother
		    # geom_smooth(aes(group=condition), method = "lm") +
		    ggtitle(paste("cluster: ", cluster_id, " num_genes: ", num_genes, sep = ""))
	print(plt)
	dev.off() #TODO ggsave()

  ## save genes of certain clusters
  if (ii  %in% c(2, 14, 15, 18, 25, 30, 33, 35, 37, 38, 42, 45, 46, 49, 51, 55, 56, 58, 59,61, 63, 68, 69,70,73, 74, 76, 86, 87, 91)){
    df <- clusters$df %>% filter(cluster == ii)
    df$gene_id <- sapply(strsplit(as.character(df$genes),'.', fixed = TRUE), "[", 1)
    cluster_file_to_Save <- paste(out_path,"genes_in_clusters/",cluster_id,".tsv", sep ="")
    write.table(df, file=cluster_file_to_Save, quote=FALSE, sep='\t', row.names=FALSE)
  }
}


ggplt<- degPlotCluster(clusters$normalized, time="treatment", col="condition", lines = FALSE, points = FALSE)
ggplt + theme(text = element_text(size=6))
ggsave(filename = paste0(out_path,"cluster_figures/cluster_rlt_sig_genes_factor_all.pdf", sep =""), plot = ggplt, width = 32, height = 32, units= "cm")
