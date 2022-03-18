# add libraries
.libPaths("/data/processing3/leily/miniconda3/envs/snakepipes_RNAseq_environment_0.1/lib/R/library")
suppressPackageStartupMessages({
library("DESeq2")
library("dplyr")
library("tibble")
library("ggplot2")
library("pheatmap")
library(splines)
})
out_path = "/data/akhtar/group2/rabbani/rna_project1904/mofLPS0Only/"
# Read count matrix
# count exons and summarize on gene level
countdata <- read.table("/data/manke/group/rabbani/rna_project1904/brb_counts.tsv",
												header=TRUE, check.names = TRUE)
# remove rows with mean <= 5
countdata = countdata[rowMeans(countdata)>5,]

# Read smaplesheet , detects conditions, generate the formula
sampleInfo <- read.table("/data/manke/group/rabbani/rna_project1904/mofOnly/samplesheet_mof0.tsv",
												 header=TRUE, check.names = TRUE, stringsAsFactor = F)
sampleInfo$condition <- factor(sampleInfo$condition, levels = c('shCTRL', 'shMOF'))
sampleInfo$treatment <- factor(sampleInfo$treatment)
sampleInfo$fullname <- paste(sampleInfo$condition,sampleInfo$treatment, sampleInfo$replicates, sep='_')
# build the matrix
# from https://www.biostars.org/p/395926/
d<-as.formula(~ 1 + condition)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata[,sampleInfo$fullname], colData = sampleInfo, design =d)
dds <- estimateSizeFactors(dds)
mod_mat <- model.matrix(design(dds), colData(dds))
dds <- DESeq2::DESeq(dds)
print(resultsNames(dds))
# lfcshrunk
res <- lfcShrink(dds, coef="condition_shMOF_vs_shCTRL", type="apeglm")
name_to_save = paste0(out_path,"mof_vsctrl_lps0_lfcshrunk.tsv", sep="")
ddr.df <- as.data.frame(res)
ddr.df$Status <- ifelse(is.na(ddr.df$padj), "None",
										ifelse(ddr.df$padj < 0.05,
												ifelse(ddr.df$log2FoldChange > 0, "UP", "DOWN"),
													"None"))
write.table(ddr.df, file=name_to_save, quote=FALSE, sep='\t')
