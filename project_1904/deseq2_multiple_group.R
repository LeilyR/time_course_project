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
out_path = "/data/akhtar/group2/rabbani/rna_project1904/pairwise_comparison/"

# Read count matrix
# count exons and summarize on gene level
countdata <- read.table("/data/manke/group/rabbani/rna_project1904/brb_counts.tsv",
												header=TRUE, check.names = TRUE)

# remove rows with mean <= 5
countdata = countdata[rowMeans(countdata)>5,]

# Read smaplesheet , detects conditions, generate the formula
sampleInfo <- read.table("/data/manke/group/rabbani/rna_project1904/samplesheet.tsv",
												 header=TRUE, check.names = TRUE, stringsAsFactor = F)
sampleInfo$condition <- factor(sampleInfo$condition, levels = c('shCTRL', 'shMOF', 'shPRDX1'))
sampleInfo$treatment <- factor(sampleInfo$treatment, levels = c("LPS0", "LPS3", "LPS12"))
treatment_nested <- unique(sampleInfo$treatment.nested)
sampleInfo$fullname <- paste(sampleInfo$condition,sampleInfo$treatment, sampleInfo$replicates, sep='_')
# typeof(sampleInfo$fullname) is character! :)
# build the matrix
# from https://www.biostars.org/p/395926/
d<-as.formula(~condition + treatment + condition:treatment)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata[,sampleInfo$fullname], colData = sampleInfo, design =d)
dds <- estimateSizeFactors(dds)
mod_mat <- model.matrix(design(dds), colData(dds))
dds <- DESeq2::DESeq(dds)
print(resultsNames(dds))


# colMeans per treatment.nested
lps0_shctrl <- colMeans(mod_mat[dds$treatment.nested == "lps0_shctrl", ])
lps0_shmof <- colMeans(mod_mat[dds$treatment.nested == "lps0_shmof", ])
lps0_shprdx1 <- colMeans(mod_mat[dds$treatment.nested == "lps0_shprdx1", ])

lps3_shctrl <- colMeans(mod_mat[dds$treatment.nested == "lps3_shctrl", ])
lps3_shmof <- colMeans(mod_mat[dds$treatment.nested == "lps3_shmof", ])
lps3_shprdx1 <- colMeans(mod_mat[dds$treatment.nested == "lps3_shprdx1", ])

lps12_shctrl <- colMeans(mod_mat[dds$treatment.nested == "lps12_shctrl", ])
lps12_shmof <- colMeans(mod_mat[dds$treatment.nested == "lps12_shmof", ])
lps12_shprdx1 <- colMeans(mod_mat[dds$treatment.nested == "lps12_shprdx1", ])

# obtain results for each pairwise contrast

mylist <- list()

# (1) in each time point (LPS0 or LP3, LP12), the DEGs in shCTRL vs shPRDX1 vs shMOF

ddr_lps0_shmof_lps0_shctrl <- DESeq2::results(dds, contrast = lps0_shmof - lps0_shctrl)
mylist[[ "ddr_lps0_shmof_lps0_shctrl" ]] <- ddr_lps0_shmof_lps0_shctrl
ddr_lps0_shprdx1_lps0_shctrl <- DESeq2::results(dds, contrast = lps0_shprdx1- lps0_shctrl)
mylist[[ "ddr_lps0_shprdx1_lps0_shctrl" ]] <- ddr_lps0_shprdx1_lps0_shctrl


ddr_lps3_shmof_lps3_shctrl <- DESeq2::results(dds, contrast = lps3_shmof - lps3_shctrl)
mylist[[ "ddr_lps3_shmof_lps3_shctrl" ]] <- ddr_lps3_shmof_lps3_shctrl
ddr_lps3_shprdx1_lps3_shctrl <- DESeq2::results(dds, contrast = lps3_shprdx1- lps3_shctrl)
mylist[[ "ddr_lps3_shprdx1_lps3_shctrl" ]] <- ddr_lps3_shprdx1_lps3_shctrl


ddr_lps12_shmof_lps12_shctrl <- DESeq2::results(dds, contrast = lps12_shmof - lps12_shctrl)
mylist[[ "ddr_lps12_shmof_lps12_shctrl" ]] <- ddr_lps12_shmof_lps12_shctrl
ddr_lps12_shprdx1_lps12_shctrl <- DESeq2::results(dds, contrast = lps12_shprdx1- lps12_shctrl)
mylist[[ "ddr_lps12_shprdx1_lps12_shctrl" ]] <- ddr_lps12_shprdx1_lps12_shctrl


# (2) in each sample (shCTRL, shMOF, shPRDX11), the DEGs in LPS0, LPS3 , PS12

ddr_lps3_shmof_lps0_shmof <- DESeq2::results(dds, contrast = lps3_shmof - lps0_shmof)
mylist[[ "ddr_lps3_shmof_lps0_shmof" ]] <- ddr_lps3_shmof_lps0_shmof
ddr_lps12_shmof_lps0_shmof <- DESeq2::results(dds, contrast = lps12_shmof - lps0_shmof)
mylist[[ "ddr_lps12_shmof_lps0_shmof" ]] <- ddr_lps12_shmof_lps0_shmof
ddr_lps12_shmof_lps3_shmof <- DESeq2::results(dds, contrast = lps12_shmof - lps3_shmof)
mylist[[ "ddr_lps12_shmof_lps3_shmof" ]] <- ddr_lps12_shmof_lps3_shmof

ddr_lps3_shprdx1_lps0_shprdx1 <- DESeq2::results(dds, contrast = lps3_shprdx1 - lps0_shprdx1)
mylist[[ "ddr_lps3_shprdx1_lps0_shprdx1" ]] <- ddr_lps3_shprdx1_lps0_shprdx1
ddr_lps12_shprdx1_lps0_shprdx1 <- DESeq2::results(dds, contrast = lps12_shprdx1 - lps0_shprdx1)
mylist[[ "ddr_lps12_shprdx1_lps0_shprdx1" ]] <- ddr_lps12_shprdx1_lps0_shprdx1
ddr_lps12_shprdx1_lps3_shprdx1 <- DESeq2::results(dds, contrast = lps12_shprdx1 - lps3_shprdx1)
mylist[[ "ddr_lps12_shprdx1_lps3_shprdx1" ]] <- ddr_lps12_shprdx1_lps3_shprdx1

ddr_lps3_shctrl_lps0_shctrl <- DESeq2::results(dds, contrast = lps3_shctrl - lps0_shctrl)
mylist[[ "ddr_lps3_shctrl_lps0_shctrl" ]] <- ddr_lps3_shctrl_lps0_shctrl
ddr_lps12_shctrl_lps0_shctrl <- DESeq2::results(dds, contrast = lps12_shctrl - lps0_shctrl)
mylist[[ "ddr_lps12_shctrl_lps0_shctrl" ]] <- ddr_lps12_shctrl_lps0_shctrl
ddr_lps12_shctrl_lps3_shctrl <- DESeq2::results(dds, contrast = lps12_shctrl - lps3_shctrl)
mylist[[ "ddr_lps12_shctrl_lps3_shctrl" ]] <- ddr_lps12_shctrl_lps3_shctrl

# Add gene symbols
gene_names <- read.table("/data/manke/group/rabbani/rna_project1904/genes.filtered.symbol", header=FALSE)
gene_names <- gene_names[!duplicated(gene_names[,1]),]
colnames(gene_names) <- c("GeneID", "external_gene_name")
print(head(gene_names))

# write all files
for(case in 1:length(mylist)){
		name_to_save = paste0(out_path,
													names(mylist)[case],".tsv", sep="")
		print(name_to_save)
		ddr.df <- as.data.frame(mylist[[case]])
		ddr.df$Status <- ifelse(is.na(ddr.df$padj), "None",
		                  	ifelse(ddr.df$padj < 0.05,
		                        ifelse(ddr.df$log2FoldChange > 0, "UP", "DOWN"),
		                        	"None"))
		ddr.df <- merge(ddr.df, gene_names, by.x = 0, by.y = "GeneID" , all.x = TRUE)
		colnames(ddr.df)[which(names(ddr.df) == "Row.names")] <- "GeneID"
		write.table(ddr.df, file=name_to_save, quote=FALSE, sep='\t')
	}



# lfcshrunk?
