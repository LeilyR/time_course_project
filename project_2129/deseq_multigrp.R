# contrast of interest
# 1. LPS0: WT(experiment)_vector (control)
# 2. LPS0: K197R(experiment)_WT (control)
# 3. LPS0: K197Q(experiment)_WT (control)
# 4. LPS0: K197Q(experiment)_K197R (control)
#
# 5. LPS3: WT(experiment)_vector (control)
# 6. LPS3: K197R(experiment)_WT (control)
# 7. LPS3: K197Q(experiment)_WT (control)
# 8. LPS3: K197Q(experiment)_K197R (control)
#
# 9. LPS3: WT(experiment)_vector (control)
# 10. LPS3: K197R(experiment)_WT (control)
# 11. LPS3: K197Q(experiment)_WT (control)
# 12. LPS3: K197Q(experiment)_K197R (control)​

# 13. WT: lps3(experiment)_vs_lps0 (control)​
# 14. WT: lps12(experiment)_vs_lps0 (control)
# 15. WT: lps12(experiment)_vs_lps3 (control)
#
# 16. K197R: lps3(experiment)_vs_lps0 (control)​
# 17. K197R: lps12(experiment)_vs_lps0 (control)​
# 18. K197R: lps12(experiment)_vs_lps3 (control)​
#
# 19. K197Q: lps3(experiment)_vs_lps0 (control)​
# 20. K197Q: lps12(experiment)_vs_lps0 (control)​
# 21. K197Q: lps12(experiment)_vs_lps3 (control)​


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
out_path = "/data/akhtar/group2/rabbani/rna_project2129/pairwise_comparison_vs_wt/"
# Read count matrix
# count exons and summarize on gene level
countdata <- read.table("/data/manke/group/rabbani/rna_project2129/brb_counts.tsv",
												header=TRUE, check.names = FALSE)

# Read smaplesheet , detects conditions, generate the formula
sampleInfo <- read.table("/data/manke/group/rabbani/rna_project2129/samplesheet_pw.tsv",
												 header=TRUE, check.names = FALSE, stringsAsFactor = F)
sampleInfo$condition <- factor(sampleInfo$condition, levels = c('WT', 'vector', 'K197R', 'K197Q'))
sampleInfo$treatment <- factor(sampleInfo$treatment, levels = c("LPS0", "LPS3", "LPS12"))
lps_treats <- unique(sampleInfo$treatment)
# build the matrix
#d<-as.formula(~1 + group) ## Originally used this
# from https://www.biostars.org/p/395926/
d<-as.formula(~condition + treatment + condition:treatment)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata[,sampleInfo$name], colData = sampleInfo, design =d)
dds <- estimateSizeFactors(dds)
mod_mat <- model.matrix(design(dds), colData(dds))
dds <- DESeq2::DESeq(dds)

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

print(resultsNames(dds))

normalized_counts <- counts(dds, normalized=T)
write.table(normalized_counts, file=paste0(out_path,'normalized_counts_vs_wt.tsv',sep =""), quote=FALSE, sep='\t', row.names=TRUE)

# calculate coefficient vectors for each group (treatment)
mylist <- list()
all_contrasts <- list()
# lps_treats
for(treat in lps_treats){
	print(treat)
	sub <- subset(sampleInfo,  sampleInfo$treatment %in% treat)
	treat_nested <- sub$treatment.nested

	for (case in unique(treat_nested)){
		mylist[[case]] <- colMeans(mod_mat[dds$treatment.nested == case, ])
	}
	# obtain results for each pairwise contrast
	for(pair in list(a = c('vector', 'WT'), b = c('K197R', 'WT'), c = c('K197Q', 'WT'), d = c('K197Q', 'K197R'))){
		name1 = paste(pair[1],treat,sep=('_'))
		name2 = paste(pair[2],treat,sep=('_'))
		this_contrast <- DESeq2::results(dds, contrast = mylist[[name1]] - mylist[[name2]])
		saving_name = paste(name1,name2,sep=('_vs_'))
		all_contrasts[[saving_name]] <- this_contrast
	}
}
# contrast for same condition betwen time points
for(cond in c('WT', 'K197Q', 'K197R')){
	time0 = paste(cond, 'LPS0', sep=('_'))
	time3 = paste(cond, 'LPS3', sep=('_'))
	time12 = paste(cond, 'LPS12', sep=('_'))
	print(time0)
	contrast3_0 <- DESeq2::results(dds, contrast = mylist[[time3]] - mylist[[time0]])
	saving_name = paste(time3,time0,sep=('_vs_'))
	all_contrasts[[saving_name]] <- contrast3_0
	contrast12_0 <- DESeq2::results(dds, contrast = mylist[[time12]] - mylist[[time0]])
	saving_name = paste(time12,time0,sep=('_vs_'))
	all_contrasts[[saving_name]] <- contrast12_0
	contrast12_3 <- DESeq2::results(dds, contrast = mylist[[time12]] - mylist[[time3]])
	saving_name = paste(time12,time3,sep=('_vs_'))
	all_contrasts[[saving_name]] <- contrast12_3
}
print(mylist)

# Add gene symbols
gene_names <- read.table("/data/manke/group/rabbani/rna_project1904/genes.filtered.symbol", header=FALSE)
gene_names <- gene_names[!duplicated(gene_names[,1]),]
colnames(gene_names) <- c("GeneID", "external_gene_name")
print(head(gene_names))

for(case in 1:length(all_contrasts)){
		print(names(all_contrasts)[case])
		name_to_save = paste0(out_path,
													names(all_contrasts)[case],".tsv", sep="")
		ddr.df <- as.data.frame(all_contrasts[[case]])
		ddr.df$Status <- ifelse(is.na(ddr.df$padj), "None",
		                  	ifelse(ddr.df$padj < 0.05,
		                        ifelse(ddr.df$log2FoldChange > 0, "UP", "DOWN"),
		                        	"None"))
		ddr.df <- merge(ddr.df, gene_names, by.x = 0, by.y = "GeneID" , all.x = TRUE)
		colnames(ddr.df)[which(names(ddr.df) == "Row.names")] <- "GeneID"
		sub_up <- subset(ddr.df, ddr.df$Status %in% c('UP'))
		sub_down <- subset(ddr.df, ddr.df$Status %in% c('DOWN'))
		print(length(sub_up$padj))
		print(length(sub_down$padj))
		write.table(ddr.df, file=name_to_save, quote=FALSE, sep='\t')
	}



# lfcshrink?
