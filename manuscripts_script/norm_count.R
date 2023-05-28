# time course analysis with deseq2


.libPaths("/path/to/lib/R/library")
suppressPackageStartupMessages({
library("DESeq2")
library("dplyr")
library("tibble")
library("ggplot2")
library(SummarizedExperiment)
library("DEGreport")
library(splines)
})
# In this case, using the likelihood ratio test with a reduced model
# which does not contain the interaction terms will test
# whether the condition induces a change in gene expression at any
# time point after the reference level time point (time 0)
out_path = "/path/to/normcount/"

# Read count matrix
# count exons and summarize on gene level
countdata <- read.table("/path/to/counts.tsv", header=TRUE, check.names = TRUE, row.names = 1)
# remove rows with mean <= 5
countdata = countdata[rowMeans(countdata)>5,]

# Read smaplesheet , detects conditions, generate the formula
sampleInfo <- read.table("/path/to/metadata.tsv", header=TRUE, check.names = TRUE,
                         stringsAsFactor = F)
sampleInfo$condition <- factor(sampleInfo$condition, levels = c('shCTRL', 'shMOF', 'shPRDX1'))
sampleInfo$time <- as.integer(sampleInfo$time)

sampleInfo$time <-ns(sampleInfo$time)[1:27,]
sampleInfo$treatment <- factor(sampleInfo$treatment, levels = c("LPS0", "LPS3", "LPS12"))

# build the matrix
d<-as.formula(~condition + time + condition:time)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata[,sampleInfo$name],
                                      colData = sampleInfo, design =d)
dds <- estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds, test="LRT", reduced = ~condition+ time)

# norm count
normalized_counts <- counts(dds, normalized=T)
write.table(normalized_counts, file=paste0(out_path,'time_course_normCount_lrt_numeric.tsv', sep = ""),
            quote=FALSE, sep='\t', row.names=TRUE)
