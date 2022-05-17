suppressPackageStartupMessages({
library("biomaRt")
})

option_list = list(
  make_option(c("-i", "--input"), default=NA, type='character',
              help="ensembl gene table, should have header, gene id column should be called 'gene_id' (tsv)"),
  make_option(c("-r", "--reference"), default=NA, type='character',
              help="reference genome (e.g. mmusculus_gene_ensembl)"),
  make_option(c("-o", "--output"), default=NA, type='character',
              help="kegg name gene table (tsv)")
)
opt = parse_args(OptionParser(option_list=option_list))
# "/data/akhtar/group2/rabbani/rna_project1904/pairwise_comparison/fdr0.05/gsea/gene_id_ensembl.tsv"
gene_names <- read.table(opt$input, header=TRUE, check.names = FALSE, stringsAsFactor = F)

gene_names$gene_id <- as.factor(gene_names$gene_id)

ensembl <- useEnsembl(biomart = "genes", version = 103)
# "mmusculus_gene_ensembl"
ensembl <- useDataset(dataset = opt$reference, mart = ensembl)


converted_gene_ids <- getBM (filters = c("ensembl_gene_id"),
                         		 attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "external_gene_name", "entrezgene_id", "entrezgene_description"),
                             values = gene_names$gene_id,
														 mart = ensembl)
# '/data/akhtar/group2/rabbani/rna_project1904/pairwise_comparison/fdr0.05/gsea/kegg_name.tsv'
write.table(converted_gene_ids, file= opt$out, quote=FALSE, sep='\t')
