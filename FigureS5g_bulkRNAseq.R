###############################################################################
# Fig. S5g — bulk RNA-seq of sorted naive CD4/CD8 T cells (SPF vs GF)
#
# Covered panels:
#   Fig. S5g: z-score heatmaps of oxidative phosphorylation (OXPHOS), TCA cycle
#             and glycolysis genes that differ between SPF and GF naive T cells
#             (bulk RNA-seq, 2 replicates per group per cell type)
#
# Workflow:
#   1. gene_id -> gene_name mapping (GENCODE vM33)
#   2. CPM filter (CPM > 2 in >= 3 samples) + protein-coding genes only
#   3. DESeq2 SPF vs GF per cell type (CD4Tn / CD8Tn)
#   4. z-score heatmaps of the curated OXPHOS / TCA / glycolysis gene lists
#
# Input:
#   RNAseq.SPF_GF_naiveT.counts.genes.matrix  (gene-level raw counts; columns:
#       CD4_GF x2, CD4_SPF x2, CD8_GF x2, CD8_SPF x2)
#   gencode.vM33.primary_assembly.annotation.gtf
#   protein_coding_genecode_vM33_primary_22157g.csv
#   GO term gene lists under ref/go_term/
#
# Source notebook: 20240806_SPF_GF_naiveT.ipynb
###############################################################################

require(DESeq2)
require(dplyr)
require(ggplot2)
require(ggrepel)
require(edgeR)
require(rtracklayer)

## z-score per gene (row)
zscore_mat <- function(mat){
  Mean <- rowMeans(mat)
  Sd <- apply(mat, 1, sd)
  mat <- sweep(mat, 1, Mean, FUN = '-')
  mat <- sweep(mat, 1, Sd, FUN = '/')
  return(mat)
}

# =============================================================================
# 1. gene_id -> gene_name mapping
# =============================================================================
gff <- rtracklayer::readGFF("ref/gtf/gencode.vM33.primary_assembly.annotation.gtf")
mapid <- gff[gff$type == "gene", c("gene_id", "gene_name")]
rm(gff)
mapid$gene_id <- sapply(
  lapply(lapply(mapid$gene_id, function(x) stringr::str_split(x, pattern = '\\.')), '[[', 1),
  '[[', 1)

# =============================================================================
# 2. expression matrix: CPM filter + protein-coding genes
# =============================================================================
mat_raw <- read.table('RNAseq.SPF_GF_naiveT.counts.genes.matrix', sep = '\t',
                      row.names = 1, header = TRUE, stringsAsFactors = FALSE)
mat_raw <- edgeR::cpm(mat_raw)

## sample columns: GF = A, SPF = B
A_name <- 'GF'; B_name <- 'SPF'
A_CD4_idx <- c(1:2); B_CD4_idx <- c(3:4)  # CD4_GF / CD4_SPF
A_CD8_idx <- c(5:6); B_CD8_idx <- c(7:8)  # CD8_GF / CD8_SPF

## keep genes with CPM > 2 in at least 3 samples
cut_count <- 2
cut_cell <- 3
mat_pc <- mat_raw[rowSums(mat_raw[, c(A_CD4_idx, B_CD4_idx, A_CD8_idx, B_CD8_idx)] > cut_count) >= cut_cell, ]
rownames(mat_pc) <- plyr::mapvalues(rownames(mat_pc), from = mapid$gene_id, to = mapid$gene_name)

## keep protein-coding genes
list_pc <- read.csv('ref/gtf/protein_coding_genecode_vM33_primary_22157g.csv', row.names = 1)
list_pc <- list_pc$gene_name
mat_pc <- mat_pc[intersect(list_pc, rownames(mat_pc)), ]

cat(paste0("get ", dim(mat_pc)[1], " filtered (cpm > ", cut_count, " in at least ",
           cut_cell, " samples) protein-coding genes for next analysis"))

# =============================================================================
# 3. DESeq2, SPF vs GF (per cell type)
# =============================================================================
run_deseq2 <- function(mat_counts, idx, celltype){
  condition <- factor(c(rep("GF", 2), rep("SPF", 2)), levels = c("GF", "SPF"))
  colData <- data.frame(condition, row.names = colnames(mat_counts)[idx])
  dds <- DESeqDataSetFromMatrix(countData = round(mat_counts[, idx]),
                                colData = colData,
                                design = ~condition)
  rownames(dds) <- plyr::mapvalues(rownames(dds), from = mapid$gene_id, to = mapid$gene_name)
  dds <- DESeq(dds)
  res <- results(dds, tidy = TRUE)
  res <- res[order(res$padj), ]
  write.table(res, sprintf("figures/SPF_vs_GF_%s_DEG_DESeq2.csv", celltype),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")
  return(res)
}

## NOTE: DESeq2 is run on the raw count matrix (not the CPM matrix)
mat_counts <- read.table('RNAseq.SPF_GF_naiveT.counts.genes.matrix', sep = '\t',
                         row.names = 1, header = TRUE, stringsAsFactors = FALSE)
res_cd4 <- run_deseq2(mat_counts, c(A_CD4_idx, B_CD4_idx), 'CD4T')
res_cd8 <- run_deseq2(mat_counts, c(A_CD8_idx, B_CD8_idx), 'CD8T')

## gene-name version of the full matrix (for heatmaps)
mat_raw2 <- mat_counts
rownames(mat_raw2) <- plyr::mapvalues(rownames(mat_raw2), from = mapid$gene_id, to = mapid$gene_name)

# =============================================================================
# 4. Fig. S5g — z-score heatmaps of OXPHOS / TCA / glycolysis genes
#    gene lists: curated DEGs of each GO term (up in SPF / up in GF)
#    column order: CD8_GF, CD4_GF, CD8_SPF, CD4_SPF
# =============================================================================
oxy_go <- read.table('ref/go_term/GO_term_0006119_oxidative_phosphorylation.txt', sep = '\t')
tca_go <- read.table('ref/go_term/GO_term_0006099_tricarboxylic_acid_cycle.txt', sep = '\t')
gly_go <- read.table('ref/go_term/GO_term_0006096_glycolytic_process.txt', sep = '\t')

plot_gene_heatmap <- function(genes, out_pdf, width = 6, height = 2, flip = FALSE){
  df <- data.frame(zscore_mat(log2(mat_raw2[genes, c(A_CD8_idx, A_CD4_idx, B_CD8_idx, B_CD4_idx)] + 1)))
  df$Gene <- rownames(df)
  df_melt <- reshape2::melt(df, id.vars = "Gene")
  p <- df_melt %>%
    ggplot(aes(x = variable, y = factor(Gene, levels = genes), fill = value)) +
    geom_tile(color = 'transparent') +
    labs(fill = "Z-score") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_blank(),
          legend.position = "right",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_gradientn(colours = c('#0d0e84', 'white', '#ce3233'), oob = scales::squish,
                         name = 'Z-Score', limits = c(-2, 2),
                         breaks = c(-2, -1, 0, 1, 2), labels = c('-2', '-1', '0', '1', '2')) +
    xlab('') + ylab('')
  if (flip) p <- p + coord_flip()
  ggsave(out_pdf, p, width = width, height = height)
  return(p)
}

## OXPHOS (GO:0006119): up in SPF / up in GF
genes_oxy <- c(  # up in SPF
  '1700066M21Rik','Ak4','Antkmt','Atp5pb','Atpsckmt','Coa6','Coq9','Cyc1','Dguok','Dld','Dnajc30','Fxn','Gadd45gip1','Iscu','Macroh2a1','Msh2','mt-Co2','Mtch2','Ndufa8','Ndufa9','Ndufa10','Ndufa11','Ndufa12','Ndufa13','Ndufab1','Ndufaf1','Ndufb2','Ndufb3','Ndufb4','Ndufb5','Ndufb6','Ndufb7','Ndufb8','Ndufb9','Ndufb10','Ndufb11','Ndufc1','Ndufc2','Ndufs1','Ndufs2','Ndufs3','Ndufs4','Ndufs5','Ndufs6','Ndufs7','Ndufs8','Ndufv1','Ndufv2','Nipsnap2','Park7','Rhoa','Sdha','Sdhaf2','Sdhb','Sdhc','Sdhd','Slc25a51','Stoml2','Tafazzin','Tefm','Uqcc3','Uqcrc1','Uqcrc2','Vcp',
  # up in GF
  'Tnf','Slc25a33','Proca1','Pink1','Ndufv3','Ndufb1','Ndufa3','Ndufa2','Ndufa1','mt-Nd5','Cox7c','Abcd1')
plot_gene_heatmap(genes_oxy, 'figures/oxy_heatmap.pdf', width = 6, height = 2)

## TCA cycle (GO:0006099): up in SPF
genes_tca <- c('Aco2','Cs','Dlat','Dlst','Fh1','Idh1','Idh2','Idh3b','Idh3g','Ireb2','Mdh1','Mdh2','Ndufs4','Ogdh','Pdha1','Sdha','Sdhaf2','Sdhb','Sdhc','Sdhd','Sucla2','Suclg1','Suclg2')
plot_gene_heatmap(genes_tca, 'figures/tca_heatmap.pdf', width = 6, height = 2, flip = TRUE)

## glycolysis (GO:0006096): up in SPF / up in GF
genes_gly <- c(  # up in SPF
  'Aldoa','Aldoc','Arl2','Ddit4','Eif6','Eno1','Eno1b','Fkrp','Foxk1','Foxk2','Galk1','Gpi1','Hdac4','Hkdc1','Jmjd8','Lipa','Mlst8','Mlx','Mtch2','Ogt','Pfkfb1','Pfkfb2','Pgam1','Pgk1','Prxl2c','Rptor','Sirt6','Slc4a1','Tkfc','Ucp2',
  # up in GF
  'Zbtb7a','Slc4a4','Sik2','Psen1','Ppp2ca','Pfkl','Pfkp','Hk2','Gapdh','Galt','Ep300')
plot_gene_heatmap(genes_gly, 'figures/gly_heatmap.pdf', width = 6, height = 1, flip = TRUE)
