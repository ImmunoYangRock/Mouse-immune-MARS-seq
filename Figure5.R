###############################################################################
# Figure 5 + Fig. S5/S6 — lymphoid functional competence (R part)
#
# Description: functional maturation of lymphoid subsets — DEG scatter,
# volcano plots and heatmaps, paired boxplots of Pdcd4/Eif4g2/Eif4a2,
# the lymphoid marker heatmap, HALLMARK GSEA, and metacell SPF-fraction
# vs metabolic pathway score correlations.
#
# Covered panels:
#   Fig. 5a  : DEG stat scatter, 8 lymphoid subsets x 3 time points (Batch_deg)
#   Fig. 5b-d: volcano plots of B_mature (b), T_naive (c), CD8T_cyto (d)
#   Fig. 5e  : paired boxplots of Pdcd4 / Eif4g2 / Eif4a2 (per-subset means,
#              paired t-test, symbols mark B / T / ILC lineages)
#   Fig. S5a : lymphoid marker heatmap (metacell top markers, minibulk z-score)
#   Fig. S5b-d: DEG heatmaps of B_mature (b), T_naive (c), CD8T_cyto (d)
#   Fig. S5e : HALLMARK GSEA dotplot of B_mature / CD8T_cyto / T_naive
#   Fig. S6a : metacell SPF fraction vs OXPHOS / glycolysis / TCA score (lymphoid)
#   Fig. S6b : same as S6a, myeloid (reconstructed, see note in that section)
#   also     : per-sample cutoff cell lists (input of Figure5.py pseudobulk),
#              DPT_Rag1 / DPT_Itm2a volcano + heatmap (related views)
#
# Input (relative to the v9/lym working directory):
#   lym/data/seurat/v9_lym_mc2seu.rds          (lymphoid Seurat object)
#   lym/output/seurat/lym_88180c_mc_info.csv   (cell-metacell info)
#   lym/output/scanpy/decoupler/*.csv          (pseudobulk DEG tables, Figure5.py)
#   lym/data/lym_annot_color.Rda               (col.21 etc.)
#   lym/data/mc1 (metacell scdb, id 'mouse_v9_lym')  for Fig. S5a and S6a
#   mye/data/mc1 (metacell scdb, id 'mouse_v9_mye'),
#   mye/output/seurat/mye_32469c_mc_info.csv         for Fig. S6b
#   GO term gene lists under ref/go_term/
#
# Source notebooks: code store/lym/lym_down_heatmap.ipynb,
#                   code store/lym/lym_R_analysis.ipynb,
#                   code store/lym/mouse_v9_lym_mc1.ipynb
###############################################################################

require(Seurat)
require(ComplexHeatmap)
require(tidyverse)
require(ggplot2)
require(ggrepel)
require(patchwork)
require(aplot)
require(foreach)
require(stringr)

`%out%` <- Negate(`%in%`)

Time.color <- c('P16' = '#DFC1DC', 'W4' = '#9685AD', 'W9' = '#785BA5')
Group.color <- setNames(c('#EDBF36', '#7CCCE7'), c('SPF', 'GF'))
Tissue.color <- setNames(c('#F3CD3E', '#6F93CC', '#8A9095', '#CEAC80', '#ED7248',
                           '#CFDD9F', '#2277B3', '#E25B3D', '#437034', '#71C6E8'),
                         c('bm', 'colon', 'heart', 'kidney', 'liver', 'lung',
                           'si', 'skin', 'spleen', 'thymus'))

setwd('lym')  # v9 root / lym

seu <- readRDS('data/seurat/v9_lym_mc2seu.rds')
load('data/lym_annot_color.Rda')  # col.21

# =============================================================================
# helper functions
# =============================================================================

## z-score per gene (row)
zscore_mat <- function(mat){
  Mean <- rowMeans(mat)
  Sd <- apply(mat, 1, sd)
  mat <- sweep(mat, 1, Mean, FUN = '-')
  mat <- sweep(mat, 1, Sd, FUN = '/')
  return(mat)
}

## minibulk (AverageExpression on the LogNormalized 'data' layer) + row z-score
## 'target' ident must be set beforehand
makePsuedoBulk <- function(seu, x){
  x <- base::intersect(x, rownames(GetAssayData(seu, layer = 'data')))
  mat <- AverageExpression(seu, features = x, slot = 'data') # genes as rows
  mat <- zscore_mat(mat$RNA)
  return(mat)
}

## minibulk without z-score (for boxplots)
makePsuedoBulk_noScale <- function(seu, x){
  x <- base::intersect(x, rownames(GetAssayData(seu, layer = 'data')))
  mat <- AverageExpression(seu, features = x, slot = 'data')
  return(data.frame(mat$RNA))
}

## remove bad features (ribosomal, mt, Rik, Gm, ...) from a gene list
rm.bad.features <- function(genelist, species = 'Mus'){
  if(species != 'Mus'){stop('Only mouse genes were required')}
  pst_nr_term <- c("^Igk", "^Rps", "^Rpl", "Rik[0-9]?$", "^AC[0-9]+\\.", "^AL[0-9]+\\.",
                   "^Gm[0-9]+", "^Snor", "-ps[0-9]?$", "-rs[0-9]?$", "\\.", '^mt-',
                   "^Dnaj", "^Hist", "^Hsp", "^Atp", "^Uqc", "^Mcm[0-9]", "^Smc[0-9]",
                   "^A[A-Z]", "^B[A-Z]", "^D[0-9]")
  pst_nr_genes <- foreach(i = pst_nr_term, .combine = c) %do% grep(i, genelist, v = T)
  genelist <- setdiff(genelist, pst_nr_genes)
  return(genelist)
}

## DEG cutoffs
LFC.cutoff <- 0.5
FDR.cutoff <- 0.05
baseMean.cutoff <- 5

spf_up <- function(data){ # SPF-upregulated genes
  data %>% filter(log2FoldChange >= LFC.cutoff & padj <= FDR.cutoff & baseMean >= baseMean.cutoff)
}
gf_up <- function(data){ # GF-upregulated genes
  data %>% filter(log2FoldChange <= -LFC.cutoff & padj <= FDR.cutoff & baseMean >= baseMean.cutoff)
}

## add SPF/GF/NS sign column to a DEG table (for volcano)
add_sign <- function(data, time_label, good_features){
  data %>%
    mutate(Time = time_label) %>%
    filter(!is.na(padj)) %>%
    tibble::rownames_to_column(var = "Gene") %>%
    filter(Gene %in% good_features) %>%
    mutate(sign = ifelse(((log2FoldChange > LFC.cutoff) & (padj <= FDR.cutoff)), 'SPF',
                  ifelse(((log2FoldChange < -LFC.cutoff) & (padj <= FDR.cutoff)), 'GF', 'NS')))
}

## faceted volcano plot (SPF up: yellow, GF up: blue)
plot_volcano <- function(merge_deg, out_pdf, width = 8, height = 5,
                         max.overlaps_spf = 15, max.overlaps_gf = 15){
  p <- ggplot() +
    geom_point(data = merge_deg %>% filter(sign == 'NS'),
               aes(x = log2FoldChange, y = -log10(padj)), shape = 19, size = 0.2, color = 'grey') +
    geom_point(data = merge_deg %>% filter(sign == 'SPF'),
               aes(x = log2FoldChange, y = -log10(padj)), shape = 19, size = 0.4, color = '#EBBE35') +
    geom_text_repel(data = merge_deg %>% filter(sign == 'SPF'),
                    aes(x = log2FoldChange, y = -log10(padj), label = Gene), size = 4,
                    max.overlaps = max.overlaps_spf) +
    geom_point(data = merge_deg %>% filter(sign == 'GF'),
               aes(x = log2FoldChange, y = -log10(padj)), shape = 19, size = 0.4, color = '#7CCAE4') +
    geom_text_repel(data = merge_deg %>% filter(sign == 'GF'),
                    aes(x = log2FoldChange, y = -log10(padj), label = Gene), size = 4,
                    max.overlaps = max.overlaps_gf) +
    facet_wrap(~Time) + xlab('Log2FoldChange') + ylab('-log10(adjusted P)') +
    theme(strip.text = element_text(size = 16))
  ggsave(out_pdf, p, width = width, height = height)
  return(p)
}

## sample annotation for a minibulk matrix, built from the object metadata
## (robust against column-name mangling by AverageExpression)
get_annotation_col <- function(seu_target, mat){
  annotation_col <- seu_target@meta.data %>%
    mutate(target = paste(Time, Tissue_mod_only, Group, Mouse, sep = ' ')) %>%
    distinct(target, Time, Tissue_mod_only, Group, Mouse) %>%
    rename(Tissue = Tissue_mod_only) %>%
    mutate(Tissue = factor(Tissue, levels = c('bm', 'spleen', 'si', 'colon', 'liver',
                                              'lung', 'heart', 'kidney', 'skin', 'thymus')),
           Group = factor(Group, levels = c('SPF', 'GF'))) %>%
    arrange(Tissue, Group, Time, Mouse)
  annotation_col <- annotation_col[match(colnames(mat), annotation_col$target), ]
  stopifnot(all(annotation_col$target == colnames(mat)))
  return(annotation_col)
}

## DEG minibulk heatmap (rows = SPF/GF DEGs, split by direction)
plot_deg_heatmap <- function(seu_target, spf.degs, gf.degs, out_pdf,
                             width, height, mark_genes = c(),
                             extra_genes = c()){
  mat <- makePsuedoBulk(seu_target, c(gf.degs, spf.degs, extra_genes))
  annotation_col <- get_annotation_col(seu_target, mat)
  mat_plot <- as.matrix(mat[, annotation_col$target])

  top_anno <- columnAnnotation(
    Tissue = annotation_col$Tissue,
    Time = annotation_col$Time,
    Group = annotation_col$Group,
    col = list(Tissue = Tissue.color, Time = Time.color, Group = Group.color))

  gene_pos <- which(rownames(mat_plot) %in% mark_genes)
  row_anno <- rowAnnotation(gene = anno_mark(at = gene_pos,
                                             labels = rownames(mat_plot)[gene_pos]))

  gene_category <- data.frame(
    gene = rownames(mat_plot),
    category = factor(ifelse(rownames(mat_plot) %in% spf.degs, 'SPF', 'GF'),
                      levels = c('SPF', 'GF')))

  p <- Heatmap(mat_plot,
               col = circlize::colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
               cluster_rows = TRUE, cluster_columns = FALSE,
               show_column_names = FALSE, show_row_names = TRUE,
               row_names_gp = gpar(fontsize = 10),
               right_annotation = row_anno,
               top_annotation = top_anno,
               column_split = annotation_col$Tissue,
               row_split = gene_category$category,
               column_names_gp = gpar(fontsize = 4),
               column_title_rot = 45)
  pdf(out_pdf, width = width, height = height)
  draw(p)
  dev.off()
  return(p)
}

## generate per-cutoff cell lists: for each cutoff in (10, 15, 20), keep cells of
## minibulk samples (Time x Tissue x Group x Mouse) with >= cutoff cells
filter_cutoff_cells <- function(seu_target){
  target_cell_list <- list()
  for (i in c(10, 15, 20)){
    target_cutoff <- data.frame(table(seu_target$target)) %>%
      dplyr::rename(target = Var1) %>%
      tidyr::separate(col = 'target', into = c('Time', 'Tissue', 'Group', 'Mouse'), sep = ' ') %>%
      filter(Freq >= i) %>%
      mutate(target = paste(Time, Tissue, Group, Mouse, sep = ' ')) %>%
      pull(target)
    seu_target_filt <- subset(seu_target, target %in% target_cutoff)
    message('cutoff', i, ': ', dim(seu_target)[2], ' -> ', dim(seu_target_filt)[2], ' cells')
    target_cell_list[[sprintf('cutoff%s', i)]] <- colnames(seu_target_filt)
  }
  return(target_cell_list)
}

make_target <- function(seu_target){
  seu_target@meta.data <- seu_target@meta.data %>%
    dplyr::mutate(target = paste(Time, Tissue_mod_only, Group, Mouse, sep = ' '))
  Idents(seu_target) <- 'target'
  return(seu_target)
}


# =============================================================================
# cutoff cell lists (input of the Figure5.py pseudobulk DEG analysis)
# =============================================================================

## T_naive (T_naive_1 + T_naive_2; bm/si/skin excluded)
seu_target <- subset(seu, (annots %in% c('T_naive_1', 'T_naive_2')) &
                       Tissue_mod_only %out% c('bm', 'si', 'skin'))
seu_target <- make_target(seu_target)
target_cell_list <- filter_cutoff_cells(seu_target)
saveRDS(target_cell_list, 'output/seurat/T_naive_subset_cutoff_cells.rds')
T_naive_cutoff15 <- target_cell_list[['cutoff15']]

## B_mature (B_mature + B_mature_B2m)
seu_target <- subset(seu, (annots %in% c('B_mature_B2m', 'B_mature')))
seu_target <- make_target(seu_target)
target_cell_list <- filter_cutoff_cells(seu_target)
saveRDS(target_cell_list, 'output/seurat/B_mature_subset_cutoff_cells.rds')
B_mature_cutoff15 <- target_cell_list[['cutoff15']]

## CD8T_cyto
seu_target <- subset(seu, (annots %in% c('CD8T_cyto')))
seu_target <- make_target(seu_target)
target_cell_list <- filter_cutoff_cells(seu_target)
saveRDS(target_cell_list, 'output/seurat/CD8T_cyto_subset_cutoff_cells.rds')
CD8T_cyto_cutoff15 <- target_cell_list[['cutoff15']]

## CD8T_resi_cyto
seu_target <- subset(seu, (annots %in% c('CD8T_resi_cyto')))
seu_target <- make_target(seu_target)
target_cell_list <- filter_cutoff_cells(seu_target)
saveRDS(target_cell_list, 'output/seurat/CD8T_resi_cyto_subset_cutoff_cells.rds')


# =============================================================================
## Fig. 5a — DEG stat scatter, 8 lymphoid subsets x 3 time points
# =============================================================================

## genes kept after bad-feature removal
good_features <- rm.bad.features(rownames(seu))

time_list <- c('P16', 'W4', 'W9')
celltype_list <- c('B_pre', 'B_mature', 'DPT_Rag1', 'DPT_Itm2a',
                   'T_naive_2', 'CD8T_cyto', 'Treg', 'NK')

## collect DEG tables (generated by Figure5.py, 300-UMI downsampled pseudobulk)
deg_total <- data.frame()
for (celltype in celltype_list){
  for (time in time_list){
    target <- read.csv(sprintf('output/scanpy/decoupler/%s_%s_pseudobulk_deg_downsample300.csv',
                               celltype, time)) %>%
      rename(gene = X) %>%
      mutate(Cluster = celltype, Time = time)
    deg_total <- rbind(deg_total, target)
  }
}

Cluster_time.level <- c('B_pre_P16', 'B_pre_W4', 'B_pre_W9',
                        'B_mature_P16', 'B_mature_W4', 'B_mature_W9',
                        'DPT_Rag1_P16', 'DPT_Rag1_W4', 'DPT_Rag1_W9',
                        'DPT_Itm2a_P16', 'DPT_Itm2a_W4', 'DPT_Itm2a_W9',
                        'Treg_P16', 'Treg_W4', 'Treg_W9',
                        'T_naive_2_P16', 'T_naive_2_W4', 'T_naive_2_W9',
                        'CD8T_cyto_P16', 'CD8T_cyto_W4', 'CD8T_cyto_W9',
                        'NK_P16', 'NK_W4', 'NK_W9')

deg_total <- deg_total %>%
  filter(!is.na(padj)) %>%
  mutate(sig = if_else(
    (log2FoldChange >= LFC.cutoff) & (padj <= FDR.cutoff) & (baseMean >= baseMean.cutoff) & (gene %in% good_features), 'SPF',
    if_else((log2FoldChange <= -LFC.cutoff) & (padj <= FDR.cutoff) & (baseMean >= baseMean.cutoff) & (gene %in% good_features), 'GF', 'NS'))) %>%
  mutate(Cluster_Time = factor(paste(Cluster, Time, sep = '_'), levels = Cluster_time.level))

## top annotation bars (celltype + time)
property <- data.frame(celltype = rep(celltype_list, each = 3),
                       time = rep(time_list, 8)) %>%
  mutate(sample = paste0(celltype, '_', time))
celltype_bar <- ggplot(property, aes(sample, 1, fill = celltype)) +
  geom_tile() + scale_fill_manual(values = col.21) + theme_void()
time_bar <- ggplot(property, aes(sample, 1, fill = time)) +
  geom_tile() +
  scale_fill_manual(values = c('P16' = '#F3C5FF', 'W4' = '#9B89B3', 'W9' = '#845EC2')) +
  theme_void()

## stat scatter
deg_point <- deg_total %>%
  ggplot(aes(Cluster_Time, stat, color = sig, size = sig)) +
  geom_jitter(shape = 19) +
  scale_color_manual(values = c('SPF' = '#EBBE35', 'GF' = '#7CCAE4', 'NS' = 'grey')) +
  scale_size_manual(values = c('SPF' = 2, 'GF' = 2, 'NS' = 0.1)) +
  xlab('') +
  theme(axis.title.y = element_text(size = 8), axis.text.x = element_blank())

deg_point %>%
  insert_top(time_bar, height = 0.1) %>%
  insert_top(celltype_bar, height = 0.1)
ggsave('figures/Batch_deg.pdf', width = 10, height = 6)


# =============================================================================
## Fig. 5b — B_mature volcano (DEG tables from Figure5.py, cutoff15)
# =============================================================================
deg_P16 <- read.csv('output/scanpy/decoupler/B_mature_P16_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W4 <- read.csv('output/scanpy/decoupler/B_mature_W4_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W9 <- read.csv('output/scanpy/decoupler/B_mature_W9_pseudobulk_cutoff15_deg.csv', row.names = 1)

merge_deg <- rbind(add_sign(deg_P16, 'P16', good_features),
                   add_sign(deg_W4, 'W4', good_features),
                   add_sign(deg_W9, 'W9', good_features))
plot_volcano(merge_deg, 'figures/Volcano/B_mature_volcano.pdf',
             width = 8, height = 5, max.overlaps_gf = 50)

# =============================================================================
## Fig. 5c — T_naive volcano
# =============================================================================
deg_P16 <- read.csv('output/scanpy/decoupler/T_naive_P16_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W4 <- read.csv('output/scanpy/decoupler/T_naive_W4_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W9 <- read.csv('output/scanpy/decoupler/T_naive_W9_pseudobulk_cutoff15_deg.csv', row.names = 1)

merge_deg <- rbind(add_sign(deg_P16, 'P16', good_features),
                   add_sign(deg_W4, 'W4', good_features),
                   add_sign(deg_W9, 'W9', good_features))
plot_volcano(merge_deg, 'figures/Volcano/T_naive_volcano.pdf', width = 8, height = 5)

# =============================================================================
## Fig. 5d — CD8T_cyto volcano (W4 / W9 only)
# =============================================================================
deg_W4 <- read.csv('output/scanpy/decoupler/CD8T_cyto_W4_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W9 <- read.csv('output/scanpy/decoupler/CD8T_cyto_W9_pseudobulk_cutoff15_deg.csv', row.names = 1)

merge_deg <- rbind(add_sign(deg_W4, 'W4', good_features),
                   add_sign(deg_W9, 'W9', good_features))
plot_volcano(merge_deg, 'figures/Volcano/CD8T_cyto_volcano.pdf', width = 8, height = 5)

# =============================================================================
## related views — DPT_Rag1 / DPT_Itm2a volcano (thymus, no cutoff15)
# =============================================================================
for (ct in c('DPT_Rag1', 'DPT_Itm2a')){
  deg_P16 <- read.csv(sprintf('output/scanpy/decoupler/%s_P16_pseudobulk_deg.csv', ct), row.names = 1)
  deg_W4 <- read.csv(sprintf('output/scanpy/decoupler/%s_W4_pseudobulk_deg.csv', ct), row.names = 1)
  deg_W9 <- read.csv(sprintf('output/scanpy/decoupler/%s_W9_pseudobulk_deg.csv', ct), row.names = 1)
  merge_deg <- rbind(add_sign(deg_P16, 'P16', good_features),
                     add_sign(deg_W4, 'W4', good_features),
                     add_sign(deg_W9, 'W9', good_features))
  plot_volcano(merge_deg, sprintf('figures/Volcano/%s_volcano.pdf', ct), width = 8, height = 5)
}


# =============================================================================
## Fig. 5e — paired boxplots of Pdcd4 / Eif4g2 / Eif4a2
#  each point = mean expression of one subset (Time x Group x celltype);
#  only subsets with both SPF and GF points are kept (group_n == 2);
#  paired t-test, point shape marks the B / T / ILC lineage
# =============================================================================
target_gl <- c('Pdcd4', 'Eif4g2', 'Eif4a2')

seu@meta.data <- seu@meta.data %>%
  dplyr::mutate(target = paste(Time, Group, annots, sep = ' '))
Idents(seu) <- 'target'

mat <- makePsuedoBulk_noScale(seu, target_gl)

## sample annotation from metadata (robust against column-name mangling)
annotation_col <- seu@meta.data %>%
  mutate(target = paste(Time, Group, annots, sep = ' ')) %>%
  distinct(target, Time, Group, celltype = annots) %>%
  mutate(Group = factor(Group, levels = c('SPF', 'GF')))
annotation_col <- annotation_col[match(colnames(mat), annotation_col$target), ]
stopifnot(all(annotation_col$target == colnames(mat)))

celltype2lineage <- c('CD8T_cyto' = 'T', 'NK' = 'ILC', 'Treg' = 'T', 'T_naive_1' = 'T',
                      'B_mature_B2m' = 'B', 'ILC' = 'ILC', 'CD8T_cycling' = 'T',
                      'CD8T_resi_cyto' = 'T', 'B_mature' = 'B', 'DPT_Itm2a' = 'T',
                      'DPT_cycling' = 'T', 'DPT_Rag1' = 'T', 'DNT_Il2ra' = 'T',
                      'B_pre' = 'B', 'T_naive_transitory' = 'T', 'T_naive_2' = 'T',
                      'B_pro_cycling' = 'B', 'B_pro' = 'B', 'Plasma' = 'B',
                      'T_gd_2' = 'T', 'T_gd_1' = 'T')

plot_paired <- function(gene_i){
  mat_plot <- mat %>%
    t() %>%
    cbind(annotation_col) %>%
    tidyr::pivot_longer(cols = all_of(target_gl), names_to = 'gene', values_to = 'expression') %>%
    mutate(time_celltype = paste0(Time, ' ', celltype)) %>%
    dplyr::filter(gene == gene_i) %>%
    group_by(Time, celltype) %>%
    mutate(group_n = n()) %>%
    filter(group_n == 2) %>%
    mutate(lineage = celltype2lineage[celltype])

  p <- mat_plot %>%
    ggplot(aes(x = Group, y = expression)) +
    geom_boxplot(aes(fill = Group), show.legend = FALSE, width = 0.6, outliers = FALSE) +
    scale_fill_manual(values = c('#DBB33B', '#77BBD2')) +
    geom_jitter(aes(shape = lineage), size = 2, width = 0.1) +
    geom_line(aes(group = time_celltype), color = 'gray', lwd = 0.5) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(colour = 'black', linewidth = 1),
          panel.background = element_blank(),
          plot.title = element_text(size = 15, hjust = 0.5, face = 'italic'),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 15, color = 'black')) +
    labs(x = '', y = 'expression', title = gene_i) +
    scale_shape_manual(values = c('B' = 0, 'T' = 1, 'ILC' = 2)) +
    ggpubr::stat_compare_means(method = "t.test", paired = TRUE,
                               comparisons = list(c("SPF", "GF"))) +
    facet_wrap(~Time)
  ggsave(sprintf('figures/%s_group_paired_plot.pdf', gene_i), p, width = 4, height = 5)
  return(p)
}

for (gene_i in target_gl){
  plot_paired(gene_i)
}


# =============================================================================
## Fig. S5a — lymphoid marker heatmap
#  top metacell markers (choose_genes_from_mc) on merged-subset minibulk,
#  rows in the manually curated order gl_new_older
# =============================================================================

require(metacell)

figpath <- 'output/mc1'; dbpath <- 'data/mc1'
scdb_init(dbpath, force_reinit = TRUE)
scfigs_init(figpath)
id <- 'mouse_v9_lym'
mat_mc <- scdb_mat(id)
mc <- scdb_mc(paste0(id, '_filt'))

get_bad_features <- function(genelist, species = 'Mus'){
  if(species != 'Mus'){stop('Only mouse genes were required')}
  pst_nr_term <- c("^Igk", "^Rps", "^Rpl", "Rik[0-9]?$", "^AC[0-9]+\\.", "^AL[0-9]+\\.",
                   "^Gm[0-9]+", "^Snor", "-ps[0-9]?$", "-rs[0-9]?$", "\\.",
                   "^Dnaj", "^Hist", "^Hsp", "^Atp", "^Uqc", "^Mcm[0-9]", "^Smc[0-9]",
                   "^A[A-Z]", "^B[A-Z]", "^D[0-9]")
  pst_nr_genes <- foreach(i = pst_nr_term, .combine = c) %do% grep(i, genelist, v = T)
  return(pst_nr_genes)
}

## choose the top marker genes of each metacell from the footprint matrix
choose_genes_from_mc <- function(mc, mat, good_mcs = colnames(mc@mc_fp),
                                 nms_per_mc = 5, nms_thresh = 5, max_num = Inf,
                                 bad_genes = c(), must_haves = c(), ord = "none"){
  lfp <- log2(mc@mc_fp[, good_mcs])
  nms <- unique(as.vector(apply(lfp, 2, function(x){ head(names(sort(x, decreasing = T)), nms_per_mc) })))
  nms <- setdiff(nms, c(bad_genes, names(which(apply(lfp[nms, ], 1, max) < nms_thresh))))
  nms <- union(must_haves, head(names(sort(apply(lfp[nms, ], 1, max), T)), max_num - length(must_haves)))
  if (ord == "hc") {
    nms <- nms[hclust(dist(cor(t(lfp[nms, ]))), "ward.D2")$order]
  } else if (ord == "max.col") {
    nms <- nms[order(max.col(lfp[nms, ]), rowSums(as.matrix(mat@mat[nms, ])))]
  }
  nms
}

bad_features <- get_bad_features(mat_mc@genes)
genes_selected <- choose_genes_from_mc(mc, mat_mc,
                                       bad_genes = bad_features,
                                       nms_per_mc = 10,
                                       nms_thresh = 2)

## merged subsets: B_mature(+B2m), T_naive(1/2/transitory)
seu@meta.data <- seu@meta.data %>% mutate(annots2 = annots) %>%
  mutate(annots2 = if_else(annots %in% c('B_mature', 'B_mature_B2m'), 'B_mature', annots2)) %>%
  mutate(annots2 = if_else(annots %in% c('T_naive_1', 'T_naive_2', 'T_naive_transitory'), 'T_naive', annots2))

## minibulk unit: annots2 x Group x Time x Mouse
seu@meta.data <- seu@meta.data %>%
  dplyr::mutate(target = paste(annots2, Group, Time, Mouse, sep = ' '))
Idents(seu) <- 'target'

celltype.levels <- c('B_pro', 'B_pro_cycling', 'B_pre',
                     'B_mature', 'Plasma',
                     'DNT_Il2ra', 'DPT_Rag1', 'DPT_cycling', 'DPT_Itm2a',
                     'T_naive',
                     'CD8T_cycling', 'CD8T_cyto', 'CD8T_resi_cyto',
                     'Treg', 'T_gd_1', 'T_gd_2', 'ILC', 'NK')

## keep minibulk samples with >= 10 cells
min.cell.num.cutoff <- 10
target_cutoff <- data.frame(table(seu$target)) %>%
  dplyr::rename(target = Var1) %>%
  tidyr::separate(col = 'target', into = c('celltype', 'Group', 'Time', 'Mouse'), sep = ' ') %>%
  dplyr::mutate(celltype = str_replace_all(celltype, "-", "_"),
                celltype = factor(celltype, levels = celltype.levels),
                Group = factor(Group, levels = c('SPF', 'GF'))) %>%
  dplyr::filter(Freq >= min.cell.num.cutoff) %>%
  mutate(target = paste(celltype, Group, Time, Mouse, sep = ' ')) %>%
  pull(target)
seu_target <- subset(seu, target %in% target_cutoff)
## remove Plasma and the 3 cycling subsets
seu_target <- subset(seu_target, annots %out% c('Plasma', 'B_pro_cycling', 'CD8T_cycling', 'DPT_cycling'))

mat <- makePsuedoBulk(seu_target, c(genes_selected, 'Itga4', 'Klrg1'))

annotation_col <- seu_target@meta.data %>%
  mutate(target = paste(annots2, Group, Time, Mouse, sep = ' ')) %>%
  distinct(target, celltype = annots2, Group, Time, Mouse) %>%
  mutate(celltype = factor(celltype, levels = celltype.levels),
         Group = factor(Group, levels = c('SPF', 'GF'))) %>%
  arrange(celltype, Group, Time, Mouse)
annotation_col <- annotation_col[match(colnames(mat), annotation_col$target), ]
stopifnot(all(annotation_col$target == colnames(mat)))

mat_plot <- as.matrix(mat[, annotation_col$target])

## column annotation
celltype_colors <- c(
  B_pro = '#7399a0', B_pre = '#006d42', B_mature = '#00417d',
  DNT_Il2ra = '#e1352a', DPT_Rag1 = '#8b5f65', DPT_Itm2a = '#b22524', T_naive = '#FF00FF',
  CD8T_cyto = '#f3762b', CD8T_resi_cyto = '#956e9e',
  Treg = '#ffad30', T_gd_1 = '#fa7d56', T_gd_2 = '#003e59',
  ILC = '#e67f98', NK = '#0065a9')
group_colors <- c("GF" = "#7BC3DC", "SPF" = "#E4B93A")
time_colors <- c('P16' = '#F3C5FF', 'W4' = '#9B89B3', 'W9' = '#845EC2')
top_anno <- columnAnnotation(
  group = annotation_col$Group,
  time = annotation_col$Time,
  col = list(group = group_colors, time = time_colors))

## manually curated gene order
gl_new_older <- c('Ccnb2','Rrm2','Plac8','Ptma','Cd24a','Pou2af1','Ebf1','Blnk','Mzb1','Tifa','Vpreb3','Smarca4','Chchd10','Pafah1b3','Cd79a','Adgrg1','Igll1','Vpreb1','Selenop','Ly6d','Iglj2','Iglj3','Iglj1','Iglv1','Cd79b','Ighm','Syk','Fcrla','Spib','Apoe','Cd83','H2-DMb2','H2-Eb1','H2-Ab1','H2-Aa','Fcmr','Ighd','Mir5107','Cd74','Txndc5','Tnfrsf13c','Ifi30','Scd1','Bank1','Ms4a1','Lars2','H2afx','Xbp1','Stmn1','Il2ra','Uck2','Cdk6','Fkbp5','Dut','Dntt','Myb','Themis','Arpp21','Aqp11','Bcl11b','Cd8a','Ldhb','Rmnd5a','Ets2','Cd8b1','Nsg2','Cd4','Rag1','Trp53inp1','Ccr9','Satb1','Trbc1','Trbc2','Tcf7','Mir8113','Cd28','Trac','Itm2a','Rgs10','Cd5','Tubb5','Lef1','Mb21d1','Dgka','Igfbp4','Dapl1','Ms4a6b','Slfn1','Sell','Klf2','Ccr7','S1pr1','Ifit1','Ifi27l2a','Stat1','Tmsb10','Mif','Ran','Ms4a4b','Klrg1','Itga4','Ly6c2','H2-Q7','Gimap4','Hcst','Ube2c','Mki67','Top2a','Tcrg-C4','Tcrg-C2','Car2','Cd160','Gzma','Ccl5','Gzmb','Jchain','Igha','Ccna2','Crip1','Ctla4','Tnfrsf4','Icos','Ikzf2','Id2','Ahnak','Tnfaip3','Junb','Dusp1','Fgl2','Bhlhe40','Hmgb2','Itgb1','Lgals1','S100a10','Fxyd5','S100a11','S100a4','S100a6','Lmna','Vim','Ccnd2','Itgae','Ctsw','Cd7','Arl4c','Lgals3','Trdc','Tcrg-C1','Ctla2a','Ctla2b','Trat1','Cd3e','Cd3g','Rgcc','Thy1','Cd3d','Pcna','Rgs1','Fos','Crem','Nr4a3','Nr4a2','Fosb','Rgs2','Klf4','Gem','Areg','Egr1','Tuba1b','Gzmc','Nr4a1','Nfkbiz','Nfkbia','Gadd45b','Ccr8','Odc1','Rora','Gata3','Emb','Ccr2','Maf','Arg1','Hilpda','Cxcr6','Il22','Tmem176a','Tmem176b','Il7r','Klrb1b','Prf1','Klrk1','Nkg7','Klrd1','Fcer1g','Xcl1','Il2rb','Klra7','Ccl3','Tyrobp','Ccl4','Klre1','Klrb1c','Ncr1','Klra4')
mat_plot <- mat_plot[intersect(gl_new_older, rownames(mat_plot)), ]

## genes highlighted on the left
gl <- c('Cd24a','Mzb1','Vpreb3','Cd79a','Igll1','Vpreb1','Ly6d','Iglj2','Ighm','Syk','H2-DMb2','H2-Eb1','H2-Ab1','H2-Aa','Fcmr','Ighd','Cd74',
        'Stmn1','Il2ra','Dntt','Themis','Bcl11b','Cd8a','Cd8b1','Cd4','Rag1','Ccr9','Tcf7','Itm2a','Cd5','Lef1','Sell','Klf2','Ccr7','S1pr1',
        'Stat1','Itga4','Ly6c2','Mki67','Top2a','Gzma','Ccl5','Gzmb','Ctla4','Tnfrsf4','Ikzf2','Id2','Itgae','Trdc','Tcrg-C1','Cd3e','Thy1','Nr4a3','Nr4a2','Nr4a1','Klf4','Areg','Egr1','Gzmc',
        'Nfkbiz','Cxcr6','Rora','Gata3','Il22','Il7r','Prf1','Nkg7','Xcl1','Ncr1')
gene_pos <- which(rownames(mat_plot) %in% gl)
row_anno <- rowAnnotation(gene = anno_mark(at = gene_pos,
                                           labels = rownames(mat_plot)[gene_pos]))

p <- Heatmap(mat_plot,
             col = circlize::colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')),
             cluster_rows = FALSE, cluster_columns = FALSE,
             show_column_names = FALSE, show_row_names = FALSE,
             left_annotation = row_anno,
             top_annotation = top_anno,
             column_split = annotation_col$celltype,
             column_title_rot = 15,
             name = 'Row\nZ-score')
pdf("./output/mc1/heatmap/Lym_heatmap_choose_top10markers_202g_minibulk_min_size_10_merged_ordered.pdf",
    width = 18, height = 12)
draw(p)
dev.off()


# =============================================================================
## Fig. S5b — B_mature DEG heatmap (W9, cutoff15 cells)
#  DEG cutoffs for the heatmap: |LFC| >= 1, FDR <= 0.05, baseMean >= 5
# =============================================================================
deg_P16 <- read.csv('output/scanpy/decoupler/B_mature_P16_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W4 <- read.csv('output/scanpy/decoupler/B_mature_W4_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W9 <- read.csv('output/scanpy/decoupler/B_mature_W9_pseudobulk_cutoff15_deg.csv', row.names = 1)

LFC.cutoff <- 1
W9_deg_spf <- rm.bad.features(rownames(spf_up(deg_W9)))
W9_deg_gf <- rm.bad.features(rownames(gf_up(deg_W9)))
spf.degs <- W9_deg_spf
gf.degs <- W9_deg_gf
LFC.cutoff <- 0.5  # restore default

## W9 cells; skin/si/kidney/lung/colon excluded (few replicates)
seu_target <- subset(seu, cellID %in% B_mature_cutoff15 &
                       Time == 'W9' &
                       Tissue_mod_only %out% c('skin', 'si', 'kidney', 'lung', 'colon'))
seu_target <- NormalizeData(seu_target, normalization.method = "LogNormalize", scale.factor = 10000)
seu_target <- make_target(seu_target)

plot_deg_heatmap(seu_target, spf.degs, gf.degs,
                 out_pdf = 'figures/heatmap/B_mature_heatmap.pdf',
                 width = 10, height = 16,
                 mark_genes = c('Ms4a1','Pax5','Bank1','Notch2','Rel','Pten','Foxn3','S1pr1','Cd24a'))

# =============================================================================
## Fig. S5c — T_naive DEG heatmap (W9, cutoff15 cells)
# =============================================================================
deg_P16 <- read.csv('output/scanpy/decoupler/T_naive_P16_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W4 <- read.csv('output/scanpy/decoupler/T_naive_W4_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W9 <- read.csv('output/scanpy/decoupler/T_naive_W9_pseudobulk_cutoff15_deg.csv', row.names = 1)

LFC.cutoff <- 1
W9_deg_spf <- rm.bad.features(rownames(spf_up(deg_W9)))
W9_deg_gf <- rm.bad.features(rownames(gf_up(deg_W9)))
spf.degs <- W9_deg_spf
gf.degs <- W9_deg_gf
LFC.cutoff <- 0.5

## W9 cells; spleen/liver/heart only
seu_target <- subset(seu, cellID %in% T_naive_cutoff15 &
                       Time == 'W9' &
                       Tissue_mod_only %in% c('spleen', 'liver', 'heart'))
seu_target <- NormalizeData(seu_target, normalization.method = "LogNormalize", scale.factor = 10000)
seu_target <- make_target(seu_target)

plot_deg_heatmap(seu_target, spf.degs, gf.degs,
                 out_pdf = 'figures/heatmap/T_naive_heatmap.pdf',
                 width = 6, height = 5,
                 mark_genes = c('Sell', 'Tcf7', 'Lef1', 'Ly6e', 'Cd52'))

# =============================================================================
## Fig. S5d — CD8T_cyto DEG heatmap (W4 + W9, cutoff15 cells)
# =============================================================================
deg_W4 <- read.csv('output/scanpy/decoupler/CD8T_cyto_W4_pseudobulk_cutoff15_deg.csv', row.names = 1)
deg_W9 <- read.csv('output/scanpy/decoupler/CD8T_cyto_W9_pseudobulk_cutoff15_deg.csv', row.names = 1)

LFC.cutoff <- 1
W4_deg_spf <- rm.bad.features(rownames(spf_up(deg_W4)))
W4_deg_gf <- rm.bad.features(rownames(gf_up(deg_W4)))
W9_deg_spf <- rm.bad.features(rownames(spf_up(deg_W9)))
W9_deg_gf <- rm.bad.features(rownames(gf_up(deg_W9)))
spf.degs <- c(W4_deg_spf, W9_deg_spf)
gf.degs <- c(W4_deg_gf, W9_deg_gf)
LFC.cutoff <- 0.5

## W4/W9 cells; bm/colon/heart/skin excluded
seu_target <- subset(seu, cellID %in% CD8T_cyto_cutoff15 &
                       Time != 'P16' &
                       Tissue_mod_only %out% c('bm', 'colon', 'heart', 'skin'))
seu_target <- NormalizeData(seu_target, normalization.method = "LogNormalize", scale.factor = 10000)
seu_target <- make_target(seu_target)

plot_deg_heatmap(seu_target, spf.degs, gf.degs,
                 out_pdf = 'figures/heatmap/CD8T_cyto_heatmap.pdf',
                 width = 12, height = 5)

# =============================================================================
## related views — DPT_Rag1 / DPT_Itm2a DEG heatmap (thymus, W9)
# =============================================================================
for (ct in c('DPT_Rag1', 'DPT_Itm2a')){
  deg_P16 <- read.csv(sprintf('output/scanpy/decoupler/%s_P16_pseudobulk_deg.csv', ct), row.names = 1)
  deg_W4 <- read.csv(sprintf('output/scanpy/decoupler/%s_W4_pseudobulk_deg.csv', ct), row.names = 1)
  deg_W9 <- read.csv(sprintf('output/scanpy/decoupler/%s_W9_pseudobulk_deg.csv', ct), row.names = 1)

  spf.degs <- rm.bad.features(c(rownames(spf_up(deg_P16)), rownames(spf_up(deg_W4)), rownames(spf_up(deg_W9))))
  gf.degs <- rm.bad.features(c(rownames(gf_up(deg_P16)), rownames(gf_up(deg_W4)), rownames(gf_up(deg_W9))))

  seu_target <- subset(seu, (annots %in% ct) & (Tissue_mod_only == 'thymus') & (Time == 'W9'))
  seu_target <- NormalizeData(seu_target, normalization.method = "LogNormalize", scale.factor = 10000)
  seu_target <- make_target(seu_target)

  plot_deg_heatmap(seu_target, spf.degs, gf.degs,
                   out_pdf = sprintf('figures/heatmap/%s_heatmap.pdf', ct),
                   width = 4.5, height = ifelse(ct == 'DPT_Rag1', 10, 4),
                   extra_genes = 'Pdcd4')
}


# =============================================================================
## Fig. S5e — HALLMARK GSEA dotplot (B_mature / CD8T_cyto / T_naive)
#  input: gsea_focus_3clusters.tsv — GSEA (desktop) output table of the three
#  subsets, columns NAME / NES / FDR.q.val / Cluster
# =============================================================================
require(tidyplots)

res <- read.csv('output/scanpy/gsea_focus_3clusters.tsv', sep = '\t')

## SPF-enriched (NES > 0) pathways
res %>% filter(FDR.q.val <= 0.25) %>%
  filter(NES > 0 & NAME %in% c('HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                               'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_TGF_BETA_SIGNALING',
                               'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_IL2_STAT5_SIGNALING') |
           NAME %in% c('HALLMARK_PROTEIN_SECRETION', 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE')) %>%
  tidyplot(y = NAME, x = NES, color = Cluster, size = -log10(FDR.q.val)) |>
  add_data_points() |>
  sort_y_axis_labels() |>
  adjust_size(width = 60, height = 80) |>
  remove_y_axis_title() |>
  adjust_font(fontsize = 12) |>
  save_plot('figures/GSEA_3clusters_SPF_points.pdf')

## GF-enriched (NES < 0) pathways
res %>% filter(FDR.q.val <= 0.25) %>%
  filter(NES < 0) %>%
  tidyplot(y = NAME, x = NES, color = Cluster, size = -log10(FDR.q.val)) |>
  add_data_points() |>
  sort_y_axis_labels(.reverse = TRUE) |>
  adjust_size(width = 60, height = 120) |>
  save_plot('figures/GSEA_3clusters_GF_points.pdf')


# =============================================================================
## Fig. S6a — lymphoid metacells: SPF fraction vs pathway score
#  pathway score = mean metacell footprint over the GO-term gene set
#  (metacell 183 excluded — annotated as 'other')
# =============================================================================

## SPF fraction per metacell
mc_info <- read.csv('output/seurat/lym_88180c_mc_info.csv', sep = ' ')
mc_info %>%
  dplyr::select(mc_id, Group) %>%
  group_by(mc_id) %>%
  count(Group) %>%
  mutate(SPF = n / sum(n)) %>%
  filter(Group == 'SPF') -> df_spf_freq

require(metacell)
scdb_init(base_dir = 'data/mc1', force_reinit = TRUE)
mc <- scdb_mc('mouse_v9_lym_filt')

## metacell footprint score -> lm fit against SPF fraction
## drop_mc: metacell columns to exclude (lymphoid: mc 183, annotated as 'other')
plot_score_lm <- function(gl, score_name, mc, df_spf_freq, out_dir = 'output/mc1/score_dot_lmfit',
                          drop_mc = 183){
  fp <- mc@mc_fp
  if (!is.null(drop_mc)) fp <- fp[, -drop_mc]
  df_plot <- data.frame(df_spf_freq,
                        pathway_score = Matrix::colSums(
                          fp[intersect(rownames(fp), gl), ]) / length(gl))
  n <- dim(df_plot)[1]
  cor_res <- cor(df_plot$pathway_score, df_plot$SPF)
  lm_res <- summary(lm(df_plot$pathway_score ~ df_plot$SPF))
  r2 <- lm_res$r.squared
  pvalue <- lm_res$coefficients[2, 4]
  p <- df_plot %>%
    ggplot(aes(SPF, pathway_score, size = n)) +
    geom_point(shape = 21) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Cor = ", round(cor_res, 3), "\nR2 = ", round(r2, 3),
                            "\nP = ", round(pvalue, 50), "\nn = ", n),
             hjust = 1.1, vjust = 2, color = "blue", size = 4) +
    geom_smooth(method = lm, color = "red", fill = "#69b3a2", se = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, sprintf('%s_lm_fit.pdf', score_name)), p, width = 8, height = 6)
  return(p)
}

gl_oxy <- read.csv('ref/go_term/GO_term_0006119_oxidative_phosphorylation.txt', sep = '\t')[, 1]
gl_tca <- read.csv('ref/go_term/GO_term_0006099_tricarboxylic_acid_cycle.txt', sep = '\t')[, 1]
gl_gly <- read.csv('ref/go_term/GO_term_0006096_glycolytic_process.txt', sep = '\t')[, 1]

plot_score_lm(gl_oxy, 'OXPHOS', mc, df_spf_freq)
plot_score_lm(gl_tca, 'TCA', mc, df_spf_freq)
plot_score_lm(gl_gly, 'Glycolysis', mc, df_spf_freq)


# =============================================================================
## Fig. S6b — myeloid metacells: SPF fraction vs pathway score
#  NOTE: the original code of this panel was not found in the notebooks;
#  it is reconstructed here following exactly the same procedure as Fig. S6a,
#  with the myeloid metacell database (mye/data/mc1, id 'mouse_v9_mye').
# =============================================================================

mc_info_mye <- read.csv('../mye/output/seurat/mye_32469c_mc_info.csv', sep = ' ')
mc_info_mye %>%
  dplyr::select(mc_id, Group) %>%
  group_by(mc_id) %>%
  count(Group) %>%
  mutate(SPF = n / sum(n)) %>%
  filter(Group == 'SPF') -> df_spf_freq_mye

scdb_init(base_dir = '../mye/data/mc1', force_reinit = TRUE)
mc_mye <- scdb_mc('mouse_v9_mye_filt')

plot_score_lm(gl_oxy, 'OXPHOS_mye', mc_mye, df_spf_freq_mye, out_dir = '../mye/output/mc1/score_dot_lmfit', drop_mc = NULL)
plot_score_lm(gl_tca, 'TCA_mye', mc_mye, df_spf_freq_mye, out_dir = '../mye/output/mc1/score_dot_lmfit', drop_mc = NULL)
plot_score_lm(gl_gly, 'Glycolysis_mye', mc_mye, df_spf_freq_mye, out_dir = '../mye/output/mc1/score_dot_lmfit', drop_mc = NULL)
