###############################################################################
# Fig. S6c — Pdcd4 expression in public scRNA-seq of SFB-colonized mouse
#            CD4 T cells (tumor and small-intestine lamina propria, siLP)
#
# Covered panel:
#   Fig. S6c: Pdcd4 violin plots per annotated CD4 T subset, split by condition
#             (SFB transfer experiment), wilcox_test with FDR adjustment
#
# Input:
#   SeuratObj_CD4_Tumor.rds  (public data re-analysis object, CD4 T, tumor)
#   SeuratObj_CD4_SILP.rds   (public data re-analysis object, CD4 T, siLP)
#
# Source notebook: sfb_tumor_siLP.ipynb
###############################################################################

require(Seurat)
require(dplyr)
require(ggplot2)
require(ggpubr)

seu <- readRDS('SeuratObj_CD4_Tumor.rds')    # tumor CD4 T cells
seu2 <- readRDS('SeuratObj_CD4_SILP.rds')    # siLP CD4 T cells

## QC view
UMAPPlot(seu)
UMAPPlot(seu2)
DotPlot(seu, group.by = 'anno',
        features = c('Pdcd4', 'Cd4', 'Cd3d', 'Cd44', 'Sell', 'Tcf7', 'Cd8a', 'Gzmb', 'Ifng'))
DotPlot(seu2, group.by = 'anno',
        features = c('Pdcd4', 'Cd4', 'Cd3d', 'Cd44', 'Sell', 'Tcf7', 'Cd8a', 'Gzmb', 'Ifng'))

## subset order
anno.levels <- c('Activated CD4+ T cells', 'TfH', 'Memory-like CD4+ T cells', 'Tregs',
                 'Gata3+ Tregs', 'Prolif. CD4+ T cells', 'NKT cells',
                 'IFN-Stim. CD4+ T cells', 'Granzyme-Bhi Tregs', 'Th1-like')
seu$anno <- factor(seu$anno, level = anno.levels)

## Fig. S6c — Pdcd4 violin per subset, split by condition (orig.ident)
VlnPlot(seu, group.by = 'anno', features = c('Pdcd4'), assay = 'RNA',
        combine = TRUE, split.by = "orig.ident") +
  ggpubr::geom_pwc(aes(group = split),
                   method = "wilcox_test", tip.length = 0,
                   hide.ns = TRUE, label = "p.adj.format",
                   p.adjust.method = "fdr") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = '')
ggsave('sfb_tumor_Pdcd4_vln.pdf', width = 8, height = 4)

## same for siLP
VlnPlot(seu2, group.by = 'anno', features = c('Pdcd4'), assay = 'RNA',
        combine = TRUE, split.by = "orig.ident") +
  ggpubr::geom_pwc(aes(group = split),
                   method = "wilcox_test", tip.length = 0,
                   hide.ns = TRUE, label = "p.adj.format",
                   p.adjust.method = "fdr") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = '')
ggsave('sfb_siLP_Pdcd4_vln.pdf', width = 8, height = 4)
