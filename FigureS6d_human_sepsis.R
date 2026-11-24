###############################################################################
# Fig. S6d — human sepsis scRNA-seq (GSE279452): plotting (R part)
#
# Covered panel:
#   Fig. S6d: PDCD4 mean expression per sample across sepsis groups (adult and
#             pediatric), and sample-mean correlation of PDCD4 with effector
#             (IFNG/TNF/GZMB/PRF1) and naive (TCF7/SELL) markers in CD8 T cells
#
# Input:
#   human_sepsis_GSE279452_887650c_processed.rds
#     (Seurat object converted from the h5ad written by FigureS6d_human_sepsis.py;
#      CD8 T cells = leiden_8 cluster '0')
#
# Source notebook: lps_human.ipynb
###############################################################################

require(Seurat)
require(tidyplots)
require(tidyverse)
require(ggpubr)
`%out%` <- Negate(`%in%`)

seu <- readRDS('human_sepsis_GSE279452_887650c_processed.rds')

## QC views
UMAPPlot(seu, group.by = 'celltype')
UMAPPlot(seu, group.by = 'leiden_8')

## PDCD4 per celltype, split by group
VlnPlot(seu, group.by = 'celltype', features = c('PDCD4'), assay = 'RNA',
        slot = 'counts', combine = TRUE, split.by = "group") +
  ggpubr::geom_pwc(aes(group = split),
                   method = "wilcox_test", tip.length = 0,
                   hide.ns = TRUE, label = "p.adj.format",
                   p.adjust.method = "fdr") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = '')

# =============================================================================
# CD8 T cells (leiden_8 == '0'), pediatric + adult groups;
# samples with < 50 cells removed
# =============================================================================
seu_cd8t <- subset(seu, leiden_8 == '0' & group %in% c('PHC', 'PIC', 'Abd-PS', 'Res-PS',
                                                       'HC', 'IC', 'Abd-S', 'Res-S', 'Uri-S', 'Skin-S'))
seu_cd8t$group <- factor(seu_cd8t$group,
                         levels = c('PHC', 'PIC', 'Abd-PS', 'Res-PS',
                                    'HC', 'IC', 'Abd-S', 'Res-S', 'Uri-S', 'Skin-S'))

tab <- table(seu_cd8t$sample)
seu_cd8t <- subset(seu_cd8t, sample %in% names(tab)[tab >= 50])

## sample-mean expression of PDCD4 and markers
df <- data.frame(PDCD4 = as.matrix(GetAssayData(seu_cd8t, layer = 'counts'))['PDCD4', ],
                 IFNG  = as.matrix(GetAssayData(seu_cd8t, layer = 'counts'))['IFNG', ],
                 TNF   = as.matrix(GetAssayData(seu_cd8t, layer = 'counts'))['TNF', ],
                 GZMB  = as.matrix(GetAssayData(seu_cd8t, layer = 'counts'))['GZMB', ],
                 PRF1  = as.matrix(GetAssayData(seu_cd8t, layer = 'counts'))['PRF1', ],
                 TCF7  = as.matrix(GetAssayData(seu_cd8t, layer = 'counts'))['TCF7', ],
                 SELL  = as.matrix(GetAssayData(seu_cd8t, layer = 'counts'))['SELL', ],
                 group = seu_cd8t$group,
                 sample = seu_cd8t$sample) |>
  group_by(group, sample) |>
  mutate(PDCD4_mean = mean(PDCD4), TCF7_mean = mean(TCF7), SELL_mean = mean(SELL),
         IFNG_mean = mean(IFNG), TNF_mean = mean(TNF), GZMB_mean = mean(GZMB),
         PRF1_mean = mean(PRF1)) |>
  distinct(group, sample, PDCD4_mean, TCF7_mean, SELL_mean, IFNG_mean, TNF_mean,
           PRF1_mean, GZMB_mean) |>
  pivot_longer(cols = c(TCF7_mean, SELL_mean, IFNG_mean, TNF_mean, PRF1_mean, GZMB_mean),
               names_to = 'gene', values_to = 'mean') |>
  mutate(gene = factor(gene, levels = c('TCF7_mean', 'SELL_mean', 'IFNG_mean',
                                        'TNF_mean', 'PRF1_mean', 'GZMB_mean')))
write.csv(df, 'human_sepsis_seu_bulk_mean_pd4_df.csv')

## Fig. S6d — PDCD4 vs marker sample-mean correlation (spearman)
df |>
  tidyplot(PDCD4_mean, mean) |>
  add_data_points() |>
  adjust_size(width = 50, height = 40) |>
  add(geom_smooth(method = "lm", se = TRUE, color = "steelblue")) |>
  add(ggpubr::stat_cor(method = "spearman",
                       aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                       label.x.npc = "left", label.y.npc = "top",
                       size = 4, color = "black")) |>
  split_plot(by = gene, ncol = 6) |>
  save_plot('pdcd4_correlation.pdf')

# =============================================================================
# PDCD4 mean expression per group (adult + pediatric)
# =============================================================================
df2 <- df |>
  pivot_wider(names_from = gene, values_from = mean) |>
  pivot_longer(cols = ends_with('_mean'), names_to = 'gene', values_to = 'mean',
               names_pattern = "(.+)_.+") |>
  mutate(gene = factor(gene, levels = c('PDCD4', 'TCF7', 'SELL', 'IFNG', 'TNF', 'PRF1', 'GZMB'))) |>
  mutate(cat = case_when(group %in% c('PHC', 'PIC', 'Abd-PS', 'Res-PS') ~ 'Pediatric',
                         group %in% c('HC', 'IC', 'Abd-S', 'Res-S', 'Uri-S', 'Skin-S') ~ 'Adult'))

## all 10 groups
df2 |>
  ungroup() |>
  filter(!is.na(mean)) |>
  tidyplot(group, mean, color = cat) |>
  add_data_points_beeswarm(size = 0.1) |>
  add_mean_bar(alpha = 0.3) |>
  add_sem_errorbar() |>
  adjust_size(width = 50, height = 30) |>
  adjust_x_axis(rotate_labels = 90) |>
  add_test_pvalue(method = 'wilcox_test', p.adjust.method = 'BH',
                  comparisons = list(c(1, 3), c(1, 4), c(2, 3), c(2, 4)),
                  hide_info = TRUE) |>
  split_plot(by = gene, ncol = 4) |>
  save_plot('PDCD4_mean_10grps.pdf')

## PDCD4 only: pediatric vs adult group comparisons
df2 |>
  filter(gene == 'PDCD4') |>
  ungroup() |>
  tidyplot(group, mean, color = cat) |>
  add_data_points_beeswarm() |>
  add_mean_bar(alpha = 0.3) |>
  add_sem_errorbar() |>
  adjust_size(width = 50, height = 50) |>
  adjust_x_axis(rotate_labels = 90) |>
  add_test_pvalue(method = 'wilcox_test', p.adjust.method = 'BH', hide.ns = TRUE,
                  comparisons = list(c(1, 3), c(1, 4), c(2, 3), c(2, 4),
                                     c(5, 7), c(5, 8), c(5, 9), c(5, 10)),
                  hide_info = TRUE) |>
  save_plot('PDCD4_mean_10grps.pdf')
