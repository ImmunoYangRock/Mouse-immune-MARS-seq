#################################################################
## R code to create figures related to Figure4
## Figs 4a-4h and Figs s4a-s4f
## Description: compositional analysis of lymphoid subsets —
## metacell 2D projection and FeaturePlots, stacked-bar frequencies,
## temporal variation, GF-vs-SPF deviation, oversampled GF%
## projection, and per-subset frequency line plots.
#################################################################

require(dplyr)
require(tidyr)
require(tibble)
require(stringr)
require(ggplot2)
require(ggpubr)
require(Seurat)

## load data
seu <- readRDS('lym/data/seurat/v9_lym_mc2seu.rds') # lymphoid seurat object, metacell 2d projection stored in reduction 'metacell'
load('lym/data/lym_annot_color.Rda') # col.21 for 21 lymphoid subsets

clusters_lymphoid <- c('B_pro','B_pro_cycling','B_pre',
                       'B_mature','B_mature_B2m','Plasma',
                       'DNT_Il2ra','DPT_Rag1','DPT_cycling','DPT_Itm2a',
                       'T_naive_1','T_naive_2','T_naive_transitory',
                       'CD8T_cycling','CD8T_cyto','CD8T_resi_cyto',
                       'Treg','T_gd_1','T_gd_2','ILC','NK')
seu$annots <- factor(seu$annots, levels = clusters_lymphoid)

## colors
Time.color <- setNames(c('#DBBED7','#9282A7','#72589C'), c('P16','W4','W9'))
Group.color <- setNames(c('#EDBF36','#7CCCE7'), c('SPF','GF'))
Tissue.color <- setNames(c('#F3CD3E','#6F93CC','#8A9095','#CEAC80','#ED7248','#CFDD9F','#2277B3','#E25B3D','#437034','#71C6E8'),
                         c('bm','colon','heart','kidney','liver','lung','si','skin','spleen','thymus'))

custom_theme <- theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
          axis.ticks.x = element_blank())

#################################################################
## Fig. 4a
## UMAP of lymphoid cells colored by 21 subsets
DimPlot(seu, reduction = 'metacell', group.by = 'annots', cols = col.21,
        label = FALSE, raster = FALSE) & xlab('') & ylab('') & ggtitle('')
ggsave('lym/figures/Fig4a_lymphoid_umap.pdf', width = 6, height = 5)

#################################################################
## Fig. 4b
## Feature plots of representative marker genes across lymphoid subtypes
## color shows log2(UMI+1)
gl <- c('Cd3d','Cd8a','Cd4','Rag1','Itm2a','Il2ra',
        'Tuba1b','Foxp3','Ccr7','Gzmk','Ccl5','Trdc',
        'Xcl1','Cd93','Igll1','Gata3','Ms4a1','Vpreb1')
FeaturePlot(seu, features = gl, reduction = 'metacell', order = TRUE,
            cols = c('lightgrey', 'red'), ncol = 6) & NoAxes() & NoLegend()
ggsave('lym/figures/Fig4b_lymphoid_featureplot.pdf', width = 15, height = 9)

#################################################################
## Figs. 4c, 4e
## Stacked bar plots showing composition of 21 lymphoid subsets
## across 3 time points and 10 tissues (SPF for Fig4c, GF for Fig4e)

## load single-cell metacell info
mc_info <- read.csv('lym/output/seurat/lym_88180c_mc_info.csv', sep = ' ')

## per-mouse frequency table of 21 lymphoid subsets
mc_info %>%
    dplyr::select(Group, Time, Mouse, Tissue_only, annots) %>%
    mutate(comb = sprintf('%s_%s_%s_%s', Group, Time, Mouse, Tissue_only)) %>%
    dplyr::select(-Group, -Time, -Mouse, -Tissue_only) -> tmp

table(tmp$annots, tmp$comb) %>%
    as.data.frame() %>%
    dplyr::rename(annots = Var1, comb = Var2, n = Freq) -> ggtable_cnt

ggtable_cnt %>%
    dplyr::group_by(comb) %>%
    transmute(prop = n/sum(n), annots = annots) %>%
    ungroup() %>%
    tidyr::separate(comb, into = c('Group','Time','Mouse','Tissue_only'), sep = '_') -> ggtable_freq

tissue_levels <- c('si','heart','lung','liver','kidney','spleen','bm','skin','colon','thymus')

## Fig. 4c: stacked bar, SPF
ggtable_freq %>%
    filter(Group == 'SPF') %>%
    group_by(annots, Tissue_only, Time) %>%
    summarise(freq = mean(prop), .groups = 'drop') %>%
ggplot(aes(Time, freq, fill = factor(annots, levels = clusters_lymphoid))) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values = col.21) + theme_void() +
    theme(axis.text = element_text(size = 10), legend.title = element_blank(),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 15),
          legend.position = 'bottom') +
    guides(fill = guide_legend(title = "", ncol = 8)) +
    facet_grid(cols = vars(factor(Tissue_only, levels = tissue_levels)))
ggsave('lym/figures/Fig4c_spf_barplot_lym.pdf', width = 12, height = 4)

## Fig. 4e: stacked bar, GF
ggtable_freq %>%
    filter(Group == 'GF') %>%
    group_by(annots, Tissue_only, Time) %>%
    summarise(freq = mean(prop), .groups = 'drop') %>%
ggplot(aes(Time, freq, fill = factor(annots, levels = clusters_lymphoid))) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values = col.21) + theme_void() +
    theme(axis.text = element_text(size = 10), legend.title = element_blank(),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 15),
          legend.position = 'bottom') +
    guides(fill = guide_legend(title = "", ncol = 8)) +
    facet_grid(cols = vars(factor(Tissue_only, levels = tissue_levels)))
ggsave('lym/figures/Fig4e_gf_barplot_lym.pdf', width = 12, height = 4)

#################################################################
## Fig. 4d
## Temporal variation (W4/P16, W9/W4) of lymphoid subset composition in SPF mice
## P16 meanprop
ggtable_freq %>%
    filter(Time == 'P16') %>%
    select(-Time) %>%
    group_by(Group, Tissue_only, annots) %>%
    summarise(mean_prop = mean(prop)) -> P16meanProp

## W4 meanprop
ggtable_freq %>%
    filter(Time == 'W4') %>%
    select(-Time) %>%
    group_by(Group, Tissue_only, annots) %>%
    summarise(mean_prop = mean(prop)) -> W4meanProp

## W4 query on P16
ggtable_freq %>%
    filter(Time == 'W4') %>%
    select(-Time) -> W4_query

## W9 query on W4
ggtable_freq %>%
    filter(Time == 'W9') %>%
    select(-Time) -> W9_query

## W4/P16
merge(W4_query, P16meanProp, by = c('Group','Tissue_only','annots')) %>%
    mutate(abs_solo = abs(prop - mean_prop)) %>%
    select(-prop, -mean_prop) %>%
    group_by(Group, Tissue_only, Mouse) %>%
    summarise(abs_total = sum(abs_solo)) %>%
    ungroup() -> W4_P16_abs_df

## W9/W4
merge(W9_query, W4meanProp, by = c('Group','Tissue_only','annots')) %>%
    mutate(abs_solo = abs(prop - mean_prop)) %>%
    select(-prop, -mean_prop) %>%
    group_by(Group, Tissue_only, Mouse) %>%
    summarise(abs_total = sum(abs_solo)) %>%
    ungroup() -> W9_W4_abs_df

## point plot shows W9/W4, W4/P16 variation in SPF
W9_W4_abs_df$stageA2B = 'W9/W4'
W4_P16_abs_df$stageA2B = 'W4/P16'

df <- rbind(W9_W4_abs_df, W4_P16_abs_df) %>%
    mutate(stageA2B = factor(stageA2B, levels = c('W4/P16','W9/W4'))) %>%
    filter(Group == 'SPF') %>% as.data.frame()

df %>%
ggplot(aes(stageA2B, abs_total, fill = stageA2B)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5) +
    geom_jitter(shape = 21, size = 2, width = 0.2, stroke = 0.2) +
    facet_wrap(~factor(Tissue_only, levels = c('si','heart','lung','liver','kidney','spleen','bm','skin','colon','thymus')), ncol = 10) +
    scale_fill_manual(values = c('#3d53a4', '#c5499b')) +
    scale_color_manual(values = c('#3d53a4', '#c5499b')) +
    custom_theme +
    theme(legend.position = 'none', axis.text.x = element_blank()) + xlab('') + ylab('') +
    ggpubr::geom_pwc(hide.ns = FALSE, label = "p.format", p.adjust.method = 'BH',
                     label.size = 4, vjust = 0.7, method = 't_test')
ggsave('lym/figures/v9_spf_lym_cell_freq_abs_stageA2B.pdf', width = 6, height = 2)

#################################################################
## Fig. 4f
## Environmental variation (GF vs SPF) of lymphoid composition across time and tissues
## calculate SPF-Mean
ggtable_freq %>%
    filter(Group == 'SPF') %>%
    select(-Group) %>%
    group_by(Time, Tissue_only, annots) %>%
    summarise(mean_prop = mean(prop)) -> SPFmeanProp

## extract GF-duplicates
ggtable_freq %>%
    filter(Group == 'GF') %>%
    select(-Group) -> GF_query

## merge
merge(GF_query, SPFmeanProp, by = c('Time','Tissue_only','annots')) %>%
    mutate(abs_solo = abs(prop - mean_prop)) %>%
    select(-prop, -mean_prop) %>%
    group_by(Time, Tissue_only, Mouse) %>%
    summarise(abs_total = sum(abs_solo)) %>%
    ungroup() %>%
ggplot(aes(Time, abs_total, fill = Time)) +
    facet_wrap(~factor(Tissue_only, levels = c('si','lung','liver','heart','spleen','bm','kidney','colon','skin','thymus')), ncol = 10) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5) +
    geom_jitter(shape = 21, size = 2, width = 0.2, stroke = 0.2) +
    scale_fill_manual(values = c('#56B4E9', '#E69F00', '#D55E00')) +
    custom_theme +
    theme(axis.text.x = element_blank(), legend.position = 'none') + xlab('') + ylab('') +
    ggpubr::geom_pwc(hide.ns = TRUE, label = "p.format", p.adjust.method = 'BH',
                     label.size = 4, vjust = 0.7, method = 't_test')
ggsave('lym/figures/v9_lym_time_freq_abs_SPF2GF.pdf', width = 6, height = 2)

#################################################################
## Fig. 4g
## 2D projections of the environmental deviation onto the lymphoid metacell landscape
## small dots: single cells; large dots: metacells (size = cell number, color = GF%)

## to compare SPF and GF fairly, cells in each pool (Time x Tissue) were oversampled
## to the same size between SPF and GF
sc_info <- read.table('lym/output/seurat/lym_88180c_mc_info.csv')
sc_info <- sc_info %>% select(-Tissue, -Phase, -colors)

## input: meta data frame and time point; output: cell IDs needed to fill by resampling
getsampleindex <- function(mtdf = sc_info, time = 'P16'){
    # cells of one time point
    sc_subset <- mtdf %>% filter(Time == time)

    # cell number per Group per tissue
    cell_counts <- mtdf %>%
        dplyr::filter(Time == time) %>%
        dplyr::group_by(Group, Tissue_only) %>%
        dplyr::summarise(tis_n_cells = n(), .groups = 'drop')

    # the larger group per tissue
    max_counts <- cell_counts %>%
        dplyr::group_by(Tissue_only) %>%
        dplyr::summarise(max_cells = max(tis_n_cells), .groups = 'drop')

    # sample size to oversample
    sample_size <- cell_counts %>%
        left_join(max_counts, by = "Tissue_only") %>%
        mutate(sample_size = max_cells - tis_n_cells)
    sample_size <- sample_size %>% filter(sample_size > 0)

    add_container <- c()
    # oversample with replacement
    for (row in 1:NROW(sample_size)){
        group = as.character(sample_size[row, 'Group'])
        tissue = as.character(sample_size[row, 'Tissue_only'])
        size = as.integer(sample_size[row, 'sample_size'])
        temp <- mtdf[mtdf$Group == group & mtdf$Tissue_only == tissue,]
        sample_index <- sample(rownames(temp), size, replace = TRUE)
        add_container <- c(add_container, sample_index)
    }
    return(add_container)
}

# cells to add per time point
p16_add <- sc_info[getsampleindex(sc_info, 'P16'),]; rownames(p16_add) <- NULL
w4_add  <- sc_info[getsampleindex(sc_info, 'W4'),];  rownames(w4_add) <- NULL
w9_add  <- sc_info[getsampleindex(sc_info, 'W9'),];  rownames(w9_add) <- NULL

# all cells after oversampling
sc_info_sample <- sc_info %>% tibble::rownames_to_column(var = 'rowname') %>%
    select(-rowname) %>%
    rbind(p16_add, w4_add, w9_add)

# plot oversampled cells for per time
plot_time <- function(mtdf = sc_info_sample, time = 'P16') {
    # calculate Group% in per mc on specific time point
    mc_group_ratio <- mtdf %>%
        dplyr::filter(Time == time) %>%
        dplyr::group_by(mc_id, Group) %>%
        dplyr::count() %>%
        dplyr::group_by(mc_id) %>%
        dplyr::mutate(prop = 100*n/sum(n)) %>%
        dplyr::ungroup() %>%
        tidyr::complete(mc_id, Group, fill = list(n = 0, prop = 0))

    # mean mc_x, mc_y
    mc_mean_xy <- mtdf %>%
        dplyr::filter(Time == time) %>%
        dplyr::select(mc_id, x, y) %>%
        dplyr::group_by(mc_id) %>%
        dplyr::mutate(x = mean(x), y = mean(y)) %>%
        dplyr::distinct()

    # merge the above two dfs, save only GF
    mc_info <- merge(mc_group_ratio, mc_mean_xy, by = 'mc_id') %>%
        dplyr::group_by(mc_id) %>%
        dplyr::mutate(n_all = sum(n)) %>%
        dplyr::distinct(mc_id, .keep_all = TRUE) %>%
        # cal GF%mc, SPF as blank
        mutate(prop = if_else(Group == 'GF', prop, 100 - prop))

    # single cell data to plot
    plot_df <- mtdf %>% dplyr::filter(Time == time)
    plot_df <- merge(plot_df, mc_info[, c('mc_id','prop')], by = 'mc_id')

    # plot
    ggplot() +
        # sc point, color: prop
        geom_point(data = plot_df, aes(x, y, color = prop), shape = 19, size = 0.01, alpha = 0.5) +
        # mc point, color: prop + size: cell number
        geom_point(data = mc_info, aes(x, y, fill = prop, size = n_all), shape = 21) +
        # fix the point size to make it comparable between 3 time groups
        scale_size_continuous(range = c(1.5, 5)) +
        scale_fill_gradientn(colours = c("#EDBF36", "white", "#7CCCE7"), limits = c(0, 100), oob = scales::squish) +
        scale_color_gradientn(colours = c("#EDBF36", "white", "#7CCCE7"), limits = c(0, 100), oob = scales::squish) +
        theme_void() +
        labs(fill = "Fraction")
}

p1 <- plot_time(mtdf = sc_info_sample, time = 'P16')
p2 <- plot_time(mtdf = sc_info_sample, time = 'W4')
p3 <- plot_time(mtdf = sc_info_sample, time = 'W9')

ggsave(filename = 'lym/output/mc1/P16_SPFvsGF_mc2d.pdf', p1, width = 5, height = 4)
ggsave(filename = 'lym/output/mc1/W4_SPFvsGF_mc2d.pdf',  p2, width = 5, height = 4)
ggsave(filename = 'lym/output/mc1/W9_SPFvsGF_mc2d.pdf',  p3, width = 5, height = 4)

# save oversampled info for other heatmaps (e.g. Fig s4a-c annotation bars)
write.csv(sc_info_sample, 'lym/output/seurat/lym_105600c_mc_info_oversample.csv')

#################################################################
## Fig. 4h
## Frequencies of selected lymphoid subsets (B_mature_B2m, T_naive1, CD8T_cyto)
## in SPF vs GF mice across indicated tissues and developmental stages
## shaded line plots show mean with SEM, each dot represents one mouse

## prepare per-mouse frequency table (ggprop_d)
sc_info %>%
    mutate(Tissue_only = recode(Tissue_only,
                                'bm' = 'BM','colon' = 'CL','heart' = 'HT','kidney' = 'KD','liver' = 'LV',
                                'lung' = 'LU','si' = 'SI','skin' = 'SK','spleen' = 'SP','thymus' = 'TM')) -> sc_info

sc_info %>%
    dplyr::select(Group, annots, Time, Tissue_only, Mouse) %>%
    dplyr::mutate(comb = sprintf('%s_%s_%s_%s', Group, Time, Tissue_only, Mouse)) %>%
    dplyr::select(-Time, -Tissue_only, -Group, -Mouse) -> tmp

ggumap_ct <- table(tmp$annots, tmp$comb) %>% as.data.frame() %>%
    dplyr::rename(annots = Var1, comb = Var2, n = Freq)

ggumap_ct %>%
    dplyr::group_by(comb) %>%
    dplyr::filter(sum(n) >= 20) %>% # filter mouse with less than 20 cells
    transmute(Fraction = 100*n/sum(n), annots = annots) -> ggprop_d

ggprop_d$Time        = sapply(ggprop_d$comb, function(x) str_split(x, '_')[[1]][2])
ggprop_d$Tissue_only = sapply(ggprop_d$comb, function(x) str_split(x, '_')[[1]][3])
ggprop_d$Group       = sapply(ggprop_d$comb, function(x) str_split(x, '_')[[1]][1])
ggprop_d$Mouse       = sapply(ggprop_d$comb, function(x) str_split(x, '_')[[1]][4])

## line layers for shaded line plots (mean +- SEM)
pd <- position_dodge(width = 0)
line_layers <- function() {
    list(
        stat_summary(fun.data = mean_se, geom = "ribbon",
                     aes(group = Group, fill = Group), alpha = 0.15, color = NA),
        geom_boxplot(width = 0.8, alpha = 0.2, position = pd, outlier.shape = NA),
        geom_jitter(size = 2, shape = 21, fill = NA, position = pd),
        # mean line
        stat_summary(fun = mean, geom = "line", aes(group = Group)),
        scale_fill_manual(values = c('#77BBD1', '#DAB33D')),
        scale_color_manual(values = c('#77BBD1', '#DAB33D')),
        # statistics: SPF vs GF within each time point
        stat_compare_means(aes(group = Group), method = "t.test", label = "p",
                           hide.ns = TRUE, label.y.npc = 0.95,
                           tip.length = 0.02, bracket.size = 0.5),
        facet_wrap(~ Tissue_only, scales = 'free_y', ncol = 8),
        theme_classic(),
        theme(axis.text.x = element_text(angle = 0, size = 10),
              strip.background = element_rect(fill = "gray95"),
              strip.text = element_text(face = "bold"),
              legend.position = "none"),
        ylab('Freq (%)'),
        xlab('Time')
    )
}

`%out%` <- Negate(`%in%`)

## Fig. 4h: T_naive1 (skin, SI, colon and BM excluded)
ggprop_d %>%
    dplyr::filter(Tissue_only %out% c('SK', 'SI', 'CL', 'BM')) %>%
    dplyr::filter(annots %in% c('T_naive_1')) %>%
    mutate(Time = factor(Time), Group = factor(Group)) %>%
    ggplot(aes(x = Time, y = Fraction, fill = Group, color = Group)) +
    line_layers()
ggsave('lym/figures/T_naive1_line_bar.pdf', width = 8, height = 2.5)

## Fig. 4h: B_mature_B2m (thymus, skin, colon and SI excluded)
ggprop_d %>%
    dplyr::filter(Tissue_only %out% c('TM', 'SK', 'CL', 'SI')) %>%
    dplyr::filter(annots %in% c('B_mature_B2m')) %>%
    mutate(Time = factor(Time),
           Group = factor(Group),
           Tissue_only = factor(Tissue_only, levels = c('HT','KD','LU','BM','LV','SP'))) %>%
    ggplot(aes(x = Time, y = Fraction, fill = Group, color = Group)) +
    line_layers()
ggsave('lym/figures/B_mature_B2m_bar.pdf', width = 8, height = 2.5)

## Fig. 4h: CD8T_cyto (thymus and skin excluded)
ggprop_d %>%
    dplyr::filter(Tissue_only %out% c('TM', 'SK')) %>%
    dplyr::filter(annots %in% c('CD8T_cyto')) %>%
    mutate(Time = factor(Time),
           Group = factor(Group),
           Tissue_only = factor(Tissue_only, levels = c('HT','KD','LU','LV','SI','CL','BM','SP'))) %>%
    ggplot(aes(x = Time, y = Fraction, fill = Group, color = Group)) +
    line_layers()
ggsave('lym/figures/CD8T_cyto_line.pdf', width = 10.8, height = 2.5)

#################################################################
## Figs. s4a-s4c
## MetaCell-based annotation heatmaps of B, T and ILC subsets
## marker genes displayed as scaled expression (footprint) per metacell
require(metacell)
figpath = 'lym/output/mc1'
dbpath  = 'lym/data/mc1'
scdb_init(dbpath, force_reinit = TRUE)
scfigs_init(figpath)
id = 'mouse_v9_lym'
mat  = scdb_mat(id)
mc   = scdb_mc(paste0(id, '_filt'))
mc2d = scdb_mc2d(paste0(id, '_filt'))
gset = scdb_gset(paste0(id, '_final'))

## metacell order per lineage (262 metacells in total)
T.mc.id <- c(226:252,      # DPT_Rag1
             258,          # DNT_Il2ra
             259:262,      # DPT_cycling
             187:189, 253:257, # DPT_Itm2a
             203, 206:209, # T_naive_1
             197:202, 204:205, 210:225, # T_naive_2
             190:196,      # T_naive_transitory
             140, 165, 166, # Treg
             113:119,      # T_gd_1
             120:129,      # T_gd_2
             184:186,      # CD8T_cycling
             141:149,      # CD8T_resi_cyto
             167:182)      # CD8T_cyto
B.mc.id <- c(105, 106,     # B_pro
             107:112,      # B_pro_cycling
             94:103,       # B_pre
             1:15, 17:46, 58:89, 92:93, # B_mature
             16, 49:57, 90:91) # B_mature_B2m
ILC.mc.id <- c(130:139,    # ILC
               150:164)    # NK

## marker genes per lineage
T.markers <- c('Cd3d','Cd3e','Cd3g', # common
               'Thy1','Rgcc','Trat1','Nsg2','Cd4','Bcl11b','Themis','Mir8113','Ccr9','Satb1','Trbc1','Ets2','Rmnd5a','Ldhb','Arpp21','Aqp11','Myb','Trp53inp1','Dntt','Rag1','Trac','Trbc2','Tcf7','Cd8a','Cd8b1', # DNT_Il2ra
               'Il2ra','Uck2','Fkbp5', # DNT_Il2ra
               'Mki67','Top2a','Ube2c','Rrm2','Hmgb2','Pcna','Dut','Tuba1b','Tubb5','Stmn1', # DPT_cycling
               'Itm2a','Lef1','Cd28','Cd5','Rgs10', # DPT_Itm2a
               'Dapl1','Igfbp4','Slfn1','Il7r','Ms4a6b','Ms4a4b','Hcst','Fxyd5','Dgka','Selenop','Emb','Mb21d1','Stat1','Ifi27l2a','Gimap4','Tmsb10','S1pr1','Ccr7','Sell','Klf2','Ifit1','Lars2', # T_naive
               'Maf','Rora','Tnfrsf4','Ctla4','Icos','Ikzf2', # Treg
               'Lmna','Gem','Klrb1b','Ctla2b','Arl4c','Crem','Nr4a3','Nr4a2','Klf4','Fos','Rgs2','Vim','Itgb1','Cxcr6','Lgals1','Ahnak','Lgals3','S100a11','Crip1','S100a10','S100a6','S100a4','Ctla2a','Ccnd2','Dusp1','Tnfaip3','Id2','Il2rb','Rgs1','Nkg7','Ctsw','Cd7','Klrd1','Tcrg-C1','Trdc','Fcer1g','Xcl1','Tcrg-C4','Tcrg-C2', # T_gd
               'Cd160','Klra7','Tyrobp','Igha', # CD8T_resi_cyto
               'H2-Q7','Ly6c2','Ccl5','Gzmb','Gzma' # CD8T_cyto
               )
B.markers <- c("Ighm","Cd79b","Cd79a","Chchd10","Spib","Vpreb3","Pou2af1","Fcrla","Ebf1","Blnk","Mzb1","Syk","Cd24a", # general
               "Pafah1b3","Tifa","Cd93","Myb","Smarca4","Rag1", # B progenitor general
               "Vpreb1","Igll1","Dntt","Vpreb2", # B_pro
               "Mki67","Top2a","Stmn1","Rrm2","Tubb5","Hmgb2","Ccna2","Ube2c","H2afx", # B_pro_cycling
               "Fcmr","Cd74","Ms4a1","Iglj3","Iglj2", # B_pre
               "Ly6d","Scd1","Iglj1","Tnfrsf13c","Ifi30",
               "H2-Eb1","H2-Aa","H2-Ab1","Ighd","H2-DMb2","Mir5107","Bank1","Cd83", # B_mature
               "B2m") # B_mature_B2m
ILC.markers <- c("Ahnak","S100a6","Id2","Tnfaip3","Dusp1",
                 "Fgl2","Rgs1","Cxcr6","S100a4","Selenop",
                 "Il22","Gzmc","Tcrg-C1","Trdc",
                 "Itgae","Lmna","Ccr8","Nr4a3",
                 "Crem","Odc1","Tmem176a","Tmem176b",
                 "Arg1","Gata3","Hilpda","Areg","Rora","Icos",
                 "Il7r","Emb","Junb","Rgs2","Fosb","Gadd45b",
                 "Fos","Nfkbia","Nr4a1","Nfkbiz",
                 "Vim","Lgals1",
                 "S100a10","S100a11","Fxyd5","Ccr2","Ctla2a", # ILC
                 "Bhlhe40","Klrb1b","Klra4","Klre1","Gzma","Ncr1","Klrd1","Tyrobp","Klrb1c","Ccl5","Ccl3","Ccl4","Gzmb","Ctsw","Fcer1g","Ms4a4b","Xcl1","Cd7","Klrk1","Nkg7","Il2rb","Prf1" # NK
                 )

## Fig. s4a: B subtypes
tgconfig::set_param("mcp_heatmap_height", 40, "metacell")
tgconfig::set_param("mcp_heatmap_width", 20, "metacell")
tgconfig::set_param("mcp_heatmap_text_cex", 1, "metacell")
mcell_mc_plot_marks2(paste0(id, '_filt'), gset_id = id, mat_id = id,
    plot_cells = FALSE, mc_ord = B.mc.id, focus_mcs = B.mc.id,
    gene_list = B.markers, reorder_marks = FALSE,
    fig_fn = sprintf('%s/mouse_v9.mc_heat_marks_B_selected_markers.pdf', figpath),
    zero_median = FALSE)

## Fig. s4b: T subtypes
tgconfig::set_param("mcp_heatmap_height", 60, "metacell")
mcell_mc_plot_marks2(paste0(id, '_filt'), gset_id = id, mat_id = id,
    plot_cells = FALSE, mc_ord = T.mc.id, focus_mcs = T.mc.id,
    gene_list = T.markers, reorder_marks = FALSE,
    fig_fn = sprintf('%s/mouse_v9.mc_heat_marks_T_selected_markers.pdf', figpath),
    zero_median = FALSE)

## Fig. s4c: ILC subtypes
tgconfig::set_param("mcp_heatmap_height", 40, "metacell")
mcell_mc_plot_marks2(paste0(id, '_filt'), gset_id = id, mat_id = id,
    plot_cells = FALSE, mc_ord = ILC.mc.id, focus_mcs = ILC.mc.id,
    gene_list = ILC.markers, reorder_marks = FALSE,
    fig_fn = sprintf('%s/mouse_v9.mc_heat_marks_ILC_selected_markers.pdf', figpath),
    zero_median = FALSE)

## annotation bars for Figs. s4a-s4c (Group/Tissue/Time composition + cell type per metacell)
require(aplot)
mc_info_os <- read.csv('lym/output/seurat/lym_105600c_mc_info_oversample.csv', row.names = 1)

plot_anno_bars <- function(mc_ids, prefix){
    ## group annotation
    group_bar <- mc_info_os %>%
        filter(mc_id %in% mc_ids) %>%
        select(mc_id, Group) %>%
        group_by(mc_id) %>%
        count(Group) %>%
        mutate(mc_id = as.character(mc_id), prop = n/sum(n)) %>%
    ggplot() +
        geom_bar(aes(factor(mc_id, level = mc_ids), prop, fill = Group), size = 0.2, stat = 'identity', width = 1) +
        scale_fill_manual(values = Group.color) +
        theme_void()

    ## tissue annotation
    tis_bar <- mc_info_os %>%
        filter(mc_id %in% mc_ids) %>%
        select(mc_id, Tissue_only) %>%
        mutate(Tissue_only = recode(Tissue_only,
                                    'bm' = 'BM','colon' = 'CL','heart' = 'HT','kidney' = 'KD','liver' = 'LV',
                                    'lung' = 'LU','si' = 'SI','skin' = 'SK','spleen' = 'SP','thymus' = 'TM')) %>%
        group_by(mc_id) %>%
        count(Tissue_only) %>%
        mutate(mc_id = as.character(mc_id), prop = n/sum(n)) %>%
    ggplot() +
        geom_bar(aes(factor(mc_id, level = mc_ids), prop, fill = Tissue_only), size = 0.2, stat = 'identity', width = 1) +
        scale_fill_manual(values = Tissue.color) + theme_void()

    ## time annotation
    time_bar <- mc_info_os %>%
        filter(mc_id %in% mc_ids) %>%
        select(mc_id, Time) %>%
        group_by(mc_id) %>%
        count(Time) %>%
        mutate(mc_id = as.character(mc_id), prop = n/sum(n)) %>%
    ggplot() +
        geom_bar(aes(factor(mc_id, level = mc_ids), prop, fill = Time), size = 0.2, stat = 'identity', width = 1) +
        scale_fill_manual(values = Time.color) + theme_void()

    ## merge
    group_bar %>% insert_top(time_bar) %>% insert_top(tis_bar)
    ggsave(sprintf('lym/output/mc1/heatmap/%s_specific_pct_bar.pdf', prefix), width = 8, height = 4)

    ## celltype annotation
    celltype_bar <- data.frame(id = 1:length(mc@annots[mc_ids]), celltype = mc@annots[mc_ids]) %>%
        mutate(celltype = factor(celltype)) %>%
    ggplot(aes(id, 1, fill = celltype)) +
        geom_tile() +
        scale_fill_manual(values = col.21) +
        theme_void()
    ggsave(sprintf('lym/output/mc1/heatmap/%s_celltype_bar.pdf', prefix), width = 8, height = 1)
}

plot_anno_bars(B.mc.id, 'B')
plot_anno_bars(T.mc.id, 'T')
plot_anno_bars(ILC.mc.id, 'ILC')

#################################################################
## Figs. s4d-s4f
## Scatter plots comparing gene expression between paired lymphoid subsets
## (log2(expression+1) of average normalized expression per subset)
## NOTE: the original plotting cells for these panels were not part of the
## provided notebooks; the code below follows the published panels
plot_pair_scatter <- function(seu, ct_x, ct_y, out_pdf, lfc_cut = 1, max.overlaps = 20){
    require(ggrepel)
    ## average expression per subset (AverageExpression returns log1p(mean(expm1(x))))
    mat <- AverageExpression(seu, group.by = 'annots', slot = 'data')$RNA
    df <- data.frame(gene = rownames(mat),
                     x = log2(expm1(mat[, ct_x]) + 1),
                     y = log2(expm1(mat[, ct_y]) + 1))
    df$sign <- ifelse(df$x - df$y > lfc_cut, ct_x,
               ifelse(df$y - df$x > lfc_cut, ct_y, 'NS'))
    ggplot(df, aes(x, y)) +
        geom_point(data = df %>% filter(sign == 'NS'), color = 'grey', size = 0.5) +
        geom_point(data = df %>% filter(sign == ct_x), color = 'red', size = 0.8) +
        geom_point(data = df %>% filter(sign == ct_y), color = 'blue', size = 0.8) +
        geom_text_repel(data = df %>% filter(sign != 'NS'), aes(label = gene),
                        size = 3, max.overlaps = max.overlaps) +
        geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'grey') +
        xlab(paste0(ct_x, ' log2(expression+1)')) +
        ylab(paste0(ct_y, ' log2(expression+1)')) +
        theme_classic()
    ggsave(out_pdf, width = 5, height = 5)
}

## Fig. s4d: B_mature_B2m vs B_mature
plot_pair_scatter(seu, 'B_mature_B2m', 'B_mature',
                  'lym/figures/B_mature_B2m_vs_B_mature_scatter.pdf')
## Fig. s4e: CD8T_cyto vs CD8T_resi_cyto
plot_pair_scatter(seu, 'CD8T_cyto', 'CD8T_resi_cyto',
                  'lym/figures/CD8T_cyto_vs_CD8T_resi_cyto_scatter.pdf')
## Fig. s4f: T_naive_1 vs T_naive_2
plot_pair_scatter(seu, 'T_naive_1', 'T_naive_2',
                  'lym/figures/T_naive_1_vs_T_naive_2_scatter.pdf')
