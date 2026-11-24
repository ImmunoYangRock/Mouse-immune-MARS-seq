#################################################################
## R code to create figures related to Figure2
## Figs 2a-2h and Figs s2a-s2i
## Description: compositional analysis of myeloid subsets —
## stacked-bar frequencies per tissue, temporal variation in SPF,
## GF-vs-SPF deviation, oversampled GF% projection on the metacell
## map, and per-subset frequency line/bar plots.
#################################################################

## load data
## dataframe ggumap from Fig1
require(dplyr)
require(tidyr)
require(tibble)
require(stringr)
require(ggplot2)
require(ggpubr)

# add lineage info
ggumap <- ggumap |> mutate(lineage = case_when(annot_r %in% c('B','ILC','T') ~ 'lym',
                                       annot_r %in% c('Granulocyte','DC','Mon-Mac') ~'mye')) 

ggumap |> 
    dplyr::select(Group,Time,Mouse,Tissue_only,lineage) |>
    dplyr::mutate(Mouse = str_replace_all(Mouse, "_", "-")) |>
    mutate(comb = sprintf('%s_%s_%s_%s',Group, Time, Tissue_only, Mouse)) |>    
    dplyr::select(-Time,-Tissue_only,-Group,-Mouse) -> tmp

## generate cross count table
table(tmp$lineage,tmp$comb) |> 
    as.data.frame() |> 
    dplyr::rename(lineage = Var1, comb = Var2, n = Freq) -> ggtable_cnt

## calculate proportions
tis_change = c('BM'='bm','TM'='thymus', 'SP'='spleen', 'LV'='liver', 'LU'='lung', 'SI'='si', 'CL'='colon', 'SK'='skin', 'KD'='kidney','HT'='heart')
ggtable_freq <- ggtable_cnt |> 
    dplyr::group_by(comb) |>
    mutate(n2 = sum(n)) %>%
    transmute(Freq = 100*n/sum(n), lineage = lineage) |>
    ungroup() |>
    # separated the crossed combinations 
    tidyr::separate(comb, into=c('Group','Time','Tissue_only','Mouse'),sep='_') |> 
    mutate(Tissue_only= plyr::mapvalues(Tissue_only, from = tis_change, to = names(tis_change)))               

#################################################################
## Fig. s2a 

ggtable_freq |>
 mutate(Tissue_only = factor(Tissue_only, levels = c('SP','TM','BM','LV','KD','SK','LU','HT','CL','SI')))|>
 ggplot(aes(Time,Freq,fill=lineage))+
 geom_jitter(shape=21,size=2,width=0.1)+
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = lineage, color = lineage),
    linewidth = 1
  ) +
 scale_fill_manual(values=c('#91C25E','#E8461E'))+
 scale_color_manual(values=c('#91C25E','#E8461E'))+
 facet_wrap(~Tissue_only+Group,nrow = 1)+
 #facet_grid(rows = vars(Group), cols = vars(Tissue_only))+
 theme_classic()+ylab('Freq(%)')
ggsave('figures/FigS2a.pdf', width=24, height=4)


#################################################################
## Fig. 2a 

## lym: across-annot_r
tmp <- ggumap |> 
    filter(lineage=='lym')%>%
    dplyr::select(Group,Time,Mouse,Tissue_only,annot_r) |>
    dplyr::mutate(Mouse = str_replace_all(Mouse, "_", "-")) |>
    mutate(comb = sprintf('%s_%s_%s_%s',Group, Time, Tissue_only, Mouse)) |>    
    dplyr::select(-Time,-Tissue_only,-Group,-Mouse) 

## generate cross count table
ggtable_cnt <- table(tmp$annot_r,tmp$comb) |> 
    as.data.frame() |> 
    dplyr::rename(annot_r = Var1, comb = Var2, n = Freq) 

## calculate proportions
ggtable_freq <- ggtable_cnt |> 
    dplyr::group_by(comb) |>
    mutate(n2 = sum(n),Freq = 100*n/sum(n)) |>
    ## separated the crossed combinations 
    tidyr::separate(comb, into=c('Group','Time','Tissue_only','Mouse'),sep='_') 


## P16 meanFreq
P16meanFreq  <- ggtable_freq |>
    filter(n2>=10) |>
    filter(Time=='P16') |>
    select(-Time) |>
    ungroup() |>
    group_by(Group, Tissue_only, annot_r) |>
    summarise(meanFreq = mean(Freq)) 

## W4 query on P16
ggtable_freq |>
    filter(n2>=10) |>
    filter(Time=='W4') |>
    select(-Time) -> W4_query #(W4, Group, Tissue_only,annots,prop)

## W4/P16
lym_W4_P16_abs_df <- merge(W4_query,P16meanFreq, by=c('Group','Tissue_only','annot_r')) |>
    mutate(abs_solo = abs(Freq-meanFreq)) |>
    select(-Freq,-meanFreq) |>
    group_by(Group,Tissue_only,Mouse) |>
    summarise(abs_total = sum(abs_solo)) |>
    ungroup() 


## mye: across-annot_r
tmp <- ggumap |> 
    filter(lineage=='mye')%>%
    dplyr::select(Group,Time,Mouse,Tissue_only,annot_r) |>
    dplyr::mutate(Mouse = str_replace_all(Mouse, "_", "-")) |>
    mutate(comb = sprintf('%s_%s_%s_%s',Group, Time, Tissue_only, Mouse)) |>    
    dplyr::select(-Time,-Tissue_only,-Group,-Mouse) 

## generate cross count table
ggtable_cnt <- table(tmp$annot_r,tmp$comb) |> 
    as.data.frame() |> 
    dplyr::rename(annot_r = Var1, comb = Var2, n = Freq) 

## calculate proportions
ggtable_freq <- ggtable_cnt |> 
    dplyr::group_by(comb) |>
    mutate(n2 = sum(n),Freq = 100*n/sum(n)) |>
    ## separated the crossed combinations 
    tidyr::separate(comb, into=c('Group','Time','Tissue_only','Mouse'),sep='_') 


## P16 meanFreq
P16meanFreq  <- ggtable_freq |>
    filter(n2>=10) |>
    filter(Time=='P16') |>
    select(-Time) |>
    ungroup() |>
    group_by(Group, Tissue_only, annot_r) |>
    summarise(meanFreq = mean(Freq)) 

## W4 query on P16
ggtable_freq |>
    filter(n2>=10) |>
    filter(Time=='W4') |>
    select(-Time) -> W4_query #(W4, Group, Tissue_only,annots,prop)

## W4/P16
mye_W4_P16_abs_df <- merge(W4_query,P16meanFreq, by=c('Group','Tissue_only','annot_r')) |>
    mutate(abs_solo = abs(Freq-meanFreq)) |>
    select(-Freq,-meanFreq) |>
    group_by(Group,Tissue_only,Mouse) |>
    summarise(abs_total = sum(abs_solo)) |>
    ungroup() 


merge(mye_W4_P16_abs_df,subset(lym_W4_P16_abs_df, Tissue_only != 'thymus'),by=c('Group','Tissue_only','Mouse')) %>%
    rename(mye=abs_total.x,lym=abs_total.y) %>%
    pivot_longer(cols = c(mye, lym), names_to = "lineage", values_to = "variation" ) %>%
ggplot(aes(lineage,variation/100,fill=lineage)) +
     geom_bar(stat = "summary", fun = "mean", width = 0.5, position = position_dodge()) +
     stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
     geom_jitter(shape = 21, size = 2, width = 0.2, stroke = 0.2) +
     facet_wrap(~Group+factor(Tissue_only,levels=c('skin','spleen','si','liver','kidney','lung','heart','colon','bm')),ncol=9)+
     scale_fill_manual(values=c('#96CB5E', '#FF410C'))+
     scale_color_manual(values=c('#96CB5E', '#FF410C'))+
    geom_pwc(hide.ns = T,label = "p.format",p.adjust.method = 'BH',label.size=4, vjust=0.7, method='t_test')
    theme_classic()+
    theme(legend.position = 'None',   
          plot.title = element_text(hjust=0.5, vjust=0.5),
          axis.ticks.x = element_blank())+xlab('')+ylab('')         

#################################################################
## Fig. 2b
## UMAP of 34,010 myeloid cells colored by 23 subsets
require(Seurat)
seu <- readRDS('mye/data/seurat/v9_mye_mc2seu.rds') # myeloid seurat object, metacell 2d projection stored in reduction 'metacell'
load('mye/data/mye_annot_color.Rda') # col.24 for 23 myeloid subsets

clusters_myeloid <- c("Eosinophil","Basophil","MastCell",
                      "Neu_bm","Neu_bm_Ccl6(Neg)","Neu_bm_Stfa2","Neu_Il1b_Cxcl2",
                      "cDC1","cDC2_Mdh2","cDC2_Ms4a4c","cDC2_skin_Retnla(high)","cDC2_skin_Retnla(low)","mDC","pDC",
                      "KupfferCell","Langerhans","Mac_gut_Mmp9","Mac_Il1b","Mac_Kidney","Mac_skin","Mac_skin_Retnla",
                      "Mon","Mon_bm")
seu$annots <- factor(seu$annots, levels = clusters_myeloid)

DimPlot(seu, reduction = 'metacell', group.by = 'annots', cols = col.24,
        label = FALSE, raster = FALSE) & xlab('') & ylab('') & ggtitle('')
ggsave('figures/Fig2b_myeloid_umap.pdf', width = 6, height = 5)

#################################################################
## Fig. 2c
## Feature plots of representative marker genes across myeloid subtypes
## color shows log2(UMI+1)
gl <- c('Csf1r','Nr4a1', # Mon
    'Lcn2','Hmgb2', # Mon_bm
    'S100a4','Apoe','Ms4a4a', # Mon total
    'Ltf','Chil3', # bm Ccl6
    'Stfa1','Stfa2','Stfa3', # bm Stfa2
    'Il1b','Cxcl2', # Il1b Cxcl2
    'Cd5l','Vsig4', # Kupffer
    'Mmp9','Dnase1l3','Apol7c', # Gut
    'Hexb','Slamf9', # Kidney
    'Cbr2','Ccl9','Clec10a', # Skin
    'Folr2','Lyve1', # skin_Folr2
    'Ccl3','Ccl4' # Mac_Il1b
     )
FeaturePlot(seu, features = gl, reduction = 'metacell', order = TRUE,
            cols = c('lightgrey', 'red'), ncol = 5) & NoAxes() & NoLegend()
ggsave('figures/Fig2c_myeloid_featureplot.pdf', width = 15, height = 12)

#################################################################
## Figs. 2d, 2f
## Stacked bar plots showing composition of 23 myeloid subsets
## across 3 time points and 10 tissues (SPF for Fig2d, GF for Fig2f)

## load single-cell metacell info
mc_info <- read.csv('mye/output/seurat/mye_32469c_mc_info.csv', sep = ' ')

## per-mouse frequency table of 23 myeloid subsets
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

## Fig. 2d: stacked bar, SPF
ggtable_freq %>%
    filter(Group == 'SPF') %>%
    group_by(annots, Tissue_only, Time) %>%
    summarise(freq = mean(prop), .groups = 'drop') %>%
ggplot(aes(Time, freq, fill = factor(annots, levels = clusters_myeloid))) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values = col.24) + theme_void() +
    theme(axis.text = element_text(size = 10), legend.title = element_blank(),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 15),
          legend.position = 'bottom') +
    guides(fill = guide_legend(title = "", ncol = 8)) +
    facet_grid(cols = vars(factor(Tissue_only, levels = c('spleen','bm','liver','thymus','lung','si','colon','heart','kidney','skin'))))
ggsave('figures/Fig2d_spf_barplot_mye.pdf', width = 12, height = 4)

## Fig. 2f: stacked bar, GF
ggtable_freq %>%
    filter(Group == 'GF') %>%
    group_by(annots, Tissue_only, Time) %>%
    summarise(freq = mean(prop), .groups = 'drop') %>%
ggplot(aes(Time, freq, fill = factor(annots, levels = clusters_myeloid))) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values = col.24) + theme_void() +
    theme(axis.text = element_text(size = 10), legend.title = element_blank(),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 15),
          legend.position = 'bottom') +
    guides(fill = guide_legend(title = "", ncol = 8)) +
    facet_grid(cols = vars(factor(Tissue_only, levels = c('spleen','bm','liver','thymus','lung','si','colon','heart','kidney','skin'))))
ggsave('figures/Fig2f_gf_barplot_mye.pdf', width = 12, height = 4)

#################################################################
## Fig. 2e
## Temporal variation (W4/P16, W9/W4) of myeloid subset composition in SPF mice
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

custom_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 0.5),
    axis.ticks.x = element_blank()
  )

## point plot shows W9/W4, W4/P16 variation in SPF
W9_W4_abs_df$stageA2B = 'W9/W4'
W4_P16_abs_df$stageA2B = 'W4/P16'

df <- rbind(W9_W4_abs_df, W4_P16_abs_df) %>%
    mutate(stageA2B = factor(stageA2B, levels = c('W4/P16','W9/W4'))) %>%
    filter(Group == 'SPF') %>% as.data.frame()

df %>%
    filter(Tissue_only != 'thymus') %>%
ggplot(aes(stageA2B, abs_total, fill = stageA2B)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5) +
    geom_jitter(shape = 21, size = 2, width = 0.2, stroke = 0.2) +
    facet_wrap(~factor(Tissue_only, levels = c('heart','kidney','spleen','lung','colon','bm','liver','si','skin')), ncol = 10) +
    scale_fill_manual(values = c('#3d53a4', '#c5499b')) +
    scale_color_manual(values = c('#3d53a4', '#c5499b')) +
    custom_theme +
    theme(axis.text.x = element_blank(), legend.position = 'none') + xlab('') + ylab('') +
    geom_pwc(hide.ns = F, label = "p.format", p.adjust.method = 'BH', label.size = 4, vjust = 0.7, method = 't_test')
ggsave('figures/v9_spf_mye_cell_freq_abs_stageA2B.pdf', width = 6, height = 2)

#################################################################
## Fig. 2g
## Environmental variation (GF vs SPF) of myeloid composition across time and tissues
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
    filter(Tissue_only != 'thymus') %>%
ggplot(aes(Time, abs_total, fill = Time)) +
    facet_wrap(~factor(Tissue_only, levels = c('spleen','colon','heart','skin','liver','kidney','si','lung','bm')), ncol = 10) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5) +
    geom_jitter(shape = 21, size = 2, width = 0.2, stroke = 0.2) +
    scale_fill_manual(values = c('#56B4E9', '#E69F00', '#D55E00')) +
    custom_theme +
    theme(axis.text.x = element_blank(), legend.position = 'none') + xlab('') + ylab('') +
    geom_pwc(hide.ns = F, label = "p.format", p.adjust.method = 'BH', label.size = 4, vjust = 0.7, method = 't_test')
ggsave('figures/v9_mye_time_freq_abs_SPF2GF.pdf', width = 6, height = 2)

#################################################################
## Fig. 2h
## 2D projections of the environmental deviation onto the myeloid metacell landscape
## small dots: single cells; large dots: metacells (size = cell number, color = GF%)

## to compare SPF and GF fairly, cells in each pool (Time x Tissue) were oversampled
## to the same size between SPF and GF
sc_info <- read.table('mye/output/seurat/mye_32469c_mc_info.csv')
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
    group = as.character(sample_size[row,'Group'])
    tissue = as.character(sample_size[row,'Tissue_only'])
    size = as.integer(sample_size[row,'sample_size'])
    temp <- mtdf[mtdf$Group == group & mtdf$Tissue_only == tissue,]
    sample_index <- sample(rownames(temp), size, replace = TRUE)
    add_container <- c(add_container, sample_index)
  }
  return(add_container)
}

# cells to add per time point
p16_add <- sc_info[getsampleindex(sc_info,'P16'),]; rownames(p16_add) <- NULL
w4_add <- sc_info[getsampleindex(sc_info,'W4'),]; rownames(w4_add) <- NULL
w9_add <- sc_info[getsampleindex(sc_info,'W9'),]; rownames(w9_add) <- NULL

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
  plot_df <- merge(plot_df, mc_info[,c('mc_id','prop')], by = 'mc_id')

  # plot
  ggplot() +
    # sc point, color: prop
    geom_point(data = plot_df, aes(x, y, color = prop), shape = 19, size = 0.01, alpha = 0.5) +
    # mc point, color: prop + size: cell number
    geom_point(data = mc_info, aes(x, y, fill = prop, size = n_all), shape = 21) +
    # fix the point size to make it comparable between 3 time groups
    scale_size_continuous(range = c(1.5, 5)) +
    scale_fill_gradientn(colours = c("#EDBF36", "white", "#7CCCE7"), limits = c(0,100), oob = scales::squish) +
    scale_color_gradientn(colours = c("#EDBF36", "white", "#7CCCE7"), limits = c(0,100), oob = scales::squish) +
    theme_void() +
    labs(fill = "Fraction")
}

p1 <- plot_time(mtdf = sc_info_sample, time = 'P16')
p2 <- plot_time(mtdf = sc_info_sample, time = 'W4')
p3 <- plot_time(mtdf = sc_info_sample, time = 'W9')

ggsave(filename = 'mye/output/mc1/P16_SPFvsGF_mc2d.pdf', p1, width = 5, height = 4)
ggsave(filename = 'mye/output/mc1/W4_SPFvsGF_mc2d.pdf', p2, width = 5, height = 4)
ggsave(filename = 'mye/output/mc1/W9_SPFvsGF_mc2d.pdf', p3, width = 5, height = 4)

# save oversampled info for other heatmaps (e.g. Fig s2b-e annotation bars)
write.csv(sc_info_sample, 'mye/output/seurat/mye_366641c_mc_info_oversample.csv')

#################################################################
## Figs. s2b-s2e
## MetaCell-based annotation heatmaps of Mon, Mac, Neu and DC subsets
## marker genes displayed as scaled expression (footprint) per metacell
require(metacell)
figpath = 'mye/output/mc1'
dbpath  = 'mye/data/mc1'
scdb_init(dbpath, force_reinit = TRUE)
scfigs_init(figpath)
id = 'mouse_v9_mye'
mat  = scdb_mat(id)
mc   = scdb_mc(paste0(id, '_filt'))
mc2d = scdb_mc2d(paste0(id, '_filt'))
gset = scdb_gset(paste0(id, '_final'))

## ranked metacell order per lineage
Neu.mc.id <- c(72, 83, 84, 91,        # Stfa2
               69, 70, 71, 82, 73, 74, 75, # Cxcl2 Il1b
               76:81,                 # Neu_bm
               85:90, 92:95           # Ltf, Ccl6-
               )
Mon.mc.id <- c(14, 15, 16, 100:119,   # Mon
               96:99)                 # Mon_bm
Mac.mc.id <- c(43,                    # liver
               54, 55, 60, 61, 62,    # gut
               45, 63, 64, 65, 66, 67, 68, # kidney
               50, 51, 52, 53, 33, 34, 35, 36, 37, 39, 40, # skin
               44, 46, 47, 48, 56, 57, 58, 59) # Mac_Il1b
DC.mc.id <- seq(1:132)[mc@annots %in% c('cDC2_Mdh2','cDC1','mDC','Langerhans','cDC2_Ms4a4c',
                                        'cDC2_skin_Retnla(high)','cDC2_skin_Retnla(low)','pDC')]

## Fig. s2b: Mon subsets
marks <- c('H2-Ab1','H2-Aa','H2-Eb1','Cd74','Mir5107', # Ag presentation
           'Ear2','Ace','Treml4','Itgal','Itga4','Cebpb','Gngt2','Tnfrsf1b','Adgre4', # Mon
           'Mir1668','Cd300a','Clec4a1','Clec4a3','Nr4a1','Csf1r', # Mon
           'Mki67','Stmn1','Top2a','S100a8','Lcn2','S100a9','Hmgb2',
           'Fn1','Ccr2','Ms4a4c','Ms4a6c','Ms4a4a','Apoe','Ly6e','Ifitm6','S100a4') # core
tgconfig::set_param("mcp_heatmap_height", 40, "metacell")
tgconfig::set_param("mcp_heatmap_width", 10, "metacell")
tgconfig::set_param("mcp_heatmap_text_cex", 1, "metacell")
mcell_mc_plot_marks2(paste0(id,'_filt'), gset_id = id, mat_id = id,
    plot_cells = FALSE, mc_ord = Mon.mc.id, focus_mcs = Mon.mc.id,
    gene_list = marks,
    fig_fn = sprintf('%s/mouse_v9.mc_heat_marks_Mon_selected.pdf', figpath),
    zero_median = FALSE)

## Fig. s2c: Mac subsets
marks <- c('C1qa','C1qb','C1qc','Mrc1','Selenop','Apoe','Clptm1','Lgmn','Mafb','Csf1r','Cd81','C3ar1','Ctsb','Ctsc','Pf4','Dab2','Ms4a7', # core
           'H2-Ab1','H2-Eb1','Cd74','H2-Aa',
           'Alb','Clec4f','Marco','Vsig4','Cd5l','Ear2','Apoc1','Slc40a1', # liver
           'Mmp9','Dnase1l3','AW112010','Ptgs1','Pla2g7','Gbp2b','Adamdec1','Jaml','Acp5','Apol7c', # gut
           'Gimap3','Ly86','B2m','Slamf9','Hexb','Mmp12','Tmem176a','Tmem176b', # kidney
           'Clec10a','Lyve1','Cbr2','Fl3a1','Vim','Ccl2','Ccl7','Ccl8','Ccl9','Cxcl2','Retnla','Folr2','Crip1','Emp3','S100a4', # skin
           'Ccl3','Ccl4','Egr1','Il1b','Btg2','Nfkbiz','Ier5','Ier2','Fos','Atf3' # Il1b
          )
mcell_mc_plot_marks2(paste0(id,'_filt'), gset_id = id, mat_id = id,
    plot_cells = FALSE, mc_ord = Mac.mc.id, focus_mcs = Mac.mc.id,
    gene_list = marks,
    fig_fn = sprintf('%s/mouse_v9.mc_heat_marks_Mac_tissue_specfic_selected.pdf', figpath),
    zero_median = FALSE)

## Fig. s2d: Neu subsets
marks <- c('Stfa1','Stfa2','Stfa3','BC100530','Stfa2l1','Tnbs1','Gm5416','Gm10913', # Stfa2
           'Cxcl2','Egr1','Fgl2','Il1b','Btg2','Trem1','Nfkbia','Jun','Junb','Cxcr2', # Cxcl2
           'Sell','H2-Q10','Csf3r','Ccl6','Ly6g','C3','Chil3','Ltf','Ly6c2','Ngp','Camp','Cd177', # Ccl6
           'S100a8','S100a9','Mmp8','Fpr2','Cd24a','Wfdc21') # core
tgconfig::set_param("mcp_heatmap_height", 30, "metacell")
mcell_mc_plot_marks2(paste0(id,'_filt'), gset_id = id, mat_id = id,
    plot_cells = FALSE, mc_ord = Neu.mc.id, focus_mcs = Neu.mc.id,
    gene_list = marks,
    fig_fn = sprintf('%s/mouse_v9.mc_heat_marks_Neu_selected_markers.pdf', figpath),
    zero_median = FALSE)

## Fig. s2e: DC subsets
marks <- c('Itgax','H2-Aa','H2-Eb1','H2-Ab1','H2-Dmb1','H2-Dma','Cd74','Mir5107', # core
           'Mdh2','Ptgs2','Dock10', # Mdh2
           'Clec9a','Xcr1','Ppt1','Ifi205','Naaa', # cDC1
           'Fscn1','Ccr7','Il4i1','Cacnb3','Ccl22','Tmem123','Basp1','Ccl5','Tspan3', # mDC
           'Rgs2','Hpgds','Epcam','Mfge8','H2-M2','Cd207','Gpr171', # Langerhans
           'Fn1','Lpl','F13a1','Retnla','Lyz1','Ear2','Ccr2', # Retnla high
           'Ms4a4c','Cd209a','Apod','Cd209e', # Ms4a4c
           'Gpr171','Irf8','Rnase6','Plac8','Ly6e','Ccnd1','Mpeg1','Cd300c',
           'Ly6c2','Bst2','Spib','Runx2','Sell','Ly6a','Iglj3','Klk1','Ctsl',
           'Bst2','Atp1b1','Siglech','Cox6a2','Ly6d','Lair1','Ccr9','Cyb561a3','Tcf4','Cd7' # pDC
           )
mcell_mc_plot_marks2(paste0(id,'_filt'), gset_id = id, mat_id = id,
    plot_cells = FALSE, mc_ord = DC.mc.id, focus_mcs = DC.mc.id,
    gene_list = marks,
    fig_fn = sprintf('%s/mouse_v9.mc_heat_marks_DC_selected_markers.pdf', figpath),
    zero_median = FALSE)

## annotation bars for Figs. s2b-s2e (Group/Tissue/Time composition + cell type per metacell)
## example shown for DC; repeat by replacing DC.mc.id with Mon/Mac/Neu.mc.id
mc_info_os <- read.csv('mye/output/seurat/mye_366641c_mc_info_oversample.csv', row.names = 1)
Group.color <- setNames(c('#EDBF36','#7CCCE7'), c('SPF','GF'))
Tissue.color <- setNames(c('#F3CD3E','#6F93CC','#8A9095','#CEAC80','#ED7248','#CFDD9F','#2277B3','#E25B3D','#437034','#71C6E8'),
                         c('bm','colon','heart','kidney','liver','lung','si','skin','spleen','thymus'))
Time.color <- c('#56B4E9', '#E69F00', '#D55E00')

plot_anno_bars <- function(mc_ids, prefix){
    require(aplot)
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
    ggsave(sprintf('mye/output/mc1/heatmap/%s_specific_pct_bar.pdf', prefix), width = 8, height = 4)

    ## celltype annotation
    celltype_bar <- data.frame(id = 1:length(mc@annots[mc_ids]), celltype = mc@annots[mc_ids]) %>%
        mutate(celltype = factor(celltype)) %>%
    ggplot(aes(id, 1, fill = celltype)) +
        geom_tile() +
        scale_fill_manual(values = col.24) +
        theme_void()
    ggsave(sprintf('mye/output/mc1/heatmap/%s_celltype_bar.pdf', prefix), width = 8, height = 1)
}

plot_anno_bars(Mon.mc.id, 'Mon')
plot_anno_bars(Mac.mc.id, 'Mac')
plot_anno_bars(Neu.mc.id, 'Neu')
plot_anno_bars(DC.mc.id, 'DC')

#################################################################
## Figs. s2f-s2i
## Frequencies of selected myeloid subsets across tissues and stages in SPF vs GF mice
## shaded line plots show mean with SEM, each dot represents one mouse

## prepare per-mouse frequency table (ggprop_d)
mc_info %>% dplyr::select(Group, annots, Time, Tissue_only, Mouse) %>%
       mutate(comb = sprintf('%s_%s_%s_%s', Group, Time, Tissue_only, Mouse)) %>%
       dplyr::select(-Time, -Tissue_only, -Group, -Mouse) -> tmp

ggumap_ct <- table(tmp$annots, tmp$comb) %>% as.data.frame() %>%
  dplyr::rename(annots = Var1, comb = Var2, n = Freq)
## remove replicates with cell numbers less than 20
ggumap_ct %>% dplyr::group_by(comb) %>% filter(sum(n) >= 20) %>%
  transmute(Fraction = 100*n/sum(n), annots = annots) -> ggprop_d

ggprop_d$Group = sapply(ggprop_d$comb, function(x) str_split(x,'_')[[1]][1])
ggprop_d$Time = sapply(ggprop_d$comb, function(x) str_split(x,'_')[[1]][2])
ggprop_d$Tissue_only = sapply(ggprop_d$comb, function(x) str_split(x,'_')[[1]][3])
ggprop_d$Mouse = sapply(ggprop_d$comb, function(x) str_split(x,'_')[[1]][4])

pd <- position_dodge(width = 0)

## line layers for shaded line plots (mean +- SEM)
line_layers <- function() {
  list(
    stat_summary(
      fun.data = mean_se,
      geom = "ribbon",
      aes(group = Group, fill = Group),
      alpha = 0.15,
      color = NA
    ),
    geom_boxplot(
      width = 0.8,
      alpha = 0.2,
      position = pd,
      outlier.shape = NA
    ),
    geom_jitter(
      size = 2,
      shape = 21,
      fill = NA,
      position = pd
    ),
    # mean line
    stat_summary(
      fun = mean,
      geom = "line",
      aes(group = Group)
    ),
    scale_fill_manual(values = c('#77BBD1', '#DAB33D')),
    scale_color_manual(values = c('#77BBD1', '#DAB33D')),
    # statistics: SPF vs GF within each time point
    stat_compare_means(
      aes(group = Group),
      method = "t.test",
      label = "p",
      hide.ns = TRUE,
      label.y.npc = 0.95,
      tip.length = 0.02,
      bracket.size = 0.5
    ),
    facet_wrap(~ Tissue_only, scales = 'free_y', ncol = 8),
    theme_classic(),
    theme(
      axis.text.x = element_text(angle = 0, size = 10),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold"),
      legend.position = "None"
    ),
    ylab('Freq (%)'),
    xlab('Time')
  )
}

## Fig. s2f: Mac_skin_Retnla (labelled Mac_skin_Folr2 in the figure) in skin at P16 and W4
ggprop_d %>%
  dplyr::filter(Tissue_only %in% c('skin')) %>%
  dplyr::filter(Time %in% c('P16','W4')) %>%
  dplyr::filter(annots %in% c('Mac_skin_Retnla')) %>%
  mutate(
    Time = factor(Time),
    Group = factor(Group)
  ) %>%
  ggplot(aes(x = Time, y = Fraction, fill = Group, color = Group)) +
    line_layers()
ggsave('figures/Mac_skin_line_bar.pdf', width = 3.2, height = 2.5)

## Fig. s2g: Mon in bone marrow, kidney, liver and lung at W4 and W9
ggprop_d %>%
  dplyr::filter(Tissue_only %in% c('bm','lung','kidney','liver')) %>%
  dplyr::filter(annots %in% c('Mon')) %>%
  mutate(
    Time = factor(Time),
    Group = factor(Group)
  ) %>%
  ggplot(aes(x = Time, y = Fraction, fill = Group, color = Group)) +
    line_layers()
ggsave('figures/Mon_line_bar.pdf', width = 6.4, height = 2.5)

## Fig. s2h: Mac_Kidney in kidney at P16, W4 and W9
ggprop_d %>%
  dplyr::filter(Tissue_only %in% c('kidney')) %>%
  dplyr::filter(annots %in% c('Mac_Kidney')) %>%
  mutate(
    Time = factor(Time),
    Group = factor(Group)
  ) %>%
  ggplot(aes(x = Time, y = Fraction, fill = Group, color = Group)) +
    line_layers()
ggsave('figures/Mac_Kidney_line_bar.pdf', width = 1.6, height = 2.5)

## Fig. s2i: Neu_Il1b_Cxcl2 across indicated tissues
ggprop_d %>%
  dplyr::filter(Tissue_only %in% c('bm','heart','kidney','lung','liver')) %>%
  dplyr::filter(annots %in% c('Neu_Il1b_Cxcl2')) %>%
  mutate(
    Time = factor(Time),
    Group = factor(Group)
  ) %>%
  ggplot(aes(x = Time, y = Fraction, fill = Group, color = Group)) +
    line_layers()
ggsave('figures/Neu_Il1b_Cxcl2_line_bar.pdf', width = 8, height = 2.5)

## Fig. s2i: Neu_bm in bone marrow
ggprop_d %>%
  dplyr::filter(Tissue_only %in% c('bm')) %>%
  dplyr::filter(annots %in% c('Neu_bm')) %>%
  mutate(
    Time = factor(Time),
    Group = factor(Group)
  ) %>%
  ggplot(aes(x = Time, y = Fraction, fill = Group, color = Group)) +
    line_layers()
ggsave('figures/Neu_bm_line_bar.pdf', width = 1.6, height = 2.5)

## Fig. s2i: Neu_bm_Ccl6(Neg) in bone marrow
ggprop_d %>%
  dplyr::filter(Tissue_only %in% c('bm')) %>%
  dplyr::filter(annots %in% c('Neu_bm_Ccl6(Neg)')) %>%
  mutate(
    Time = factor(Time),
    Group = factor(Group)
  ) %>%
  ggplot(aes(x = Time, y = Fraction, fill = Group, color = Group)) +
    line_layers()
ggsave('figures/Neu_bm_Ccl6(Neg)_line_bar.pdf', width = 1.6, height = 2.5)
