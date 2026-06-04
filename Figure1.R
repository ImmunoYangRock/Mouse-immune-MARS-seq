#################################################################
## R code to create figures related to Figure1
## Figs 1c, 1d, 1e, 1f, 1g and Figs s1b, s1c, s1f, s1g, s1i
#################################################################

## load data
seurat <- readRDS('all/data/mousev9_20230407_120649c.rds')
annot_d = c('B_pro','B_pre','B_mature','plasma','T_pre','T_naive','T_mature','T_gd','ILC2','NK','Mon','Mac','KupfferCell','Langerhans','cDC2','cDC1','pDC','mDC','Neu','Basophil','Eos','MastCell')
seurat$annot_d <- factor(seurat$annot_d, level = annot_d)
## select columns to plot
ggumap = data.frame(seu@reductions$umap@cell.embeddings,
                    annot_r = factor(seu$annot_r), 
                    annot_d = factor(seu$annot_d),
                    Tissue = seu$Tissue_only,
                    Group = seu$Group, 
                    Time = seu$Time, 
                    Mouse = seu$Mouse, 
                    Plate = seu$Plate,
                    Phrase = seu$Phase, 
                    Gate = seu$Gate)
#################################################################
## Fig. 1c Cell-type specific DotPlot with bar showing time/tissue/group 
gl <- c('Igll1','Vpreb1','Vpreb2', # B_pro # 'Bach2','Cd93','Tcf3',
'Cd37','Ly6d', # B_pre
'Ms4a1','Fcmr','Scd1','Dntt', # B_mature        
'Jchain','Igha','Sdc1',# Plasma
'Rag1','Ccr9', # T_pre
'Lef1','S1pr1', # T_naive
'Ccl5','Gzmb','Gzma','Gzmk', # T_mature
'Trdc','Tcrg-C1', # T_gd
'Cd3d',
'Rora','Cxcr6','Gata3', # ILC
'Gzmc','Xcl1','Klrc1', # NK
'Ms4a6c','Clec4a3','Sirpb1c', # Mon
'Pf4','C1qa','C1qb','C1qc', # Mac
'Cd5l','Clec4f','Vsig4', #KupfferCell
'Cd207','H2-M2','Epcam', # Langerhans
'Cd209a','Clec4b1','H2-DMb1', # cDC2
'Xcr1','Clec9a', # cDC1
'Siglech','Cox6a2', # pDC
'Fscn1','Cacnb3', # mDC
'Retnlg','S100a8','S100a9', # Neu
'Mcpt8','Car1', # Basophil
'Ear1','Car4','Ear2',  #Eos
'Mcpt4','Cma1','Cpa3') # MastCell

## dotplot
plot.Dot <- DotPlot(seurat,group.by = 'annot_d', 
        features = gl,cols = c("white", "red"),
        dot.min = 0, cluster.idents = F) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.direction = "vertical", #horizontal
        legend.position = "bottom")

## bar, time
time.pct <- ggumap %>% 
  dplyr::select(Time,annot_d) %>% 
  group_by(Time,annot_d) %>% 
  count() %>%
ggplot(aes(annot_d, freq, fill = Time)) + 
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values = c('#bcbddcff','#9e9ac8ff','#756bb1ff'))+coord_flip()+theme_void()+
  theme(legend.text = element_text(size=12),legend.title =element_text(size=12) )+
  guides(fill = guide_legend(title = "Time", override.aes = list(size = 5), ncol = 1)) 

## bar, tissue
tissue.pct <-ggumap %>% 
  dplyr::select(Tissue,annot_d) %>% 
  group_by(Tissue,annot_d) %>% 
  count() %>%
ggplot(aes(annot_d, freq, fill = Tissue)) + 
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values =pal_simpsons()(10) )+coord_flip()+theme_void()+
  theme(legend.text = element_text(size=12),legend.title =element_text(size=12) )+
  guides(fill = guide_legend(title = "Tissue", override.aes = list(size = 5), ncol = 1)) 

## bar, group
group.pct <-ggumap %>% 
  dplyr::select(Group,annot_d) %>% 
  group_by(Group,annot_d) %>% 
  count() %>%
ggplot(aes(annot_d, freq, fill = Group)) + 
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_manual(values =c('#6ee2ffff','#f7c530ff') )+coord_flip()+theme_void()+
  theme(legend.text = element_text(size=12),legend.title =element_text(size=12) )+
  guides(fill = guide_legend(title = "Group", override.aes = list(size = 5), ncol = 1)) 

## merge
plot.Dot %>% 
  insert_left(time.pct,width = .05) %>% 
  insert_left(group.pct,width=.05) %>% 
  insert_right(tissue.pct,width=.08) 
ggsave('all/output/seurat/marker_gene_plot.pdf',dotplot,width=20,height=8)

#################################################################
## Fig1d-Fig1f
## bar plot to show celltype frequenceis across time-tissue among SPF mice
options(repr.plot.width = 20, repr.plot.height=10)
ggumap %>% dplyr::select(annot_d, Group, Time, Tissue) %>% 
dplyr::group_by(annot_d,Tissue, Time, Group) %>% count(name='freq') %>%
ggplot(aes(Time, freq, fill = annot_d)) + geom_bar(stat = 'identity', position='fill') +
scale_fill_manual(values = col.22) + theme_void()+ #coord_flip()+ 
theme(axis.text = element_text(size=10), legend.title = element_blank(),
      legend.text =element_text(size=10), 
      strip.text =element_text(size=15),
      legend.position = 'bottom')+
guides(fill = guide_legend(title = "", ncol=22)) +
facet_grid(rows = vars(factor(Group,level=c('SPF','GF'))), cols = vars(factor(Tissue,level=c('spleen','bm','liver','thymus','lung',
                                                                                             'si','colon','heart','kidney','skin')))) -> barplot
ggsave('all/output/seurat/barplot.pdf',barplot,width=12,height=4)

#################################################################
## Fig 1e-Fig1g
## dot plot to show time/group variation
custom_theme <- theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5, vjust=0.5),
    axis.ticks.x = element_blank()
  )

## Firstly, prepare a summary table ggtable_freq
ggumap %>% 
    dplyr::select(Group, Time, Tissue, annot_d) %>%  
    mutate(comb = sprintf('%s_%s_%s',Group, Time, Tissue)) %>%    
    dplyr::select(-Group,-Time,-Tissue) -> tmp
## count
table(tmp$annot_d,tmp$comb) %>% as.data.frame() %>% 
    dplyr::rename(annot_d = Var1, comb = Var2, n = Freq) -> ggtable_cnt
## freq
ggtable_cnt %>% 
    dplyr::group_by(comb) %>% 
    transmute(prop = n/sum(n), annot_d = annot_d) -> ggtable_freq
## reshape
ggtable_freq <- ggtable_freq %>% tidyr::separate(comb, into = c("Group", "Time", "Tissue"),sep = "_")

#################################################################
## Fig 1e
## W4/P16, also W9/W4, across Group and Tissue
## P16 meanprop
ggtable_freq %>%
    filter(Time=='P16') %>%
    select(-Time) %>%
    group_by(Group, Tissue, annot_d) %>%
    summarise(mean_prop = mean(prop)) -> P16meanProp #(P16-> Group,Tissue,annot_d,mean_prop)

## W4 meanprop
ggtable_freq %>%
    filter(Time=='W4') %>%
    select(-Time) %>%
    group_by(Group, Tissue, annot_d) %>%
    summarise(mean_prop = mean(prop)) -> W4meanProp #(W4-> Group,Tissue,annots,mean_prop)

## W4 query on P16
ggtable_freq %>%
    filter(Time=='W4') %>%
    select(-Time) -> W4_query #(W4, Group, Tissue,annots,prop)

## W9 query on W4
ggtable_freq %>%
    filter(Time=='W9') %>%
    select(-Time) -> W9_query #(W4, Group, Tissue,annots,prop)

## W4/P16
merge(W4_query,P16meanProp, by=c('Group','Tissue','annot_d')) %>%
    mutate(abs_solo = abs(prop-mean_prop)) %>%
    select(-prop,-mean_prop) %>%
    group_by(Group,Tissue,Mouse) %>%
    summarise(abs_total = sum(abs_solo)) %>%
    ungroup() -> W4_P16_abs_df

# W9/W4
merge(W9_query,W4meanProp, by=c('Group','Tissue','annot_d')) %>%
    mutate(abs_solo = abs(prop-mean_prop)) %>%
    select(-prop,-mean_prop) %>%
    group_by(Group,Tissue,Mouse) %>%
    summarise(abs_total = sum(abs_solo)) %>%
  ungroup() -> W9_W4_abs_df

# W9/W4 & W4/P16 only SPF
W9_W4_abs_df$stageA2B = 'W9/W4'
W4_P16_abs_df$stageA2B = 'W4/P16'

SPF_df <- rbind(W9_W4_abs_df,W4_P16_abs_df) %>%
mutate(stageA2B = factor(stageA2B,levels=c('W4/P16','W9/W4'))) %>%
filter(Group=='SPF') %>% as.data.frame() 

means <- SPF_df %>%
  group_by(Tissue,stageA2B) %>%
  summarize(mean_abs_total = mean(abs_total))

SPF_df %>% 
  #filter(Tissue!='skin') %>%
ggplot(aes(stageA2B,abs_total,fill=stageA2B)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5) + # fun.data = mean_se, geom = "errorbar", width = 0.2
  geom_jitter(shape=21,size=2,width = 0.25) +
  facet_wrap(~factor(Tissue,levels=c('si','heart','kidney','colon','liver','bm','lung','spleen','skin','thymus')),ncol=10)+
  scale_fill_manual(values=c('#374F9A','#B74791'))+
  custom_theme+
  theme(strip.text.x = element_text(size=8),legend.position = 'None')
ggsave('figures/v9_spf_all_cell_freq_abs_stageA2B.pdf',width=6,height=2)


#################################################################
## Fig 1g
## condition
Time <- c('P16','W4','W9')
Tissue <- c('bm','thymus', 'spleen', 'liver', 'lung', 'si', 'colon', 'kidney','heart','skin')
Annots <- unique(ggtable_freq$annot_d)
## build table frame
ggabs <- ggtable_freq %>% 
    ungroup() %>% 
    dplyr::select(Time,Tissue) %>% 
    distinct(Time,Tissue) %>%
    mutate(dif = numeric(n()))

## loop
for (time in Time){
    df_time = ggtable_freq[ggtable_freq$Time==time,]
    for (tissue in Tissue) {
        df_time_tissue = df_time[df_time$Tissue==tissue,]
        ## Sum up the differences corresponding to each time and each organization
        abs_total = 0
        for (annot in Annots) {
            df_time_tissue_annots = df_time_tissue[df_time_tissue$annot_d==annot,]
            if (nrow(df_time_tissue_annots) > 0) {
                ##  verify two groups(SPF and GF) all existed
                spf_value = df_time_tissue_annots[df_time_tissue_annots$Group == 'SPF', 'prop']
                gf_value = df_time_tissue_annots[df_time_tissue_annots$Group == 'GF', 'prop']
                if (length(spf_value) > 0 && length(gf_value) > 0) {
                    abs_solo = abs(spf_value - gf_value)
                    abs_total = abs_total + abs_solo
                }
            }
        }
        message(time,tissue)
        ggabs$dif[ggabs$Time == time & ggabs$Tissue == tissue] <- abs_total$prop
    }
}

## mutate a col to save all var cross-combinations
ggumap %>% 
    dplyr::select(Group, Time, Mouse, Tissue, annot_d) %>%  
    mutate(comb = sprintf('%s_%s_%s_%s', Group, Time, Tissue, Mouse)) %>%  
    dplyr::select(-Time, -Tissue, -Group, -Mouse) -> tmp

## generate cross count table
table(tmp$annot_d, tmp$comb) %>% 
    as.data.frame() %>% 
    dplyr::rename(annot_d = Var1, comb = Var2, n = Freq) -> ggtable_cnt

## calculate proportions
ggtable_cnt %>% 
    dplyr::group_by(comb) %>%
    transmute(prop = n/sum(n), annot_d = annot_d) %>%
    ungroup() -> ggtable_freq

## split variables from the crossed combinations
ggtable_freq <- ggtable_freq %>%
    tidyr::separate(comb, into = c("Group", "Time", "Tissue", "Mouse"), sep = "_")      

## calculate SPF-Mean                               
ggtable_freq %>%
    filter(Group=='SPF') %>%
    select(-Group) %>%
    group_by(Time, Tissue, annot_d) %>%
    summarise(SPF_prop = mean(prop),.groups = 'drop') -> SPFmeanProp #(Time,Tissue,annots,SPF_prop)
                            
## extract GF-duplicates                            
ggtable_freq %>%
    filter(Group=='GF') %>%
    select(-Group) -> GF_query #(Time,Tissue,annots,prop)
                            
## merge 
options(repr.plot.width =10, repr.plot.height=4)
merge(GF_query,SPFmeanProp, by=c('Time','Tissue','annot_d'))   %>%
    mutate(abs_solo = abs(prop-SPF_prop)) %>% 
    select(-prop,-SPF_prop) %>%
    group_by(Time,Tissue,Mouse) %>%
    summarise(abs_total = sum(abs_solo),.groups = 'drop') %>%
    ungroup() %>%
ggplot(aes(Time,abs_total,fill=Time)) + 
    facet_wrap(~factor(Tissue,levels=c('si','heart','skin','kidney','liver','colon','lung','bm','spleen','thymus')),ncol = 10)+
    geom_jitter(shape=21,size=2,width = 0.2) + 
    scale_fill_manual(values=c('#56B4E9', '#E69F00', '#D55E00'))+
    stat_summary(fun = mean, geom = "crossbar", width = 0.5) + # fun.data = mean_se, geom = "errorbar", width = 0.2                            
    custom_theme+
    theme(strip.text.x = element_text(size=8)) + 
    Seurat::NoLegend()
ggsave('figures/v9_all_time_freq_abs_SPF2GF.pdf',width=6,height=2)          

#################################################################
## Fig s1b
## bar plot to show cell counts across group-time-tissue
SPF.freq <- ggumap %>% 
  filter(Group == 'SPF') %>% 
  dplyr::select( Time, Tissue) %>% 
  table() %>% 
  as.data.frame.matrix() %>% rownames_to_column('Time') %>% 
  reshape2::melt(variable.name = 'Tissue')

GF.freq <- ggumap %>% 
  filter(Group == 'GF') %>% 
  dplyr::select( Time, Tissue) %>% 
  table() %>% 
  as.data.frame.matrix() %>% rownames_to_column('Time') %>% 
  reshape2::melt(variable.name = 'Tissue')

### X axis: time, Y axis: tissue
Time.level <- c('P16','W4','W9')
Tissue.level <- c('bm','thymus', 'spleen', 'liver', 'lung', 'si', 'colon', 'skin', 'kidney','heart')
SPF.bar <- ggplot(SPF.freq, aes(x=factor(Time, level = Time.level), y=value, fill = Time)) + 
  scale_fill_manual(values = c('#87ceeb','#48d1cc','#008080'))+
  geom_bar(stat = 'identity')+ 
  ylim(0,5000)+
  geom_text(aes(x = factor(Time), y = value+200, label = value), size = 0.8, fontface = "bold") +
  theme_void()+
  guides(fill = guide_legend(title = "", override.aes = list(size = 3)))+
  theme(axis.text.x = element_blank(),
        strip.text = element_blank(),
       legend.position = "None")+
  facet_wrap(factor(Tissue,level = Tissue.level)~.,ncol = 10,strip.position="left")

GF.bar <- ggplot(GF.freq, aes(x=factor(Time, level = Time.level), y=value, fill = Time)) + 
  scale_fill_manual(values = c('#87ceeb','#48d1cc','#008080'))+
  geom_bar(stat = 'identity')+  
  ylim(0,5000)+
  geom_text(aes(x = factor(Time), y = value+200, label = value), size = 0.8, fontface = "bold") +
  theme_void()+
  guides(fill = guide_legend(title = "", override.aes = list(size = 3)))+
  theme(axis.text = element_blank(),
        strip.text = element_blank(),
       legend.position = "None")+
  facet_wrap(factor(Tissue,level = Tissue.level)~.,ncol = 10,strip.position="left")

options(repr.plot.width = 6, repr.plot.height=4)
require(aplot)
SPF.bar %>% insert_bottom(GF.bar,height=1)
ggsave('all/output/seurat/barplot/cell_count_barplot.pdf',width=3,height=1.5)


#################################################################
## Fig s1c
## qc checek
v9.qc.data = data.frame(nCount_RNA = Matrix::colSums(seu@assays$RNA@counts),
                       nFeature_RNA = Matrix::colSums(seu@assays$RNA@counts>0),
                       percent_mt = Matrix::colSums(seu@assays$RNA@counts[grep('^mt-', rownames(seu), value=TRUE),])/
                                                    Matrix::colSums(seu@assays$RNA@counts),
                       Group = seu@meta.data$Group)
ggpubr::ggviolin(v9.qc.data, x='Group',y = c('percent_mt'), color = 'Group',
         palette = c('royalblue','orangered'),add = 'boxplot',add.params = list(fill = 'white'))
ggsave('all/output/seurat/qc_plot/v9_mt.pdf',width=4,height=2)

ggpubr::ggviolin(v9.qc.data[v9.qc.data$nCount_RNA<=10000,], x='Group',y = c('nCount_RNA'), color = 'Group',
         palette = c('royalblue','orangered'), add = 'boxplot',trim = TRUE, add.params = list(fill = 'white'))
ggsave('all/output/seurat/qc_plot/v9_count.pdf',width=4,height=2)

ggpubr::ggviolin(v9.qc.data[v9.qc.data$nCount_RNA<=10000,], x='Group',y = c('nFeature_RNA'), color = 'Group',
         palette = c('royalblue','orangered'),add = 'boxplot',trim = TRUE,add.params = list(fill = 'white'))
ggsave('all/output/seurat/qc_plot/v9_feature.pdf',width=4,height=2)


#################################################################
## Fig s1f, s1g, s1i
## statistical charts

## prepare dataframe to plot, ggprop_d, ggprop_r
ggumap %>% 
  dplyr::select(Group, annot_d, Time, Tissue, Mouse) %>%  
  mutate(Mouse = sprintf('%s_%s_%s',Group, Time, Mouse)) %>%  
  dplyr::select(-Time) -> tmp
ggprop_d <- table(tmp$annot_d,tmp$Mouse,tmp$Tissue) %>% 
  as.data.frame() %>% 
  dplyr::rename(annot_d = Var1, Mouse = Var2, Tissue = Var3, n = Freq)  %>%
  dplyr::group_by(Mouse,Tissue) %>% 
  transmute(prop = 100* n/sum(n), annot_d = annot_d) 
ggprop_d$Group = sapply(ggprop_d$Mouse, function(x) str_split(x,'_')[[1]][1]) 
ggprop_d$Time = sapply(ggprop_d$Mouse, function(x) str_split(x,'_')[[1]][2]) 

ggumap %>% 
  dplyr::select(Group, annot_r, Time, Tissue, Mouse) %>%  
  mutate(Mouse = sprintf('%s_%s_%s',Group, Time, Mouse)) %>%  
  dplyr::select(-Time) -> tmp
ggprop_r <- table(tmp$annot_r,tmp$Mouse,tmp$Tissue) %>% 
  as.data.frame() %>% 
  dplyr::rename(annot_r = Var1, Mouse = Var2, Tissue = Var3, n = Freq)  %>%
  dplyr::group_by(Mouse,Tissue) %>% 
  transmute(prop =100* n/sum(n), annot_r = annot_r) 
ggprop_r$Group = sapply(ggprop_r$Mouse, function(x) str_split(x,'_')[[1]][1])
ggprop_r$Time = sapply(ggprop_r$Mouse, function(x) str_split(x,'_')[[1]][2]) 

## bar_theme
bar_layers <- function() {
  list(
    geom_bar(stat = "summary", fun = "mean", width = 0.5, position = position_dodge()),
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2),
    scale_fill_manual(values = c('#dadaebff','#9e9ac8ff','#756bb1ff')),
    geom_jitter(shape = 21, size = 1.5, width = 0.2, stroke = 0.2, fill = "white"),
    theme_classic(),
    theme(axis.text.x = element_blank(), panel.background = element_blank(), plot.background = element_blank()),
    ylab('Freq(%)'),
    xlab(''),
    NoLegend()
  )
}

#################################################################
## Fig s1f      
## si&colon: T_mature, B_mature                        
ggprop_d %>% na.omit %>% dplyr::filter(Tissue %in% c('si','colon') & annot_d %in% c('T_mature','B_mature')) %>%
mutate(annot_d = factor(annot_d,levels = c('T_mature','B_mature'))) %>%
ggplot(aes(x =  paste0(Group,Time), y = prop,fill=Time)) + 
 bar_layers()+
 ggsignif::geom_signif(comparisons = list(c("GFP16", "GFW4"), c("GFW4", "GFW9"),
                                          c("SPFP16", "SPFW4"), c("SPFW4", "SPFW9"), 
                                          c("SPFP16", "GFP16"), c("SPFW4", "GFW4"),  c("SPFW9", "GFW9")), 
         margin_top = 0.05, test = "t.test", map_signif_level = F,  family = "serif", textsize = 2.5,vjust=1.5, step_increase = 0.05,
         test.args = list(alternative = "two.side", var.equal = F, paired=F)) + 
 facet_grid(rows = vars(Tissue), cols = vars(annot_d), scales='free') 
ggsave('figures/B_T_mature_gut_group_time_freq.pdf',width=3.5, height=4)  


#################################################################
## Fig s1g
## bm: T_naive, T_mature, nK, B_mature                      
ggprop_d %>% na.omit %>% dplyr::filter(Tissue == 'bm' & annot_d %in% c('T_naive','T_mature','NK','B_mature')) %>%
mutate(annot_d = factor(annot_d, levels = c('NK','T_naive','T_mature','B_mature'))) %>%
ggplot(aes(x = paste0(Group,Time), y = prop, fill=Time)) + bar_layers() + 
 ggsignif::geom_signif(comparisons = list(c("GFP16", "GFW4"), c("GFW4", "GFW9"),
                                          c("SPFP16", "SPFW4"), c("SPFW4", "SPFW9"), 
                                          c("SPFP16", "GFP16"), c("SPFW4", "GFW4"),  c("SPFW9", "GFW9")), 
         margin_top = 0.05, test = "t.test", map_signif_level = F,  family = "serif", textsize = 1.5, step_increase = 0.15, 
         test.args = list(alternative = "two.side", var.equal = F, paired=F)) + 
 facet_wrap(~annot_d,nrow=1,scales='free')
ggsave('figures/BM_lym_group_time_freq.pdf',width=6,height=2)


#################################################################
## Fig s1i 
## bm: cDC1, cDC2, pDC                       
ggprop_d %>% na.omit %>% dplyr::filter(Tissue == 'bm' & annot_d %in% c('cDC1','cDC2','pDC')) %>%
ggplot(aes(x = paste0(Group,Time), y = prop, fill=Time)) + 
 bar_layers()+
 ggsignif::geom_signif(comparisons = list(c("GFP16", "GFW4"), c("GFW4", "GFW9"),
                                          c("SPFP16", "SPFW4"), c("SPFW4", "SPFW9"), 
                                          c("SPFP16", "GFP16"), c("SPFW4", "GFW4"),  c("SPFW9", "GFW9")), 
         margin_top = 0.05, test = "t.test", map_signif_level = F,  family = "serif", textsize = 1.5, step_increase = 0.15, 
         test.args = list(alternative = "two.side", var.equal = F, paired=F)) + 
 facet_wrap(~annot_d,scales='free')
ggsave('figures/BM_DC_group_time_freq.pdf',width=4.5,height=2)                       
                       
