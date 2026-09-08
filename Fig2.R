#################################################################
## R code to create figures related to Figure2
## Figs 2a-2h and Figs s2a-s2i
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





#################################################################
## Fig. 2c




#################################################################
## Fig. 2d




#################################################################
## Fig. 2e





#################################################################
## Fig. 2f




#################################################################
## Fig. 2g





#################################################################
## Fig. 2h



#################################################################
## Fig. s2b



#################################################################
## Fig. s2c



#################################################################
## Fig. s2d



#################################################################
## Fig. s2e



#################################################################
## Fig. s2g


#################################################################
## Fig. s2h


#################################################################
## Fig. s2i
