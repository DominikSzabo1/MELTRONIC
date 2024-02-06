library(tidyverse)
`%notin%` <- Negate(`%in%`)

#list files in directory and read in
score_files <- list.files('../../data/melting_scores/', pattern = 'n574.tsv$',full.names = TRUE)
score_files_tbl <- tibble(path=score_files) %>% 
  tidyr::separate(path, into=c('1', '2', '3', '4', '5' ,'comparison'),sep='/', remove=FALSE) %>% 
  dplyr::mutate(comparison = str_remove(comparison, '_long_genes_n574.tsv')) %>% 
  dplyr::select(path, comparison) %>% 
  dplyr::mutate(path_number = as.character(row_number()))
read_in <- map_dfr(score_files_tbl$path, read_tsv, .id = 'path_number') 

#modify read in files
scores_allcomp <- left_join(read_in, score_files_tbl) %>% 
  dplyr::select(-path, -path_number) %>% 
  dplyr::mutate(score = case_when(melting_score > 0 ~ melting_score,
                                  melting_score == 0 ~ -condensation_score)) %>% 
  tidyr::separate(comparison, into=c('DN', 'treatment_comp', 'timepoint_comp', 'genotype_comp', 'vs', 'DN2', 'treatment_ref', 'genotype_ref')) %>% 
  dplyr::mutate(comp_short = paste(treatment_comp, timepoint_comp, 'vs', treatment_ref, genotype_ref,  sep = '_')) %>% 
  dplyr::mutate(comp_short = factor(comp_short, levels = c('cocaine_1d_vs_saline_wt', 'cocaine_14dR1_vs_saline_wt', 'cocaine_14dR2_vs_saline_wt'))) %>%  
  dplyr::select(chrom, start, end, ID, perc_na_or_missing_ref, perc_na_or_missing_inp, comp_short, score) %>% 
  dplyr::rename(gene_id = ID, start_bin=start, end_bin=end)

#USE a melting /condensing threshold of 10 and -10 for the long genes
scores_allcomp_wide <- scores_allcomp %>% 
  dplyr::select(-perc_na_or_missing_ref, -perc_na_or_missing_inp) %>% 
  pivot_wider(names_from=comp_short, values_from = score, names_prefix = 'score_') %>% 
  dplyr::mutate(across(starts_with('score'), ~  case_when(. > 10 ~ 'm', 
                                                          .< -10 ~ 'c',
                                                          TRUE ~ 'n'),
                       .names='state_{.col}'),
                ms_dyn_R1=paste0(state_score_cocaine_1d_vs_saline_wt, state_score_cocaine_14dR1_vs_saline_wt),
                ms_dyn_R2=paste0(state_score_cocaine_1d_vs_saline_wt, state_score_cocaine_14dR2_vs_saline_wt),
                state_score_cocaine_1d_vs_saline_wt = factor(state_score_cocaine_1d_vs_saline_wt,levels=c('m','n','c')),
                state_score_cocaine_14dR1_vs_saline_wt = factor(state_score_cocaine_14dR1_vs_saline_wt,levels=c('m','n','c')),
                state_score_cocaine_14dR2_vs_saline_wt = factor(state_score_cocaine_14dR2_vs_saline_wt,levels=c('m','n','c')),
                ms_dyn_R1 = factor(ms_dyn_R1, levels=c('mm','nm','cm', 'mn','nn','cn','mc','nc','cc')),
                ms_dyn_R2 = factor(ms_dyn_R2, levels=c('mm','nm','cm', 'mn','nn','cn','mc','nc','cc'))) %>% 
  dplyr::mutate(state_repr = if_else(state_score_cocaine_14dR1_vs_saline_wt==state_score_cocaine_14dR2_vs_saline_wt, TRUE, FALSE))


#read in CPMs and DEG of long genes 
library(readxl)
deg_cpm <- read_xlsx('../../data/Table_S3_RNA_cpm_deg.xlsx', sheet = 2)

scores_allcomp_wide_repr <- left_join(scores_allcomp_wide, deg_cpm %>% dplyr::select(-chrom)) %>% 
  dplyr::filter(state_repr==TRUE)

#Fig 3A
library(ggrepel)
theme_set(theme_classic())

scores_allcomp_wide_toshow <- scores_allcomp_wide_repr %>% 
  dplyr::mutate(highlight_left = if_else(gene_name %in% c('Rbfox1', 'Ptprd', 'Il1rapl2', 'Cdh18', 'Csmd1', 'Lrrc7', 'Alk', 'Erbb4', 'Lrp1b',  'Adgrl3', 'Zfhx3', 'Rbms3', 'Grin2a', 'Syt1', 'Galntl6'), TRUE, FALSE))

scores_allcomp_wide_toshow %>% 
  ggplot()+
  geom_hline(yintercept=c(-10,10), color='grey50', linetype='dotted')+
  geom_point(aes(x='coc_1d_score', y=score_cocaine_1d_vs_saline_wt, color=state_score_cocaine_1d_vs_saline_wt), alpha=0.5, position=position_jitter(width=0.05, seed=2))+
  geom_point(aes(x='coc_14d_score', y=score_cocaine_14dR1_vs_saline_wt, color=state_score_cocaine_14dR1_vs_saline_wt), alpha=0.5, position=position_jitter(width=0.05, seed=2))+
  geom_segment(aes(x='coc_1d_score', xend='coc_14d_score', y=score_cocaine_1d_vs_saline_wt, yend=score_cocaine_14dR1_vs_saline_wt, color=state_score_cocaine_1d_vs_saline_wt, alpha = highlight_left), position = position_jitter(width=0.05, seed=2))+
  scale_alpha_discrete(c(0.6,1))+
  scale_x_discrete(limits = c('coc_1d_score', 'coc_14d_score')) +
  geom_text_repel(data = scores_allcomp_wide_toshow %>% dplyr::filter(gene_name %in% c('Rbfox1', 'Ptprd', 'Lrc7', 'Syt1', 'Grin2a', 'Alk', 'Cdh18', 'Il1rapl2')),aes(x='coc_14d_score', y=score_cocaine_14dR1_vs_saline_wt, label=gene_name), nudge_x = 0.2, min.segment.length = 10)+
  geom_text_repel(data = scores_allcomp_wide_toshow %>% dplyr::filter(gene_name %in% c('Lrp1b', 'Galntl6', 'Adgrl3', 'Csmd1', 'Erbb4', 'Rbms3', 'Zfhx3')),aes(x='coc_1d_score', y=score_cocaine_1d_vs_saline_wt, label=gene_name), nudge_x = -0.2, min.segment.length = 10)+
  scale_color_manual(values = c('n'='grey70', 'm'='#4256DE', 'c' = '#FF6699'))+
  theme(axis.title = element_blank(), legend.position = 'none')


#Fig 3B
scores_allcomp_wide_repr %>% 
  ggplot()+
  geom_point(aes(x=score_cocaine_1d_vs_saline_wt, y=CPM_sal_d01, color=state_score_cocaine_1d_vs_saline_wt))+
  scale_y_log10()+
  scale_color_manual(values=c('n'='grey70', 'm'='#4256DE', 'c' = '#FF6699'))+
  coord_cartesian(ylim=c(1, 10100))+
  geom_text_repel(data=scores_allcomp_wide_repr %>% dplyr::filter(gene_name %in% c('Rbfox1', 'Ptprd', 'Csmd1', 'Grin2a', 'Zfhx3', 'Galntl6', 'Alk', 'Lrp1b', 'Lsamp')), aes(x=score_cocaine_1d_vs_saline_wt, y=CPM_sal_d01, label=gene_name), min.segment.length = 0)+
  theme(legend.position = 'none')

medians_d01 <- scores_allcomp_wide_repr %>% 
  group_by(state_score_cocaine_1d_vs_saline_wt) %>%
  summarise(median=median(CPM_sal_d01, na.rm=TRUE))

scores_allcomp_wide_repr %>% 
  ggplot()+
  geom_density(aes(y=CPM_sal_d01, fill=state_score_cocaine_1d_vs_saline_wt), alpha=0.5)+
  scale_y_log10()+
  coord_cartesian(ylim=c(1, 10100))+
  geom_hline(data=medians_d01, aes(yintercept=median, color=state_score_cocaine_1d_vs_saline_wt), linewidth=1.5, color=c('m'='#4256DE', 'c' = '#FF6699', 'n'='grey70'), linetype='dashed')+ 
  scale_fill_manual(values=c('n'='grey70', 'm'='#4256DE', 'c' = '#FF6699'))

wilcox.test((scores_allcomp_wide_repr %>% dplyr::filter(state_score_cocaine_1d_vs_saline_wt=='n'))$CPM_sal_d01,
            (scores_allcomp_wide_repr %>% dplyr::filter(state_score_cocaine_1d_vs_saline_wt=='c'))$CPM_sal_d01)
wilcox.test((scores_allcomp_wide_repr %>% dplyr::filter(state_score_cocaine_1d_vs_saline_wt=='n'))$CPM_sal_d01,
            (scores_allcomp_wide_repr %>% dplyr::filter(state_score_cocaine_1d_vs_saline_wt=='m'))$CPM_sal_d01)

#fig S18A
scores_allcomp_wide_repr %>% 
  ggplot()+
  geom_point(aes(x=score_cocaine_14dR1_vs_saline_wt, y=CPM_sal_d14, color=state_score_cocaine_14dR1_vs_saline_wt))+
  scale_y_log10()+
  scale_color_manual(values=c('n'='grey70', 'm'='#4256DE', 'c' = '#FF6699'))+
  coord_cartesian(ylim=c(1, 10100))+
  geom_text_repel(data=scores_allcomp_wide_repr %>% dplyr::filter(gene_name %in% c('Rbfox1', 'Ptprd', 'Csmd1', 'Grin2a', 'Zfhx3', 'Galntl6', 'Alk', 'Lrp1b', 'Lsamp')), aes(x=score_cocaine_14dR1_vs_saline_wt, y=CPM_sal_d14, label=gene_name), min.segment.length = 0)+
  theme(legend.position = 'none')

medians_d14 <- scores_allcomp_wide_repr %>% 
  group_by(state_score_cocaine_14dR1_vs_saline_wt) %>%
  summarise(median=median(CPM_sal_d14, na.rm=TRUE))

scores_allcomp_wide_repr %>% 
  ggplot()+
  geom_density(aes(y=CPM_sal_d14, fill=state_score_cocaine_14dR1_vs_saline_wt), alpha=0.5)+
  scale_y_log10()+
  coord_cartesian(ylim=c(1, 10100))+
  geom_hline(data=medians_d14, aes(yintercept=median, color=state_score_cocaine_14dR1_vs_saline_wt), linewidth=1.5, color=c('m'='#4256DE', 'n'='grey70',  'c' = '#FF6699'), linetype='dashed')+ 
  scale_fill_manual(values=c('n'='grey70', 'm'='#4256DE', 'c' = '#FF6699'))

wilcox.test((scores_allcomp_wide_repr %>% dplyr::filter(state_score_cocaine_14dR1_vs_saline_wt=='n'))$CPM_sal_d14,
            (scores_allcomp_wide_repr %>% dplyr::filter(state_score_cocaine_14dR1_vs_saline_wt=='c'))$CPM_sal_d14)
wilcox.test((scores_allcomp_wide_repr %>% dplyr::filter(state_score_cocaine_14dR1_vs_saline_wt=='n'))$CPM_sal_d14,
            (scores_allcomp_wide_repr %>% dplyr::filter(state_score_cocaine_14dR1_vs_saline_wt=='m'))$CPM_sal_d14)

#fig S18B
scores_allcomp_wide_repr %>% group_by(state_score_cocaine_1d_vs_saline_wt) %>% summarise(median(logFC_d01))
scores_allcomp_wide_repr %>% 
  ggplot()+
  geom_vline(xintercept=-0.0164, color='grey50', linewidth=0.5)+
  geom_violin(aes(y=state_score_cocaine_1d_vs_saline_wt, x=logFC_d01, fill=state_score_cocaine_1d_vs_saline_wt))+
  geom_boxplot(aes(y=state_score_cocaine_1d_vs_saline_wt, x=logFC_d01), outlier.shape=NA, width=0.2)+
  geom_jitter(aes(y=state_score_cocaine_1d_vs_saline_wt, x=logFC_d01), alpha=0.2, size=0.2)+
  scale_fill_manual(values=c('n'='grey70', 'm'='#4256DE', 'c' = '#FF6699'))+
  theme(legend.position = 'none')

#fig S18C
scores_allcomp_wide_repr %>% group_by(state_score_cocaine_14dR1_vs_saline_wt) %>% summarise(median(logFC_d14))
scores_allcomp_wide_repr %>% 
  ggplot()+
  geom_vline(xintercept=0.0488, color='grey50', linewidth=0.5)+
  geom_violin(aes(y=state_score_cocaine_14dR1_vs_saline_wt, x=logFC_d14, fill=state_score_cocaine_14dR1_vs_saline_wt))+
  geom_boxplot(aes(y=state_score_cocaine_14dR1_vs_saline_wt, x=logFC_d14), outlier.shape=NA, width=0.2)+
  geom_jitter(aes(y=state_score_cocaine_14dR1_vs_saline_wt, x=logFC_d14), alpha=0.2, size=0.2)+
  scale_fill_manual(values=c('n'='grey70', 'm'='#4256DE', 'c' = '#FF6699'))+
  theme(legend.position = 'none')
