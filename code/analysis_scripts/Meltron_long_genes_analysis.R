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
