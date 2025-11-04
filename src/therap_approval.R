library(tidyverse)
library(RColorBrewer)
mytheme <- theme(panel.grid.major = element_line(colour=brewer.pal(9,"Pastel1")[9],
                                                 linetype = "solid"),
                 panel.background = element_rect(fill='transparent', color="#000000"),
                 panel.border=element_rect(fill='transparent', color='black'),
                 axis.text = element_text(size = 7),
                 axis.title = element_text(size = 8),
                 legend.text = element_text(size = 7),
                 legend.title = element_text(size = 8),
                 legend.background = element_blank())

setwd('../result/')

### Data loading
## The probabilities of success (PoSs) of clinical trials for vaccines and drugs (2025-09-25)
library(readxl)
pos_indu_spon_vacc <- read_excel('../data/PoS_vaccine&drug/PoS_Lo2020.xlsx',
                                 sheet = 'industry-sponsored_vaccine',
                                 skip = 1)

pos_non_indu_spon_vacc <- read_excel('../data/PoS_vaccine&drug/PoS_Lo2020.xlsx',
                                 sheet = 'non-industry-sponsored_vaccine',
                                 skip = 1)

pos_indu_spon_drug <- read_excel('../data/PoS_vaccine&drug/PoS_Lo2020.xlsx',
                                 sheet = 'industry-sponsored_drug',
                                 skip = 1)

pos_non_indu_spon_drug <- read_excel('../data/PoS_vaccine&drug/PoS_Lo2020.xlsx',
                                 sheet = 'non-industry-sponsored_drug',
                                 skip = 1)

pos_list <- list(
  indu_spon_vacc = pos_indu_spon_vacc,
  non_indu_spon_vacc = pos_non_indu_spon_vacc,
  indu_spon_drug = pos_indu_spon_drug,
  non_indu_spon_drug = pos_non_indu_spon_drug
)

## The correspondence between viruses and diseases (2025-09-25)
v_d <- read_excel('../data/PoS_vaccine&drug/viruses_disease.xlsx')

pos1A <- v_d %>%
  left_join(pos_indu_spon_vacc, by = 'Disease') %>%
  select(Virus, Abbre, Disease, PoS1A) %>%
  rename(`industry-sponsored vaccine` = PoS1A) %>%
  left_join(pos_non_indu_spon_vacc, by = 'Disease') %>%
  select(Virus, Abbre, Disease, `industry-sponsored vaccine`, PoS1A) %>%
  rename(`non-industry-sponsored vaccine` = PoS1A) %>%
  left_join(pos_indu_spon_drug, by = 'Disease') %>%
  select(Virus, Abbre, Disease, `industry-sponsored vaccine`,`non-industry-sponsored vaccine`,
         PoS1A) %>%
  rename(`industry-sponsored drug` = PoS1A) %>%
  left_join(pos_non_indu_spon_drug, by = 'Disease') %>%
  select(Virus, Abbre, Disease, `industry-sponsored vaccine`,`non-industry-sponsored vaccine`,
         `industry-sponsored drug`, PoS1A) %>%
  rename(`non-industry-sponsored drug` = PoS1A)

library(openxlsx)
write.xlsx(pos1A, 'therapeutics_approval/PoS1A.xlsx')

### The correlation between relative IC and PoS1A
## Relative IC between H. sapiens and M. musculus (2025-09-26)
ic <- read_delim('IC/IC_ortho1to1_with_human.txt',
                 delim = '\t',
                 skip_empty_rows = FALSE)

relative_ic <- ic %>%
  select(study, species, study_gene, gene, IC, mean, sd, Z) %>%
  mutate(relative_IC = log2(mean/IC)) %>%
  filter(!is.na(relative_IC))

vtg <- read_delim('../data/VTGs/virally-targeted_genes29.txt',
                  delim = '\t',
                  skip_empty_rows = FALSE)

relat_ic_vtg <- relative_ic %>%
  inner_join(vtg, by = c("study_gene" = "Gene"), relationship = "many-to-many")

relat_ic_vtg_Mm <- relat_ic_vtg %>%
  filter(species == 'M. musculus') %>%
  select(Abbre, relative_IC)

ic_vtg_Mm_summary <- Rmisc::summarySE(relat_ic_vtg_Mm, measurevar="relative_IC", 
                                     groupvars = c("Abbre"))

## The correlation
for (i in names(pos_list)) {
  pos <- pos_list[[i]] %>%
    select(Disease, PoS1A) %>%
    inner_join(v_d, by = 'Disease')
  
  pos <- pos %>%
    select(-Virus) %>%
    inner_join(ic_vtg_Mm_summary, by = 'Abbre')
  
  p <- ggplot(pos, aes(relative_IC, PoS1A)) +
    geom_smooth(method = "lm", se = TRUE, color = '#EB8588', fill = '#EB8588') +
    geom_point(color = '#399385', shape = 1) +
    ggpubr::stat_cor(method = "spearman", 
                     cor.coef.name = "rho") +
    labs(x = 'Relative IC', y = 'The probabilities of success (PoSs)') +
    mytheme
  
  ggsave(p,
         filename = paste0('therapeutics_approval/IC_PoS_', i, '.pdf'),
         width = 3.6,
         height = 3.6,
         units = c("cm"))
}

### The correlation between LCC proportion and PoS1A
## LCC proportion between H. sapiens and M. musculus (2025-09-26)
rand_lcc_summary <- read.table('../intermediate/viral_module/random_network_LCC_summary.txt',
                               header = TRUE,
                               sep = '\t')
lcc_compar <- rand_lcc_summary %>%
  select(species, abbre, proportion, sd, Density) %>%
  rename(Hs_proportion = proportion)

all_lcc_statistic <- data.frame()
for (f in list.files(path = 'viral_module/', pattern = 'LCC_statistic.txt')) {
  s <- sub("_LCC_statistic.txt", "", f)
  lcc_statistic <- read.table(paste0('viral_module/', s, '_LCC_statistic.txt'),
                              header = TRUE,
                              sep = '\t')
  lcc_statistic$Species <- s
  all_lcc_statistic <- rbind(all_lcc_statistic, lcc_statistic)
}

lcc_compar <- all_lcc_statistic %>%
  select(Abbre, Proportion, Species) %>%
  rename(abbre = Abbre, sp_proportion = Proportion, species = Species) %>%
  inner_join(lcc_compar, by = c('abbre', 'species'))  %>%
  mutate(proportion_diff = sp_proportion - Hs_proportion,
         z.score = (sp_proportion - Hs_proportion)/sd,
         diff_dire = ifelse(proportion_diff > 0, 'up', 'down'),
         p.value = 2 * pnorm(abs(z.score), lower.tail = FALSE))

lcc_compar_Mm <- lcc_compar %>%
  filter(species == 'Mm') %>%
  select(abbre, proportion_diff) %>%
  rename(Abbre = abbre)

## The correlation
for (i in names(pos_list)) {
  pos <- pos_list[[i]] %>%
    select(Disease, PoS1A) %>%
    inner_join(v_d, by = 'Disease')
  
  pos <- pos %>%
    select(-Virus) %>%
    inner_join(lcc_compar_Mm, by = 'Abbre')
  
  p <- ggplot(pos, aes(proportion_diff, PoS1A)) +
    geom_smooth(method = "lm", se = TRUE, color = '#EB8588', fill = '#EB8588') +
    geom_point(color = '#399385', shape = 1) +
    ggpubr::stat_cor(method = "spearman", 
                     cor.coef.name = "rho") +
    labs(x = 'LCC proportion (%) diff.', y = 'The probabilities of success (PoSs)') +
    mytheme
  
  ggsave(p,
         filename = paste0('therapeutics_approval/LCC_proportion_PoS_', i, '.pdf'),
         width = 3.6,
         height = 3.6,
         units = c("cm"))
}

### The correlation between sequence similarity and PoS1A
## Sequence similarity between H. sapiens and M. musculus (2025-09-26)
seq_simi <- read_delim('sequence_similarity/sequence_similarity.txt',
                       delim = '\t',
                       skip_empty_rows = FALSE)
seq_simi_Mm <- seq_simi %>%
  filter(species == 'M. musculus') %>%
  select(study_gene, similarity) %>%
  inner_join(vtg, by = c("study_gene" = "Gene"), relationship = "many-to-many")

ss_vtg_Mm_summary <- Rmisc::summarySE(seq_simi_Mm, measurevar="similarity", 
                                      groupvars = c("Abbre"))

## The correlation
for (i in names(pos_list)) {
  pos <- pos_list[[i]] %>%
    select(Disease, PoS1A) %>%
    inner_join(v_d, by = 'Disease')
  
  pos <- pos %>%
    select(-Virus) %>%
    inner_join(ss_vtg_Mm_summary, by = 'Abbre')
  
  p <- ggplot(pos, aes(similarity, PoS1A)) +
    geom_smooth(method = "lm", se = TRUE, color = '#EB8588', fill = '#EB8588') +
    geom_point(color = '#399385', shape = 1) +
    ggpubr::stat_cor(method = "spearman", 
                     cor.coef.name = "rho") +
    labs(x = 'Sequence similarity', y = 'The probabilities of success (PoSs)') +
    mytheme
  
  ggsave(p,
         filename = paste0('therapeutics_approval/SS_PoS_', i, '.pdf'),
         width = 3.6,
         height = 3.6,
         units = c("cm"))
}






