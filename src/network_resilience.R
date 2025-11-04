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

net_evo_info <- data.frame(Species = c("Ce", "Dm", "Hs", "Mm", "Rn", "Sc", "Xl"),
                           Edge_num = c(32411, 56260, 698482, 187081, 5766, 151766, 1528),
                           Node_num = c(7502, 8859, 18101, 11297, 2837, 4982, 1131),
                           Diver_time = c(567, 567, 0, 87, 87, 1275, 352))
net_evo_info <- net_evo_info %>%
  mutate(Density = (2*Edge_num)/(Node_num*(Node_num-1)))

#### The ICs between human and other six species
ic <- read_delim('IC/IC_ortho1to1_with_human.txt',
                 delim = '\t',
                 skip_empty_rows = FALSE)

### Relative ICs (2025-08-18)
relative_ic <- ic %>%
  select(study, species, study_gene, gene, IC, mean, sd, Z) %>%
  mutate(relative_IC = log2(mean/IC)) %>%
  filter(!is.na(relative_IC))

## The relationship between IC and sequence similarity (2025-08-18)
seq_simi <- read_delim('sequence_similarity/sequence_similarity.txt',
                       delim = '\t',
                       skip_empty_rows = FALSE)

ss_ic <- relative_ic %>%
  select(study, species, study_gene, gene, relative_IC) %>%
  inner_join(seq_simi, by = c('study', 'species', 'study_gene', 'gene'))

# The relationship between IC and sequence similarity by VTP categories (2025-08-26)
vtp_cate <- read_delim('../data/VTGs/viralNum_target_a_gene.txt',
                       delim = '\t',
                       skip_empty_rows = FALSE)

ss_ic_vtp <- ss_ic %>%
  left_join(vtp_cate, by = c('study_gene' = 'Gene')) %>%
  mutate(Virally_targeted = replace_na(Virally_targeted, 'no-VTP'))

color_values <- c('#AEE1D5', '#399385', '#224D48')
names(color_values) <- c('no-VTP', 'Specific', 'Pan')

p_s12 <- ggplot(ss_ic_vtp, aes(similarity, relative_IC)) +
  geom_point(aes(color = Virally_targeted), shape = 1) +
  geom_smooth(aes(color = Virally_targeted),
              method = "lm", 
              se = TRUE) +
  ggpubr::stat_cor(aes(color = Virally_targeted),
                   method = "spearman", 
                   cor.coef.name = "rho") +
  scale_color_manual(values = color_values)  +
  labs(x = 'Sequence similarity', y = 'Relative IC') +
  mytheme +
  facet_wrap(vars(species), ncol = 3, scales = 'free')

ggsave(p_s12,
       filename = "Fig_S12.pdf",
       width = 12.2,
       height = 7.2,
       units = c("cm"))

## The relationship between IC and ICC (2025-07-18)
icc <- data.frame()
for (f in list.files(path = 'ICC/', pattern = '.txt')) {
  file <- paste0('ICC/', f)
  df <- read.table(file, header = TRUE, sep = '\t')
  icc <- rbind(icc, df)
}
icc <- icc %>%
  rename(ICC = mean)

icc_ic <- relative_ic %>%
  select(study, species, study_gene, gene, relative_IC) %>%
  inner_join(icc, by = c('study', 'species', 'study_gene', 'gene'))

# The relationship between IC and ICC by VTP categories (2025-08-26)
icc_ic_vtp <- icc_ic %>%
  left_join(vtp_cate, by = c('study_gene' = 'Gene')) %>%
  mutate(Virally_targeted = replace_na(Virally_targeted, 'no-VTP'))


p_s13 <- ggplot(icc_ic_vtp, aes(ICC, relative_IC)) +
  geom_point(aes(color = Virally_targeted), shape = 1) +
  geom_smooth(aes(color = Virally_targeted),
              method = "lm", 
              se = TRUE) +
  ggpubr::stat_cor(aes(color = Virally_targeted),
                   method = "spearman", 
                   cor.coef.name = "rho") +
  scale_color_manual(values = color_values)  +
  labs(x = 'ICC', y = 'Relative IC') +
  mytheme +
  facet_wrap(vars(species), ncol = 3, scales = 'free')

ggsave(p_s13,
       filename = "Fig_S13.pdf",
       width = 12.2,
       height = 7.2,
       units = c("cm"))

### Viral perturbation and local network struction (2025-08-09) 
vtg <- read_delim('../data/VTGs/virally-targeted_genes29.txt',
                  delim = '\t',
                  skip_empty_rows = FALSE)

relat_ic_vtg <- relative_ic %>%
  inner_join(vtg, by = c("study_gene" = "Gene"), relationship = "many-to-many")

for (s in unique(relat_ic_vtg$species)) {
  relat_ic_vtg_s <- relat_ic_vtg %>%
    filter(species == s) %>%
    select(Abbre, relative_IC)
  
  ic_vtg_s_summary <- Rmisc::summarySE(relat_ic_vtg_s, measurevar="relative_IC", 
                                          groupvars = c("Abbre"))
  ic_vtg_s_summary <- ic_vtg_s_summary %>%
    mutate(fragmentation = ifelse(relative_IC < 0, 'decreased', 'increased'))
  
  relat_ic_s <- relative_ic %>%
    filter(species == s)
  
  Mean <- c()
  Sd <- c()
  for (v in ic_vtg_s_summary$Abbre) {
    rand_mean_ic <- c()
    for (i in 1:1000) {
      set.seed(i)
      ic_samp <- sample(relat_ic_s$relative_IC, 
                        size = ic_vtg_s_summary[ic_vtg_s_summary$Abbre==v, 2], 
                        replace = TRUE)
      rand_mean_ic <- c(rand_mean_ic, mean(ic_samp))
    }
    
    Mean <- c(Mean, mean(rand_mean_ic))
    Sd <- c(Sd, sd(rand_mean_ic))
  }
  
  ic_vtg_s_summary$Mean <- Mean
  ic_vtg_s_summary$Sd <- Sd
  
  ic_vtg_s_summary <- ic_vtg_s_summary %>%
    mutate(Z = (relative_IC - Mean)/Sd,
           P = 2*pnorm(-abs(Z)),
           signif = ifelse(P < 0.001, '***',
                           ifelse(P < 0.01, '**',
                                  ifelse(P < 0.05, '*', 'ns'))))
  
  p <- ggplot(ic_vtg_s_summary) +
    geom_point(aes(Abbre, relative_IC, colour = fragmentation)) +
    geom_errorbar(aes(Abbre, relative_IC, ymin=relative_IC-se, ymax=relative_IC+se, colour = fragmentation), 
                  width=.4) +
    scale_color_manual(values = c('#1080E7', '#EB8588')) +
    geom_text(aes(x = Abbre, y = Mean, label = signif)) +
    scale_x_discrete(limits = ic_vtg_s_summary$Abbre[order(ic_vtg_s_summary$relative_IC)]) +
    mytheme +
    theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1),
          legend.position = 'top')
  
  sp <- gsub("(\\w)\\.\\s(\\w).*", "\\1\\2", s)
  ggsave(p,
         filename = paste0('IC/', sp, "_relative_IC.pdf"),
         width = 10,
         height = 7,
         units = c("cm"))
}


# The correlation between relative IC and host range (2025-08-14)
vh2 <- read_delim('../data/virus_host/virus_host_summary.txt',
                  delim = '\t',
                  skip_empty_rows = FALSE)
carrier_num <- vh2 %>%
  group_by(Abbre) %>%
  summarise(Carrier_num = n())

relat_ic_carrier <- relat_ic_vtg %>%
  select(species, relative_IC, Abbre) %>%
  inner_join(carrier_num, by = 'Abbre')

relat_ic_carrier_summary <- Rmisc::summarySE(relat_ic_carrier, measurevar="relative_IC",
                                             groupvars = c("species", "Abbre", "Carrier_num"))

relat_ic_carrier_summary <- net_evo_info %>%
  filter(Species != 'Hs') %>%
  mutate(Species = unique(relat_ic_carrier_summary$species)) %>%
  inner_join(relat_ic_carrier_summary, by = c('Species' = 'species'))

# The mean rank according to relative ICs (2025-08-18)
ranked_ric_carrier_num <- relat_ic_carrier_summary %>%
  group_by(Species) %>%
  mutate(Rank = dense_rank(relative_IC)) %>%
  arrange(Species, relative_IC)

ranked_ric_summary <- Rmisc::summarySE(ranked_ric_carrier_num, measurevar="Rank", 
                                       groupvars = c("Abbre"))

p_4b <- ggplot(ranked_ric_carrier_num, aes(Abbre, Rank)) +
  geom_boxplot() +
  geom_point(aes(colour = Species)) +
  scale_color_manual(values = c('#3F4C8C', '#6A79B0','#67C1EC', '#EB8588', '#DE6736', '#925E48'),
                     limits = c("M. musculus", "R. norvegicus", "X. laevis",
                                "C. elegans", "D. melanogaster", "S. cerevisiae")) +
  scale_x_discrete(limits = ranked_ric_summary$Abbre[order(ranked_ric_summary$Rank)]) +
  mytheme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(p_4b,
       filename = "Fig_4B.pdf",
       width = 16,
       height = 6,
       units = c("cm"))

### The correlation between relative IC and LCC proportion (2025-08-15)
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

lcc_rIC <- lcc_compar %>%
  mutate(species = relat_ic_carrier_summary$Species) %>%
  select(abbre, species, proportion_diff, diff_dire) %>%
  rename(Species = species, Abbre = abbre) %>%
  inner_join(relat_ic_carrier_summary, by = c('Species', 'Abbre'))

p_s14 <- ggplot(lcc_rIC, aes(proportion_diff, relative_IC)) +
  geom_smooth(method = "lm", se = TRUE, color = '#EB8588', fill = '#EB8588') +
  geom_point(color = '#399385', shape = 1) +
  #geom_text(aes(label = Abbre), size = 2) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho") +
  labs(x = 'LCC proportion (%) diff.', y = 'Relative IC') +
  facet_wrap(vars(Species), ncol = 3, scales = "free") +
  mytheme

ggsave(p_s14,
       filename = "Fig_S14.pdf",
       width = 10.8,
       height = 7.2,
       units = c("cm"))

### The correlation between relative IC and network density (2025-08-28)
ic_vtg_summary <- Rmisc::summarySE(relat_ic_vtg, measurevar="IC",
                                   groupvars = c("species", "Abbre"))

rand_net <- read_delim('../intermediate/viral_module/random_network_LCC_summary.txt',
                       delim = '\t',
                       skip_empty_rows = FALSE)

rand_net <- rand_net %>%
  select(species, Density) %>%
  distinct() %>%
  mutate(species = case_match(
    species,
    "Ce" ~ "C. elegans",
    "Dm" ~ "D. melanogaster", 
    "Mm" ~ "M. musculus",
    "Rn" ~ "R. norvegicus",
    "Sc" ~ "S. cerevisiae",
    "Xl" ~ "X. laevis",
    .default = species
  ))

ic_vtg_summary <- ic_vtg_summary %>%
  inner_join(rand_net, by = 'species')


library(ggpmisc)
p_4f <- ggplot(ic_vtg_summary, aes(Density, IC)) +
  stat_smooth(formula = y ~ I(1/x), 
              method = "glm", 
              method.args = list(family = gaussian(link = "inverse")), 
              se=TRUE, color = '#EB8588', fill = '#EB8588') +
  stat_poly_eq(formula = y ~ I(1/x), 
               method = "glm",
               method.args = list(family = gaussian(link = "inverse")), 
               aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")),
               parse = TRUE,
               rr.digits = 3, 
               size = 2) +
  geom_point(color = '#399385') +
  mytheme

ggsave(p_4f,
       filename = "Fig_4F.pdf",
       width = 5,
       height = 4.5,
       units = c("cm"))

p_s15 <- ggplot(ic_vtg_summary, aes(Density, IC)) +
  stat_smooth(formula = y ~ I(1/x), 
              method = "glm", 
              method.args = list(family = gaussian(link = "inverse")), 
              se=TRUE, color = '#EB8588', fill = '#EB8588') +
  stat_poly_eq(formula = y ~ I(1/x), 
               method = "glm",
               method.args = list(family = gaussian(link = "inverse")), 
               aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")),
               parse = TRUE,
               rr.digits = 3, 
               size = 2) +
  geom_point(color = '#399385') +
  facet_wrap(vars(Abbre), ncol = 5, scales = "free_y") +
  mytheme

ggsave(p_s15,
       filename = "Fig_S15.pdf",
       width = 18,
       height = 21.6,
       units = c("cm"))

# GLM: parameter estimation (2025-08-29)
ic_summary_test <- ic_vtg_summary %>%
  filter(Abbre == 'ZIKV') %>%
  select(Density, IC) %>%
  rename(x=Density, y = IC)

glm.sol <- glm(y ~ I(1/x), family = gaussian(link = inverse), data = ic_summary_test)
summary(glm.sol)

