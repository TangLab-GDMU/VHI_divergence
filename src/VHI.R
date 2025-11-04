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

setwd('../data/')

#########################################################
### The viral module (the large connected component, LCC) 
### LCC size and % targets formed LCC (2025-07-18)

for (s in list.files(path = 'VTGs/')) {
  vtg <- data.frame()
  for (f in list.files(path = paste0('VTGs/', s, '/'))) {
    df <- read.table(paste0('VTGs/', s, '/', f),
                     header = TRUE, sep = '\t')
    vtg <- rbind(vtg, df)
  }
  
  lcc_statistic <- data.frame()
  for (f in list.files(path = paste0('../intermediate/viral_module/', s, '/'), pattern = 'lcc_statistic.txt')) {
    v <- sub("_lcc_statistic.txt", "", f)
    df <- read.table(paste0('../intermediate/viral_module/', s, '/', f),
                     header = TRUE, sep = '\t')
    df$Abbre <- v
    lcc_statistic <- rbind(lcc_statistic, df)
  }
  
  lcc_statistic <- vtg %>%
    group_by(Abbre) %>%
    summarise(Number = n()) %>%
    inner_join(lcc_statistic, by = 'Abbre') %>%
    mutate(Proportion = size/Number*100)
  
  write.table(lcc_statistic,
              paste0('../result/viral_module/', s, "_LCC_statistic.txt"),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
  
  tr <- max(lcc_statistic$Proportion)/max(lcc_statistic$size)
  
  p <- ggplot() +
    geom_bar(data = lcc_statistic, aes(Abbre, size), 
             stat="identity", 
             fill = '#399385', 
             alpha = .75) +
    geom_line(data = lcc_statistic, aes(Abbre, Proportion/tr), 
              group = 1, 
              color = '#EB8588', 
              size = 1.2) +
    scale_y_continuous(name = "LCC size",
                       sec.axis = sec_axis(trans=~.*tr, name="LCC proportion (%)"),
                       expand = c(0.02, 0.02)) +
    mytheme +
    theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1))
  
  ggsave(p,
         filename = paste0('../result/viral_module/', s, "_LCC_statistic.pdf"),
         width = 9,
         height = 6,
         units = c("cm"))
}

## Relationship between % targets formed LCC and network scale (2025-07-18)
setwd('../result/')
net_evo_info <- data.frame(Species = c("Ce", "Dm", "Hs", "Mm", "Rn", "Sc", "Xl"),
                           Edge_num = c(32411, 56260, 698482, 187081, 5766, 151766, 1528),
                           Node_num = c(7502, 8859, 18101, 11297, 2837, 4982, 1131),
                           Diver_time = c(567, 567, 0, 87, 87, 1275, 352))
net_evo_info <- net_evo_info %>%
  mutate(Density = (2*Edge_num)/(Node_num*(Node_num-1)))


all_lcc_statistic <- data.frame()
for (f in list.files(path = 'viral_module/', pattern = '_LCC_statistic.txt')) {
  s <- sub("_LCC_statistic.txt", "", f)
  lcc_statistic <- read.table(paste0('viral_module/', f),
                              header = TRUE,
                              sep = '\t')
  lcc_statistic$Species <- s
  all_lcc_statistic <- rbind(all_lcc_statistic, lcc_statistic)
}

lcc_summary <- Rmisc::summarySE(all_lcc_statistic, measurevar="Proportion", 
                                groupvars = c("Species"))
lcc_summary <- lcc_summary %>%
  inner_join(net_evo_info, by = 'Species')
lcc_summary2 <- lcc_summary %>%
  filter(Species != 'Hs')

p_3a <- ggplot(lcc_summary2, aes(Density, Proportion)) +
  geom_smooth(method = "lm", se = TRUE, color = '#EB8588', fill = '#EB8588') +
  geom_point(color = '#399385') +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho",
                   label.x = min(lcc_summary2$Density), label.y = 80) +
  labs(x = 'Network density', y = 'LCC proportion (%)') +
  mytheme

p_3b <- ggplot(lcc_summary2, aes(Diver_time, Proportion)) +
  geom_smooth(method = "lm", se = TRUE, color = '#EB8588', fill = '#EB8588') +
  geom_point(color = '#399385') +
  scale_x_reverse() +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho",
                   label.y = 80) +
  labs(x = 'Divergence time (Million Year Ago)', y = 'LCC proportion (%)') +
  mytheme

p_prop <- ggpubr::ggarrange(p_3a, p_3b,
                            nrow=1, ncol = 2, labels = c("A", "B"),
                            font.label = list(size = 12))
ggsave(p_prop,
       filename = "Fig_3AB.pdf",
       width = 10,
       height = 4.5,
       units = c("cm"))

## Comparison of % targets formed LCC across species (2025-07-19)
library(igraph)
rand_lcc_statistic <- data.frame()
for (s in c('Ce', 'Dm', 'Mm', 'Rn', 'Sc', 'Xl')) {
  for (vf in list.files('../data/VTGs/Hs/')) {
    v <- sub(".txt", "", vf)
    targets <- read.table(paste0('../data/VTGs/Hs/', vf),
                          header = TRUE,
                          sep = '\t')
    for (i in 0:999) {
      ppi <- read.table(paste0('../intermediate/random_network/', s, '/edge_list_', i, '.txt'),
                        header = FALSE,
                        sep = '\t',
                        quote = '')
      ppin <- graph_from_data_frame(ppi, directed = FALSE)
      nodes <- vertex_attr(ppin)
      nodes <- nodes$name
      
      df <- read.table(paste0('../intermediate/viral_module/', s, '/intermediate/edge_list_', i, '_', v, '_lcc_statistic.txt'),
                       header = TRUE,
                       sep = '\t')
      df$VTG_num <- length(intersect(nodes, targets$Gene))
      df$netID <- i
      df$abbre <- v
      df$species <- s
      
      rand_lcc_statistic <- rbind(rand_lcc_statistic, df)
    }
  }
}

rand_lcc_statistic <- rand_lcc_statistic %>%
  mutate(proportion = size/VTG_num*100)
write.table(rand_lcc_statistic,
            '../intermediate/viral_module/random_network_LCC_statistic.txt',
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

rand_lcc_statistic <- read.table('../intermediate/viral_module/random_network_LCC_statistic.txt',
                                 sep = '\t',
                                 header = TRUE)
rand_lcc_summary <- Rmisc::summarySE(rand_lcc_statistic, measurevar="proportion", 
                                     groupvars = c("species", 'abbre'))

rand_lcc_summary2 <- Rmisc::summarySE(rand_lcc_statistic, measurevar="proportion", 
                                     groupvars = c("species"))

# Density of random network (2025-07-19)
mean_N <- c()
mean_M <- c()
for (s in c('Ce', 'Dm', 'Mm', 'Rn', 'Sc', 'Xl')) {
  Ns <- c()
  Ms <- c()
  for (i in 0:999) {
    ppi <- read.table(paste0('../intermediate/random_network/', s, '/edge_list_', i, '.txt'),
                      header = FALSE,
                      sep = '\t',
                      quote = '')
    ppin <- graph_from_data_frame(ppi, directed = FALSE)
    nodes <- vertex_attr(ppin)
    nodes <- nodes$name
    N <- length(nodes)
    edges <- E(ppin)
    M <- length(edges)
    Ns <- c(Ns, N)
    Ms <- c(Ms, M)
  }
  mean_N <- c(mean_N, mean(Ns))
  mean_M <- c(mean_M, mean(Ms))
}

rand_net_dens <- data.frame(species = c('Ce', 'Dm', 'Mm', 'Rn', 'Sc', 'Xl'),
                            mean_N = mean_N,
                            mean_M = mean_M)
rand_net_dens <- rand_net_dens %>%
  mutate(Density = (2*mean_M)/(mean_N*(mean_N-1)))

rand_lcc_summary <- rand_lcc_summary %>%
  inner_join(rand_net_dens, by = 'species')
write.table(rand_lcc_summary,
            '../intermediate/viral_module/random_network_LCC_summary.txt',
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)
rand_lcc_summary <- read.table('../intermediate/viral_module/random_network_LCC_summary.txt',
                               header = TRUE,
                               sep = '\t')

library(ggpmisc)
p_3c <- ggplot(rand_lcc_summary, aes(Density, proportion)) +
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
  geom_errorbar(aes(ymin=proportion-se, ymax=proportion+se), 
                color = '#399385') +
  mytheme

ggsave(p_3c,
       filename = "Fig_3C.pdf",
       width = 5,
       height = 4.5,
       units = c("cm"))

p_s8 <- ggplot(rand_lcc_summary, aes(Density, proportion)) +
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
               size = 1) +
  geom_point(color = '#399385') +
  geom_errorbar(aes(ymin=proportion-se, ymax=proportion+se), 
                color = '#399385') +
  facet_wrap(vars(abbre), ncol = 5, scales = "free_y") +
  mytheme

ggsave(p_s8,
       filename = "Fig_S8.pdf",
       width = 18,
       height = 21.6,
       units = c("cm"))

# GLM: parameter estimation (2025-07-20)
lcc_summary_test <- rand_lcc_summary %>%
  filter(abbre == 'ZIKV') %>%
  select(Density, proportion) %>%
  rename(x=Density, y = proportion)

glm.sol <- glm(y ~ I(1/x), family = gaussian(link = inverse), data = lcc_summary_test)
summary(glm.sol)

# Comparison of LCC proportion across species under same network scale (2025-07-21) 

lcc_compar <- rand_lcc_summary %>%
  select(species, abbre, proportion, sd, Density) %>%
  rename(Hs_proportion = proportion)

lcc_compar <- all_lcc_statistic %>%
  select(Abbre, Proportion, Species) %>%
  rename(abbre = Abbre, sp_proportion = Proportion, species = Species) %>%
  inner_join(lcc_compar, by = c('abbre', 'species'))  %>%
  mutate(proportion_diff = sp_proportion - Hs_proportion,
         z.score = (sp_proportion - Hs_proportion)/sd,
         diff_dire = ifelse(proportion_diff > 0, 'up', 'down'),
         p.value = 2 * pnorm(abs(z.score), lower.tail = FALSE))

for (s in c('Ce', 'Dm', 'Mm', 'Rn', 'Sc', 'Xl')) {
  lcc_diff <- lcc_compar %>%
    filter(species == s)
  
  p <- ggplot(lcc_diff) +
    geom_segment(aes(x = abbre, xend = abbre, y = Hs_proportion, yend = sp_proportion,
                     colour = diff_dire),
                 arrow = arrow(length = unit(0.2, "cm")),
                 linewidth = 1.5) +
    geom_text(aes(x = abbre, y = sp_proportion, label = round(z.score, 3)), 
              vjust = -0.5, size = 2) +
    scale_color_manual(values = c('#1080E7', '#EB8588')) +
    scale_x_discrete(limits = lcc_diff$abbre[order(lcc_diff$proportion_diff)]) +
    labs(x = '', y = 'LCC proportion (%) diff.') +
    mytheme +
    theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1),
          legend.position = 'top')
  
  ggsave(p,
         filename = paste0('viral_module/', s, "_LCC_proportion_diff.pdf"),
         width = 10,
         height = 7,
         units = c("cm"))
}

lcc_proportion_sp <- lcc_compar %>%
  select(sp_proportion, species, Hs_proportion, proportion_diff) %>%
  inner_join(net_evo_info, by = c('species' = 'Species'))

# The host range of each virus (2025-08-14)
vh <- read_delim('../data/virus_host/virus_host_eid2',
                 delim = '\t',
                 skip_empty_rows = FALSE)
vh <- vh %>%
  filter(`Carrier Rank` == 'species') %>%
  select(Carrier, `Carrier TaxId`, `Carrier Rank`, `Carrier Taxa`, Cargo, `Cargo TaxId`) %>%
  distinct()

virus <- read_delim('../data/virus_host/viruses_eid2',
                    delim = '\t',
                    skip_empty_rows = FALSE)

vh2 <- inner_join(vh, virus, by = c('Cargo' = 'eid2')) %>%
  select(Virus, Abbre, Carrier, `Carrier TaxId`, `Carrier Taxa`) %>%
  add_row(Virus = "Human papillomavirus 8", Abbre = 'HPV-8', 
          Carrier = "homo sapiens", `Carrier TaxId` = 9606,
          `Carrier Taxa` = 'primates') %>%
  distinct()

write.table(vh2,
            '../data/virus_host/virus_host_summary.txt',
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

carrier_num <- vh2 %>%
  group_by(Abbre) %>%
  summarise(Carrier_num = n())
  
# The correlation between LCC proportion and carrier number of virus (2025-08-14)
lcc_carrier <- lcc_compar %>%
  inner_join(carrier_num, by = c('abbre' = 'Abbre')) %>%
  filter(abbre != 'H5N1')
  
p_s10 <- ggplot(lcc_carrier, aes(Carrier_num, proportion_diff)) +
  geom_smooth(method = "lm", se = TRUE, color = '#EB8588', fill = '#EB8588') +
  geom_point(color = '#399385', shape = 1) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho",
                   label.y = -10) +
  labs(x = 'Carrier number', y = 'LCC proportion (%) diff.') +
  facet_wrap(vars(species), ncol = 3, scales = "free_y") +
  mytheme

ggsave(p_s10,
       filename = "Fig_S10.pdf",
       width = 10.8,
       height = 7.2,
       units = c("cm"))

