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

#### Size of estimated human interactomes (2025-05-22)
rand_net_size <- read.table('random_network_size.txt',
                            header = FALSE,
                            sep = '\t')
colnames(rand_net_size) <- c('Species', 'Node_num', 'Edge_num')

rand_net_stat <- rand_net_size %>%
  group_by(Species) %>%
  summarise(Mean_node_num = round(mean(Node_num), 3),
            Sd_node_num = round(sd(Node_num), 3),
            Mean_edge_num = round(mean(Edge_num),3))

p_1 <- ggplot(rand_net_stat, aes(Mean_node_num, log10(Mean_edge_num))) +
  geom_point(aes(size = log10(Mean_edge_num)), color = '#399385') +
  geom_text(aes(label = Species), 
            size = 1) +
  theme(legend.position.inside = c(.8, .2)) +
  mytheme

ggsave(p_1,
       filename = 'Fig_1.pdf',
       width = 10,
       height = 6,
       units = c('cm'))

#### The relative number of 1:1 orthologs between human and other six species (2025-05-26)
library(readr)
library(igraph)
library(ggpmisc)
hs_ppi <- read.table('integrated_ppi_Hs.txt',
                      header = FALSE)
hs_ppin <- graph_from_data_frame(hs_ppi, directed = FALSE)

ortho_Hs <- read_delim('../data/orthology/ortho1to1_with_human.txt',
                       delim = '\t',
                       skip_empty_rows = FALSE)
ortho_Hs_stat <- ortho_Hs %>%
  filter(Human_gene %in% V(hs_ppin)$name) %>%
  group_by(Species) %>%
  summarise(Ortho_num = n())

ortho_Hs_stat <- ortho_Hs_stat %>%
  mutate(Rela_ortho_num = Ortho_num/length(V(hs_ppin)$name),
         Diver_time = c(567, 567, 87, 87, 1275, 352))

p_2a <- ggplot(ortho_Hs_stat, aes(Diver_time, Rela_ortho_num)) +
  geom_point(aes(size = Ortho_num), color = '#399385') +
  labs(x = 'Divergence time (Million Year Ago)',
       y = 'Relative 1:1 orthologs number with human\n(1:1 orthologs number/Node number in human interactome)') +
  scale_x_reverse() +
  geom_text(aes(label = Species), size = 1) +
  mytheme + 
  theme(legend.position = c(0.2, 0.7)) +
  stat_smooth(formula = y ~ x, method = "lm", se=TRUE, color = '#EB8588', fill = '#EB8588') +
  stat_fit_glance(aes(label = paste("R-squared = ", signif(..r.squared.., digits = 3),
                                    ", p = ", signif(..p.value.., digits = 4), sep = "")),
                  label.x = 0.1, size = 4)

ggsave(p_2a,
       filename = 'Fig_2A.pdf',
       width = 7,
       height = 6,
       units = c('cm'))

#### The relationship between ICC and sequence similarity (2025-07-02)
icc <- data.frame()
for (f in list.files(path = 'ICC/', pattern = '.txt')) {
  file <- paste0('ICC/', f)
  df <- read.table(file, header = TRUE, sep = '\t')
  icc <- rbind(icc, df)
}
icc <- icc %>%
  rename(ICC = mean)

seq_simi <- read_delim('sequence_similarity/sequence_similarity.txt',
                       delim = '\t',
                       skip_empty_rows = FALSE)

ss_icc <- icc %>%
  select(study, species, study_gene, gene, ICC) %>%
  inner_join(seq_simi, by = c('study', 'species', 'study_gene', 'gene'))

p_s4 <- ggplot(ss_icc, aes(similarity, ICC)) +
  geom_point(color = '#399385', shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = '#EB8588') +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho") +
  xlim(0, 1) +
  ylim(0, 1) +
  mytheme +
  facet_wrap(vars(species), ncol = 3)

ggsave(p_s4,
       filename = "Fig_S4.pdf",
       width = 10.8,
       height = 7.2,
       units = c("cm"))


#### Sequence- and network-based function difference of most conserved and divergent orthologs
### GSEA (2025-07-05)
## ICC
library(clusterProfiler)
top_gsea_icc <- data.frame()
for (sp in unique(icc$species)) {
  sp_icc <- icc %>%
    filter(species == sp) %>%
    mutate(transformed_ICC = 2*normalized_ICC - 1) %>%
    arrange(desc(transformed_ICC))
  
  gene_list <- sp_icc$transformed_ICC
  names(gene_list) <- sp_icc$study_gene
  
  res <- GSEA(
    gene_list,
    TERM2GENE = human_pathways[,1:2],
    pAdjustMethod = "BH",
    nPermSimple = 100000,
    seed = 2914,
    by = "fgsea"
  )
  
  sp2 <- gsub("(\\w)\\.\\s(\\w).*", "\\1\\2", sp)
  res <- res@result
  
  if(nrow(res) > 0) {
    write.table(res[,2:11],
                paste0('ICC/GSEA/GSEA_Reactome_', sp2, ".txt"),
                sep = '\t',
                row.names = FALSE,
                quote = FALSE)
    
    top_conser <- res %>%
      filter(NES > 0) %>%
      arrange(desc(NES)) %>%
      head(10)
    
    top_diverg <- res %>%
      filter(NES < 0) %>%
      arrange(NES) %>%
      head(10)
    
    top_all <- bind_rows(
      "Conserved" = top_conser,
      "Divergent" = top_diverg,
      .id = "Group"
    )
    
    top_all <- top_all %>%
      dplyr::select(Group, Description, setSize, p.adjust, leading_edge) %>%
      mutate(GeneRatio = as.numeric(str_extract(leading_edge, "(?<=tags=)\\d+"))/100) %>%
      arrange(desc(Group), GeneRatio)
    
    top_all$Species <- sp
    top_gsea_icc <- rbind(top_gsea_icc, top_all)
    
    p <- ggplot(top_all, aes(GeneRatio, Description)) +
      geom_point(aes(size = -log10(p.adjust), colour = Group)) +
      scale_y_discrete(limits = top_all$Description) +
      scale_color_manual(values = c('#EB8588', '#1080E7')) +
      #facet_wrap(vars(Group), ncol = 2, scales = "free_x") +
      mytheme
    ggsave(p,
           filename = paste0('ICC/GSEA/GSEA_Reactome_', sp2, ".pdf"),
           width = 20,
           height = 12,
           units = c("cm"))
  }
}

# Function statics (2025-07-08)
top_gsea_icc2 <- human_pathway_cata %>%
  dplyr::select(pathway, category) %>%
  distinct() %>%
  inner_join(top_gsea_icc, by = c("pathway" = "Description")) %>%
  group_by(category, Group, Species) %>%
  summarise(number = n())

p_2e <- ggplot(top_gsea_icc2, aes(category, Species)) +
  geom_point(aes(size = number, colour = Group)) +
  scale_color_manual(values = c('#EB8588', '#1080E7')) +
  scale_y_discrete(limits = rev(c('M. musculus', 'R. norvegicus', 'X. laevis',
                              'D. melanogaster', 'C. elegans', 'S. cerevisiae'))) +
  facet_wrap(vars(Group), ncol = 1) +
  mytheme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(p_2e,
       filename = paste0('Fig_2E.pdf'),
       width = 12,
       height = 9,
       units = c("cm"))

## sequence similarity (2025-07-05)
top_gsea_ss <- data.frame()
for (sp in unique(seq_simi$species)) {
  sp_ss <- seq_simi %>%
    filter(species == sp) %>%
    mutate(normalized_SS = min_max_normalize(similarity)) %>%
    mutate(transformed_SS = 2*normalized_SS - 1) %>%
    arrange(desc(transformed_SS))
  
  gene_list <- sp_ss$transformed_SS
  names(gene_list) <- sp_ss$study_gene
  
  res <- GSEA(
    gene_list,
    TERM2GENE = human_pathways[,1:2],
    pAdjustMethod = "BH",
    nPermSimple = 100000,
    seed = 2914,
    by = "fgsea"
  )
  
  sp2 <- gsub("(\\w)\\.\\s(\\w).*", "\\1\\2", sp)
  res <- res@result
  if(nrow(res) > 0) {
    write.table(res[,2:11],
                paste0('sequence_similarity/GSEA/GSEA_Reactome_', sp2, ".txt"),
                sep = '\t',
                row.names = FALSE,
                quote = FALSE)
    
    top_conser <- res %>%
      filter(NES > 0) %>%
      arrange(desc(NES)) %>%
      head(10)
    
    top_diverg <- res %>%
      filter(NES < 0) %>%
      arrange(NES) %>%
      head(10)
    
    top_all <- bind_rows(
      "Conserved" = top_conser,
      "Divergent" = top_diverg,
      .id = "Group"
    )
    
    top_all <- top_all %>%
      dplyr::select(Group, Description, setSize, p.adjust, leading_edge) %>%
      mutate(GeneRatio = as.numeric(str_extract(leading_edge, "(?<=tags=)\\d+"))/100) %>%
      arrange(desc(Group), GeneRatio)
    
    top_all$Species <- sp
    top_gsea_ss <- rbind(top_gsea_ss, top_all)
    
    p <- ggplot(top_all, aes(GeneRatio, Description)) +
      geom_point(aes(size = -log10(p.adjust), colour = Group)) +
      scale_y_discrete(limits = top_all$Description) +
      scale_color_manual(values = c('#EB8588', '#1080E7')) +
      #facet_wrap(vars(Group), ncol = 2, scales = "free_x") +
      mytheme
    ggsave(p,
           filename = paste0('sequence_similarity/GSEA/GSEA_Reactome_', sp2, ".pdf"),
           width = 20,
           height = 12,
           units = c("cm"))
  }
}

# Function statics (2025-07-08)
top_gsea_ss2 <- human_pathway_cata %>%
  dplyr::select(pathway, category) %>%
  distinct() %>%
  inner_join(top_gsea_ss, by = c("pathway" = "Description")) %>%
  group_by(category, Group, Species) %>%
  summarise(number = n())

p_2d <- ggplot(top_gsea_ss2, aes(category, Species)) +
  geom_point(aes(size = number, colour = Group)) +
  scale_color_manual(values = c('#EB8588', '#1080E7')) +
  scale_y_discrete(limits = rev(c('M. musculus', 'R. norvegicus', 'X. laevis',
                                  'D. melanogaster', 'C. elegans', 'S. cerevisiae'))) +
  facet_wrap(vars(Group), ncol = 1) +
  mytheme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(p_2d,
       filename = paste0('Fig_2D.pdf'),
       width = 12,
       height = 9,
       units = c("cm"))








