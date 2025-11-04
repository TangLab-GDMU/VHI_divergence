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

### Orthologs among 7 species (2025-04-22)
group_out <- read.table("../../data/orthology/gene_groups_out.txt", 
                        header = TRUE, sep = '\t', row.names= 1)

## The evolutionary relationship (2025-04-22)
library(ggtree)
library(treeio)
library(ggimage)

tree <- read.newick("../../data/orthology/species.nwk")

# figure S3A
p_s3a <- ggtree(tree, size = 1.5) +
  xlim(NA, 2500) +
  geom_tiplab(aes(label = paste0("italic('", label, "')")), 
              parse = TRUE, linewidth = 4, offset=10) +
  geom_tiplab(aes(image = paste0("../../data/orthology/species_images/", label, '.png')),
              geom = "image", offset=900,
              size = .1) +
  geom_rootedge(rootedge = 150,size=1.5,color="black")

ggsave(p_s3a,
       filename = '../../result/Fig_S3A.pdf',
       width = 6,
       height = 3,
       units = c('cm'))


## Gene family size (2025-04-22)
group_out2 <- group_out
group_out2[group_out2 >= 3] <- ">= 3"
group_out2 <- group_out2 %>%
  pivot_longer(Hs:Sc, names_to = "Species", values_to = "Group_size") %>%
  filter(Group_size != 0) %>%
  group_by(Species, Group_size) %>%
  dplyr::summarise(Group_number = n())
group_out2$Species <- factor(group_out2$Species, 
                             levels = rev(c("Hs", "Mm", 'Rn', 'Xl', "Dm", 'Ce', "Sc")),
                             labels = rev(c("H. sapiens", "M. musculus", "R. norvegicus", "X. laevis",
                                            "D. melanogaster", "C. elegans", "S. cerevisiae")))
group_out2$Group_size <- factor(group_out2$Group_size, 
                             levels = rev(c(1, 2, ">= 3")))
group_out2 %>% 
  group_by(Species) %>% 
  mutate(prop = Group_number/ sum(Group_number)) %>% 
  View

# figure S3B
p_s3b <- ggplot(group_out2, aes(Species, Group_number, fill = Group_size)) +
  geom_bar(stat = 'identity', position="fill", alpha = 0.9, width = .95) +
  scale_fill_manual(values = c('#AEE1D5', '#53AE9F', '#2C756C')) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  labs(y = 'Frequency', x = '', fill = 'Family size') +
  coord_flip() +
  mytheme +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = 'top')

ggsave(p_s3b,
       filename = '../../result/Fig_S3B.pdf',
       width = 5,
       height = 3)

## One-to-one orthologs (2025-04-22)
ortho1to1 <- group_out
colnames(ortho1to1) <- c("H. sapiens", "M. musculus", "R. norvegicus", "X. laevis",
                         "D. melanogaster", "C. elegans", "S. cerevisiae")
ortho1to1$Group <- rownames(ortho1to1)
ortho1to1 <- ortho1to1 %>%
  rowwise() %>%
  mutate(Max = max(c_across(`H. sapiens`:`S. cerevisiae`))) %>%
  filter(Max == 1) %>%
  dplyr::select(Group, `H. sapiens`:`S. cerevisiae`)

ortho1to1 <- as.data.frame(ortho1to1)

ortho1to1_lt <- as.list(ortho1to1)

for (i in 2:8) {
  index <- which(ortho1to1_lt[names(ortho1to1_lt)[i]][[1]] != 0)
  ortho1to1_lt[[names(ortho1to1_lt)[i]]] <- ortho1to1_lt$Group[index]
}

ortho1to1_lt <- ortho1to1_lt[-1]
library(ComplexHeatmap)
ortho1to1_m <- make_comb_mat(ortho1to1_lt)

# figure S3C
p_s3c <- UpSet(ortho1to1_m)

pdf('../../result/Fig_S3C.pdf',
    width = 16,
    height = 7)
draw(p_s3c)
dev.off()

## The common one-to-one orthologs with human (2025-04-22)
ortho_hs <- group_out
colnames(ortho_hs) <- c("H. sapiens", "M. musculus", "R. norvegicus", "X. laevis",
                        "D. melanogaster", "C. elegans", "S. cerevisiae")
ortho_hs$OG <- rownames(ortho_hs)

ortho_hs <- ortho_hs %>%
  rowwise() %>%
  mutate(Max = max(c_across(`H. sapiens`:`S. cerevisiae`))) %>%
  filter(Max == 1 & `H. sapiens` == 1) %>%
  dplyr::select(OG, `H. sapiens`:`S. cerevisiae`)

ortho_hs2 <- ortho_hs %>%
  pivot_longer(`M. musculus`:`S. cerevisiae`, names_to = "Species", values_to = "Gene_num") %>%
  filter(Gene_num != 0)

# figure S3D
p_s3d <- ggplot(ortho_hs2, aes(Species)) +
  geom_bar(stat = "count") +
  scale_x_discrete(limits = rev(c("M. musculus", "R. norvegicus", "X. laevis",
                                  "D. melanogaster", "C. elegans", "S. cerevisiae"))) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  labs(x = "", y = "Number of one-to-one orthologs\nwith H. sapiens") +
  coord_flip() +
  mytheme +
  theme(plot.margin = unit(c(.5, .5, .5, .5), "cm"))

ggsave(p_s3d,
       filename = '../../result/Fig_S3D.pdf',
       width = 3,
       height = 2)

