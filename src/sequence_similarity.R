setwd('../intermediate/sequence_similarity/trimmed_aligned_fasta/')

library(tidyverse)
library(ape)

study <- c()
species <- c()
study_gene <- c()
gene <- c()
similarity <- c()
for (f in list.files(pattern = '.fasta')) {
  seq <- read.dna(f, format="fasta")
  simi <- 1-dist.dna(seq, model="TN93")
  simi_matr <- as.matrix(simi)
  similarity <- c(similarity, simi_matr[1,2])
  
  name <- str_remove(f, "\\.fasta$")
  stu <- str_split(name, '_')[[1]][3]
  sp <- str_split(name, '_')[[1]][4]
  stu_g <- str_split(name, '_')[[1]][5]
  g <- str_split(name, '_')[[1]][6]
  study <- c(study, stu)
  species <- c(species, sp)
  study_gene <- c(study_gene, stu_g)
  gene <- c(gene, g)
  print(name)
}

seq_simi <- data.frame(study = study,
                       species = species,
                       study_gene = study_gene,
                       gene = gene,
                       similarity = similarity)

seq_simi <- seq_simi %>%
  mutate(similarity = ifelse(similarity > 0, similarity, 0),
         study = 'H. sapiens',
         species = ifelse(species == 'Ce', 'C. elegans', 
                          ifelse(species == 'Dm', 'D. melanogaster',
                                 ifelse(species == 'Mm', 'M. musculus',
                                        ifelse(species == 'Rn', 'R. norvegicus',
                                               ifelse(species == 'Sc', 'S. cerevisiae', 'X. laevis'))))))

write.table(seq_simi,
            '../result/sequence_similarity/sequence_similarity.txt',
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)






