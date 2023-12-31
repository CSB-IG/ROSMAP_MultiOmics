#Script to convert RDS file to matrix
#Script by Aidee Lashmi

#libraries ---

library(tidyverse)
library(Matrix)
mi_list <- read_rds("ROSMAP_RNAseq_MutualInfo_allMCI_matrix.rds")

#Unneest RDS ---

coexpre_mat_unnested <-as.data.frame(do.call(cbind, mi_list))

coexpre_mat_unnested %>% as_tibble()

mm <- coexpre_mat_unnested %>% as.matrix() 

mm2 <- mm

mm2[upper.tri(mm2)]  <- mm2[lower.tri(mm2)]

#mm2   %>% as_tibble() %>% pull(1) %>% unlist() %>% as.numeric() %>% length() #vroom::vroom_write(x = "/datos/rosmap/matriz_coexpre_20231011.txt")

mm3 <- mm2   %>% as_tibble() %>% mutate(V1 = colnames(.)) %>% pivot_longer(cols = -V1, names_to = "V2", values_to = "value") 
mm3 <- mm3 %>% mutate(value = mm3$value %>% sapply(FUN = function(i){ifelse(is.null(i), 0, i)}))

#Write mew matrix

mm3 %>% vroom::vroom_write("matriz_coexpre_allMCI_11052023.txt") #change name file every time

mm3 %>% 
  rename(gene = V1) %>% 
  pivot_wider(id_cols = gene, names_from = V2, values_from = value) %>% 
  vroom::vroom_write("matriz_coexpre_allMCI_11052023.txt")
