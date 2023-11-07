##Grafo Rosmap
library(igraph)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
##Graph inference

##Load MI mat
mat_mi_rosmap<- vroom::vroom(file="matriz_coexpre_allAD_11052023_zero.txt") #cargar matriz de coexpre

##Create network from matrix
#Option 1)Specific threshold 
  #g<-mat_mi_rosmap%>% pivot_longer(cols = -gene,names_to='gene_to',values_to='mi') %>% 
  #  filter(mi >= 0.5) %>% 
  #  graph_from_data_frame(d=.,directed=F)

#Option 2) Top 10k edges
#Obtain top 10k edges:
matrix<-subset(mat_mi_rosmap, select = -gene)
mat_colapse<-unlist(matrix, use.names=FALSE) 
x<-tail(sort(mat_colapse),10000) 
x<-as_tibble(x) #From x extract lowest value to use for threshold
g<-mat_mi_rosmap%>% tidyr::pivot_longer(cols = -gene,names_to='gene_to',values_to='mi') %>% 
  filter(mi >= 0.8250788) %>% 
  graph_from_data_frame(d=.,directed=F)


##Olfactory genes subgraph
olfatory_genes<- vroom::vroom(file="ROSMAP_olfaction_genes.csv", delim = ',')

