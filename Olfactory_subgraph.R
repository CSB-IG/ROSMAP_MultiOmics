##Grafo Rosmap
library(igraph)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
##Graph inference

##Load MI mat
mat_mi_rosmap<- vroom::vroom(file="matriz_coexpre_allAD_11052023_zero.txt") #cargar matriz de coexpre

##Network from matrix
g<-mat_mi_rosmap%>% tidyr::pivot_longer(cols = -gene,
  names_to='gene_to',values_to='mi') %>% 
  graph_from_data_frame(d=.,directed=F)

#Subset edges
valor_corte<-0.8
borrar_enlaces <-E(g)[mi <valor_corte]
g<- delete_edges(g,edges = borrar_enlaces)

                      ##Olfactory genes subgraph##
#Load olfactory genes from dataset
olfatory_genes<- vroom::vroom(file="ROSMAP_olfaction_genes.csv", delim = ',')

##Graph subset
##1) Tidying ensemble ids, adding number to each node
vertices<-vertex_attr(g) %>% as_tibble() #Nombres nodos
numero_vertices<-row.names(vertices) #Número nodos
vertices$number <-numero_vertices
identificadores<-sapply(strsplit(pull(vertices,name),".",fixed=T),function(x) x[1])
vertices$'ENSEMBL' <- identificadores

##Search olfaction genes inside the network
vertices_olfato<-inner_join(vertices,olfatory_genes,by='ENSEMBL')

##First neighbors of olfaction genes
primeros_vecinos <-adjacent_vertices(g,as.numeric(vertices_olfato$number),mode = 'all')
primeros_vecinos<-unlist(primeros_vecinos) %>% as_tibble()

#number of vertex of interest to do graph subset
vertices_interes_olf<-vertices_olfato %>% pull('number') %>% as.numeric()
vertices_interes<- c(vertices_interes_olf, primeros_vecinos %>% pull(value))

#Subgraph from vertex of interest
sub_g<-subgraph(g,vertices_interes)

#subgraph vertexs with biofun anotation
nodos_subg<-vertex_attr(sub_g) %>% as_tibble() 
nodos_subg<-inner_join(vertices,nodos_subg)
nodos_subg$bio_func <- ifelse(nodos_subg$number %in% vertices_interes_olf, "olfaction", "other")

#Subgraph centrality
node_degree <- degree(sub_g) 
nodos_subg$degree <-node_degree

Betweenness <- betweenness(sub_g)
nodos_subg$Betweenness <-Betweenness

##Subgraph global description
subg_size<-gsize(sub_g)
componentes<-components(sub_g)
cc<-transitivity(sub_g,type = "global")
cc_local<-transitivity(sub_g,type = "local")