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

##Graph subset
##Paso 1: Tidying ensemble ids, adding number to each node
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

plot(sub_g,edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
     
     vertex.color='pink',vertex.size=2,vertex.label=NA)

#Medidas de centralidad
Indegree <- degree(sub_g, mode="in")

vertices_olfato$Indegree <-Indegree##Hay un error aquí porque vertices olfato no tiene los primvecinos

vertices_olfato %>% dplyr::group_by(Indegree) %>% 
  tally()%>% ggplot(aes(x=Indegree,y=n))+
  geom_bar(stat="identity")