#Análisis de enriquecimiento para saber funciones de primeros vecinos de nodos olf.
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(tidyverse)

#Load subgraph nodes

nodos_subg<-vroom::vroom(file="/datos/rosmap/metadata/nodos_subg_olf.csv", delim = ',')
#Identify no olfaction nodes
no_olfaction <- subset(nodos_subg, bio_func == "other")

#Enrichment analysis
go_enrich <- enrichGO(gene = no_olfaction$ENSEMBL,
                      universe = nodos_subg$ENSEMBL,
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

#Analysis visualization
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)