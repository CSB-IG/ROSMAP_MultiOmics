library(magrittr)
library(dplyr)

#RNAseq data load
path1<-'ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv'
path2<-'ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv'
genes_expre_p1_p6 <- vroom::vroom(file = path1,delim="\t")
genes_expre_p7_p8 <- vroom::vroom(file = path2,delim="\t")

##Join tables of RNAseq data into one
genes_expre_p1_p6$tracking_id=NULL
genes_expre_p7_p8$tracking_id=NULL
expre_junto<-dplyr::left_join(x=genes_expre_p1_p6,y=genes_expre_p7_p8,
                              by=c('gene_id'='gene_id'))

#gene_id column cleaning. (Quito los números después del punto de los ENSEMBLE IDs)
identificadores<-expre_junto %>% 
  pull(gene_id) 
identificadores<-sapply(strsplit(identificadores,".",fixed=T),function(x) x[1])
expre_junto<-expre_junto %>% add_column(identificadores)

#Save file with expression data joined and clean.
vroom::vroom_write(expre_junto, 'expre_junto.csv', delim = ',')
