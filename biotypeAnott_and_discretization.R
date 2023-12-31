library(magrittr)
library(biomaRt)

#Load expression data
expre_junto <- vroom::vroom(file = expre_junto.csv,delim=",")
identificadores<-expre_junto %>% 
  pull(identificadores)

##Annotate gene_biotype
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", 
                             "percentage_gene_gc_content", "gene_biotype"),
              filters = "ensembl_gene_id", 
              values=identificadores,mart=mart)

##Join biomaRt annotations with gene expression data
expre <-dplyr::left_join(x=expre_junto,y=myannot,
                         by=c('identificadores'='ensembl_gene_id'))

##Protein coding data extraction
expression_mat <- filter(expre, gene_biotype=='protein_coding')
##Values only
expression_mat_values <-expression_mat %>% 
  dplyr::select(-gene_id,-gene_biotype,-identificadores,-percentage_gene_gc_content)

#Data discretization
mat_dis<-infotheo::discretize(t(expression_mat_values))

#Save important data
vroom::vroom_write(expre, 'RNAseq_expre_annot_data.csv', delim = ',')
vroom::vroom_write(expression_mat, 'RNAseq_protcod_expremat.csv', delim = ',')
vroom::vroom_write(mat_dis, 'RNAseq_discretized_mat.csv', delim = ',')