library(GO.db)
library(magrittr)
library(tidyverse)
library(org.Hs.eg.db)

#Part 1: GO anottation search

##GO database
all_go_terms <- keys(org.Hs.eg.db, keytype = "GO")

##GO id + physiological related process
mis_go <- AnnotationDbi::select(GO.db, keys = all_go_terms, columns = c("GOID", "TERM"), keytype = "GOID") %>% 
  as_tibble() 

##Search fo specific biological process or processes
mis_go_olfatorios <- 
  mis_go %>% 
  filter(grepl(x = TERM, "olfactory|olfaction|hyposmia"))

# get list of ensemble id genes for each GO term
mis_genes_olfatorios <- 
  map(.x = mis_go_olfatorios$GOID, .f = function(i){
    associated_genes <- AnnotationDbi::select(org.Hs.eg.db, 
                                              keys = i, 
                                              keytype = "GO", 
                                              columns = "ENSEMBL"
    ) %>% 
      as_tibble()
  })

##List to table: ENSEMBLE ids + GO + physio process
tabla_olf <- 
  mis_genes_olfatorios %>% 
  bind_rows() %>% 
  left_join(mis_go_olfatorios, by=c("GO"="GOID")) 

#Tally
genes_olf_anotados <-
  tabla_olf %>% 
  group_by(ENSEMBL) %>% tally()


##Part 2:Inner join with data ENSEMBLE ids

##Load ensmble ids from data
nombres<-vroom::vroom(file='ROSMAP_ensembl_id.csv',delim=',')#ENSMBLE IDS from ROSMAP RNAseq data

#Clean ids
identificadores<-pull(nombres,identificadores)
identificadores<-sapply(strsplit(identificadores,".",fixed=T),function(x) x[1])
identificadores<-as.data.frame(identificadores)
colnames(identificadores)[1] <- "ENSEMBL"

##Join known GO olfactory ensemble IDS with project data ensemble ids
ROSMAP_olfaction_genes<-dplyr::inner_join(genes_olf_anotados,identificadores, by='ENSEMBL')

#Save
vroom::vroom_write(ROSMAP_olfaction_genes, 'ROSMAP_olfaction_genes.csv', delim = ',')