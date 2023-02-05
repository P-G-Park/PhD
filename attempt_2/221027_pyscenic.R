pacman::p_load(tidyverse, Seurat, SCENIC, ggsci)
options(stringsAsFactors = FALSE)
# 

load('./raw_data/attempt_2/processed_raw_v3_lr.RData')

future::plan('multisession')


#SeuratDisk::as.loom(Seurat_scenic_old, filename = 'old.loom')
#SeuratDisk::as.loom(Seurat_scenic_young, filename = 'young.loom')

# COL <- pal_jco(alpha = 0.8)(9)
# col_vars <- list(cell_type=c("KRM"=COL[1], 
#                              "CDI"=COL[2], 
#                              "PODO"=COL[3], 
#                              "DCT/CDP"=COL[4], 
#                              "PT"=COL[5],
#                              "SMC/PERI"=COL[6],
#                              "EC"=COL[7],
#                              "DL/tAL"=COL[8],
#                              "TAL"=COL[9]))
# 
# dbDir <- "./raw_data/attempt_2/221011_scenic" # RcisTarget databases location
# myDatasetTitle <- "SCENIC" 
# dbs <- list.files(dbDir)[c(1,3)]
# names(dbs) <- c('10kb','500bp')
# scenicOptions <- initializeScenic(org='hgnc', dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
# 
# scenicOptions@inputDatasetInfo$cellInfo <- scenic_meta_old
# scenicOptions@inputDatasetInfo$colVars <- col_vars

###
pacman::p_load(tidyverse, Seurat, SCENIC, ggsci, SCopeLoomR)
pyscenic_old <- SCopeLoomR::open_loom('./python/pyscenic_old.loom')

regulons_incidMat <- get_regulons(pyscenic_old, column.attr.name = 'Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(pyscenic_old, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(pyscenic_old)
embeddings <- get_embeddings(pyscenic_old)

cellInfo <- data.frame(CellType=Idents(Seurat_EC_cc_old))

regulonActivity_byCellType_old <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType_Scaled_old <- t(scale(t(regulonActivity_byCellType_old), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

o <- reshape2::melt(regulonActivity_byCellType_Scaled_old) %>% 
  `colnames<-`(c("Regulon", "cell_type", "RelativeActivity_old")) %>% 
  filter(cell_type == 'KRM') %>% 
  arrange(desc(RelativeActivity_old))

pyscenic_young <- SCopeLoomR::open_loom('./python/pyscenic_young.loom')

regulons_incidMat <- get_regulons(pyscenic_young, column.attr.name = 'Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(pyscenic_young, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(pyscenic_young)
embeddings <- get_embeddings(pyscenic_young)

cellInfo <- data.frame(CellType=Idents(Seurat_EC_cc_young))

regulonActivity_byCellType_young <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType_Scaled_young <- t(scale(t(regulonActivity_byCellType_young), center = T, scale=T))



ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled_young, name="Regulon activity")

y <- reshape2::melt(regulonActivity_byCellType_Scaled_young) %>% 
  `colnames<-`(c("Regulon", "cell_type", "RelativeActivity_young")) %>% 
  filter(cell_type == 'KRM') %>% 
  arrange(desc(RelativeActivity_young))

oo <- tibble(regulon = rownames(regulonActivity_byCellType_old),
       old_activity = regulonActivity_byCellType_old[,'KRM'])

yy <- tibble(regulon = rownames(regulonActivity_byCellType_young),
       young_activity = regulonActivity_byCellType_young[,'KRM'])


anv <- full_join(oo,yy) %>% mutate(old_activity = if_else(is.na(old_activity), 0, old_activity),
                            young_activity = if_else(is.na(young_activity), 0, young_activity),
                            relative_activity = old_activity - young_activity) %>% 
  arrange(desc(relative_activity))

anv[anv$regulon == 'CEBPD(+)',]
str_sort(regulons$`CEBPD(+)`)
