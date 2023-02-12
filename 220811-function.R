##################################################################################################

`%notin%` <- Negate(`%in%`)

##################################################################################################

RawMat_To_SoupXChannel_Using_DropletUtils <- function(raw_matrix, FDR_cut = 0.05){
  require(DropletUtils)
  require(SoupX)
  set.seed(100)
  tod <- raw_matrix
  message("Running emptyDrops")
  e.out <- DropletUtils::emptyDrops(raw_matrix)
  is.cell <- ((e.out$FDR <= FDR_cut) & !is.na(e.out$FDR))
  SoupData <- SoupChannel(tod = tod, toc = tod[,is.cell])
  return(SoupData)
}

##################################################################################################

SoupXChannel_To_Adjust_Matrix_Removing_Amb_mRNA_Using_SoupX <- function(SoupData){
  require(scran)
  require(scMerge)
  require(tidyverse)
  require(SoupX)
  message("Running SoupX")
  SoupData <- setClusters(SoupData, quickCluster(SoupData[['toc']]))
  out <- autoEstCont(SoupData) %>% adjustCounts(roundToInt = TRUE)
  return(out)
}

##################################################################################################

From_Adjust_Matrix_to_Seurat_Normalization <- function(x, MinCell = 3, MinFeature = 200){
  require(Seurat)
  require(scran)
  require(SeuratWrappers)
  require(tidyverse)
  CreateSeuratObject(counts = x, min.cells = MinCell, project = deparse(substitute(x)), min.features = MinFeature) %>% 
    AddMetaData(., metadata = PercentageFeatureSet(., pattern = "^(?i)MT-"), col.name = 'percent.mt') %>% 
    subset(subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15) %>% 
    as.SingleCellExperiment() %>% 
    computeSumFactors(., clusters = quickCluster(.)) %>% 
    logNormCounts() %>% 
    as.Seurat()
}
####  Define QC parameters  #######
min_cells_per_features <- 3 
MinFeature <- 200
MaxFeature <- 2500
MaxMTpercent <- 30
num_HVG <- 2000

##################################################################################################
From_Gz_to_sparseMatrix <- function(pti.data){
  require(tidyverse)
  pti.matrix <- pti.data %>% select(-1) %>% as.matrix %>% 
    as(Class = 'sparseMatrix')
  rownames(pti.matrix) <- pti.data[[1]]
  return(pti.matrix)
}
##################################################################################################
# 220913. clustering
clustering <- function(seurat_obj, resol = 0.5, reduction_method = 'harmony'){
  require(Seurat)
  seurat_obj <- seurat_obj %>% 
    RunUMAP(reduction = reduction_method, dims = 1:30) %>% 
    FindNeighbors(reduction = reduction_method,  dims = 1:30) %>% 
    FindClusters(resolution = resol)
  return(seurat_obj)
}

##################################################################################################
# 220903. define_markers

markers_JCI <- rev(c('CUBN','NRP2','SLC34A1','SLC22A8','SLC5A2',
                     'ALDOB','CFH','CLDN1','SLC12A1','SLC12A3','SLC12A2',
                     'SLC8A1','AQP2','AQP6','KIT','UMOD','ATP6V0D2',
                     'SLC26A4','NPHS1','NPHS2','WT1','PECAM1', 'PLVAP',
                     'FLT1','PDGFRB','ITGA8','PTPRC'))

markers_immune_JCI <- rev(c('PTPRC','CD14','CD68', 'C1QC','LYZ','CD163','FCGR1A','FCGR3A',
                            'S100A8','FCER1A','CD3E','TRAC','CD4','CD8B','CCR7','SELL','FOXP3','GNLY',
                            'TRDC','CD79A','IGHD'))

markers_park <- rev(c('ALDOB','CRYAA', 'ACAA2','PLIN2', 'SLC3A1', 
                         'CUBN', 'SPP1', 'TACSTD2', 'CLDN4',
                         'SLC12A1','UMOD','AQP2','SLC8A1',
                         'CLCNKB','SPINK1', 'ATP6V0D2', 
                         'ACTA2','TAGLN', 'PDGFRB',
                         'FLT1','PECAM1','EMCN',
                         'NPHS1','NPHS2','WT1', 'PTPRC'))

markers_immune_park <- rev(c('PTPRC','CD68', 'CD14','LYZ', 'C1QC', 
                         'FCER1A', 'S100A8', 'CD3E', 'CD40LG',
                         'TRAC', 'CD8B','GNLY','TRDC','CD79A','SELL'))

##################################################################################################

# extract DEGs from seurat object

DEGs <- list()
DEG <- function(Seurat_Obj){
  require(writexl)
  for (i in 0:max(as.integer(levels(Seurat_Obj@active.ident)))){
  Markers_i <- FindMarkers(Seurat_Obj, ident.1 = as.character(i),
                           only.pos = TRUE)
  Markers_i$cluster <- i
  Markers_i$Gene_Name <- row.names(Markers_i)
  Markers_i <- bind_cols(Markers_i[,6:7], Markers_i[,1:5])
  DEGs[[str_c('Cluster ',i)]] <- Markers_i
  }
  write_xlsx(DEGs, path = str_c('./',deparse(substitute(Seurat_Obj)),'.xlsx'))
  return(DEGs)
}

DEG_ident <- function(Seurat_Obj){
  require(writexl)
  for (i in levels(Idents(Seurat_Obj))){
    Markers_i <- FindMarkers(Seurat_Obj, ident.1 = i, only.pos = TRUE)
    Markers_i$cluster <- i
    Markers_i$Gene_Name <- row.names(Markers_i)
    Markers_i <- bind_cols(Markers_i[,6:7], Markers_i[,1:5])
    DEGs[[i]] <- Markers_i
  }
  write_xlsx(DEGs, path = str_c('./',deparse(substitute(Seurat_Obj)),'.xlsx'))
  return(DEGs)
}

##################################################################################################
# Seurat plotting function

quick <- function(x){DimPlot(x, reduction = "umap", label = TRUE)+
    theme(plot.title = element_blank(), aspect.ratio = 1
          ) & NoLegend()}
quickdot <- function(x, feat = markers_park_221025){DotPlot(x, features = feat) +
    theme(axis.text.x = element_text(hjust = 1,angle = 30))+
    coord_flip()}
quickplot <- function(x){CombinePlots(list(quick(x), quickdot(x, feat = c(markers_immune_JCI, 'HBB','HBA1'))))}


###########################################################################

# 220928 cellphone DB
#' cell-cell interaction analysis using method CellPhoneDB
#'
#' This function allows you to use an CellPhoneDB on Seurat Object
#'
#' @param Seurat_obj Seurat Object
#'
#' @return output tables for cell-cell interction inference
#'
#' @keywords Seurat, single cell sequencing, RNA-seq, cell-cell interactions
#'
#' @examples
#'
#'
#'
#'
#' @export
#'


cellphone_for_seurat <- function(Seurat_obj){
  
  counts <- as.data.frame(
    as.matrix(
      Seurat_obj@assays$RNA@data
    )
  )
  
  #colnames(counts) <- paste('d-pos_', colnames(counts), sep = '')
  
  library("biomaRt")
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  genes  <-  getBM(filters='hgnc_symbol',
                   attributes = c('ensembl_gene_id','hgnc_symbol'),
                   values = rownames(counts),
                   mart = ensembl)
  
  counts <- counts[rownames(counts) %in% genes$hgnc_symbol,]
  
  counts <- tibble::rownames_to_column(
    as.data.frame(counts), var = 'hgnc_symbol')
  
  counts <- plyr::join(counts, genes)
  
  counts$hgnc_symbol <- NULL
  
  counts <- cbind(counts[,which(colnames(counts) == 'ensembl_gene_id')], counts)
  
  colnames(counts)[1] <- 'Gene'
  counts$ensembl_gene_id <- NULL
  
  metadata <- data.frame(Cell = rownames(Seurat_obj@meta.data),
                         cell_type = Idents(Seurat_obj)
  )
  
  #metadata$Cell <- paste('d-pos_', metadata$Cell, sep = '')
  
  write.table(counts,
              file = 'counts.txt',
              quote = F,
              col.names = T,
              row.names = F,
              sep = '\t')
  
  write.table(metadata,
              file = 'metadata.txt',
              quote = F,
              col.names = T,
              row.names = F,
              sep = '\t')
  
  system('cellphonedb method statistical_analysis metadata.txt counts.txt --iterations=10 --threads=2')
  
  system('cellphonedb plot dot_plot')
  
  system('cellphonedb plot heatmap_plot metadata.txt')
}

###########################################################################
#221011

Seurat_to_10x <- function(seurat_obj, path){
  require(Matrix)
  writeMM(seurat_obj@assays$RNA@data, file = str_c(path, '/matrix.mtx'))
  write(x = rownames(seurat_obj@assays$RNA@data), file = str_c(path, '/features.tsv'))
  write(x = colnames(seurat_obj@assays$RNA@data), file = str_c(path, '/barcodes.tsv'))
  seurat_obj@meta.data$Cell = rownames(seurat_obj@meta.data)
  write.table(seurat_obj@meta.data[, c('Cell', 'cell_type')], 
              file =str_c(path, '/metadata.tsv'), sep = '\t', quote = F, row.names = F)
}

###########################################################################
#221025
markers_park_221025 <- rev(c('ALDOB','ACAA2','PLIN2', 'SLC3A1', 
                      'CUBN', 'SPP1', 'TACSTD2', 'CLDN4',
                      'SLC12A1','UMOD','AQP2','SLC8A1',
                      'CLCNKB','SPINK1', 'ATP6V0D2', 
                      'ACTA2','TAGLN', 'PDGFRB',
                      'FLT1','PECAM1','EHD3','PLVAP',
                      'NPHS1','NPHS2','PTPRC'))

Seurat_to_cpdb <- function(seurat_obj, path){
  require(Matrix)
  writeMM(seurat_obj@assays$RNA@data, file = str_c(path, '/matrix.mtx'))
  write(x = rownames(seurat_obj@assays$RNA@data), file = str_c(path, '/features.tsv'))
  write(x = colnames(seurat_obj@assays$RNA@data), file = str_c(path, '/barcodes.tsv'))
  seurat_obj@meta.data$Cell = rownames(seurat_obj@meta.data)
  write.table(seurat_obj@meta.data[, c('Cell', 'cell_type')], 
              file =str_c(path, '/metadata.tsv'), sep = '\t', quote = F, row.names = F)
}
require(SeuratDisk)

SeuratDisk::as.h5Seurat



cpdb_dot <- function(pvalues_path, means_path, y.size = 12,
                     return_list = FALSE){
  all_pval = read_tsv(pvalues_path)
  all_means = read_tsv(means_path) 
  
  sel_pval = all_pval[,-c(1:11)]
  
  selected_columns = colnames(sel_pval)
  
  df_names = expand.grid(all_pval$interacting_pair, selected_columns)
  df_names1 = expand.grid(all_pval$gene_a, selected_columns)
  df_names2 = expand.grid(all_pval$gene_b, selected_columns)
  
  pval = unlist(sel_pval) ; pr = unlist(all_means[,-c(1:11)])
  plot_data <- cbind(df_names,df_names1$Var1, df_names2$Var1, pval, pr) %>% 
    mutate (Var2 = as.character(Var2),
            L = str_replace(Var2, '\\|.*',""),
            R = str_replace(Var2, '.*\\|',""))  %>% 
    filter(pval <0.05) %>% 
    select(-Var2) %>% 
    as_tibble()  # %>% 
  # filter(L == 'KRM')
  
  colnames(plot_data) = c('pair','sender','receiver', 'pvalue', 'mean', 'L', 'R')
  plot_data <- plot_data %>% 
    mutate(sender = as.character(sender),
           receiver = as.character(receiver))
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  cpdbdot <- ggplot(plot_data,aes(x=R,y=pair)) +
    geom_point(aes(size=mean)) +
    # scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=y.size, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  if (return_list){
    return(plot_data)
  }
  else (return(cpdbdot))
}
