pacman::p_load(Seurat, tidyverse)

setwd('./raw_data/attempt_2/221108_st/kpmp_data')

st_21_019 <- Load10X_Spatial(data.dir = './raw_data/st/kpmp_data',
                filename = '21-019.h5',
                image = new(
                  Class = 'VisiumV1',
                  image = png::readPNG('./raw_data/st/kpmp_data/21-019.png'),
                  scale.factors = scalefactors(1,2000,1,0.1),
                  coordinates = read.csv('./raw_data/st/kpmp_data/21-019.csv', row.names = 1),
                  spot.radius = 0.03
                )) %>% 
  SCTransform(assay = "Spatial", verbose = FALSE)

subset(st_21_019, anterior1_imagerow > 500)


st_M61 <- Load10X_Spatial(data.dir = getwd(),
                             filename = 'M61.h5',
                             image = new(
                               Class = 'VisiumV1',
                               image = png::readPNG('M61.png'),
                               scale.factors = scalefactors(1,2000,1,0.1),
                               coordinates = read.csv('M61.csv', row.names = 1),
                               spot.radius = 0.014
                             )) %>% 
  SCTransform(assay = "Spatial", verbose = FALSE)

st_M32 <- Load10X_Spatial(data.dir = getwd(),
                          filename = 'M32.h5',
                          image = new(
                            Class = 'VisiumV1',
                            image = png::readPNG('M32.png'),
                            scale.factors = scalefactors(1,2000,1,0.1),
                            coordinates = read.csv('M32.csv', row.names = 1),
                            spot.radius = 0.014
                          )) %>% 
  SCTransform(assay = "Spatial", verbose = FALSE)

SpatialFeaturePlot(st_21_019, 'PLVAP')+ theme(aspect.ratio = 2.3)
SpatialFeaturePlot(st_21_019, 'CEBPD')+ theme(aspect.ratio = 2.3)
SpatialFeaturePlot(st_21_019, 'PLVAP')+ theme(aspect.ratio = 2.3)

SpatialFeaturePlot(st_M61, 'C1QC')+ theme(aspect.ratio = 1)
SpatialFeaturePlot(st_M61, 'CEBPD')+ theme(aspect.ratio = 1)
SpatialFeaturePlot(st_M32, 'CEBPD')+ theme(aspect.ratio = 1)

getMarkers <- function(ST, value){
  spatial_data <- ST@assays$SCT@data
  KRM_data <- spatial_data[,spatial_data['C1QC',] > value]
  CEBPD_hi <- KRM_data['CEBPD',] > 1 
  KRM_seurat <- CreateSeuratObject(KRM_data) 
  KRM_seurat$CEBPD_hi <- CEBPD_hi
  Idents(KRM_seurat) <- KRM_seurat$CEBPD_hi
  return(FindMarkers(KRM_seurat, ident.1 = 'TRUE'))
}

SpatialFeaturePlot(st_21_019, 'C1QC')+ theme(aspect.ratio = 2.3)
getMarkers(st_21_019, 0.5)

SpatialFeaturePlot(st_M32, 'C1QC')+ theme(aspect.ratio = 1)
getMarkers(st_M32, 0.5)
