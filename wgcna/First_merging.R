pacman::p_load(data.table, Seurat, tidyverse, GEOquery, R.utils)

####  Define QC parameters  #######
min_cells_per_features <- 3 
MinFeature <- 200
MaxFeature <- 2500
MaxMTpercent <- 30
num_HVG <- 2000

all_list <- list()

#1. Arazi_2019
arazi_dir <- './raw_data/wgcna/scRNA_data/sdy997'
arazi_meta <- fread(str_c(arazi_dir, '/SDY997_EXP15176_celseq_meta.tsv.725704.gz')) %>% 
  mutate(disease = recode(disease, `Control` = 'normal', `SLE` = 'ln'))
arazi_exp <- fread(str_c(arazi_dir, '/SDY997_EXP15176_celseq_matrix_ru1_molecules.tsv.725705.gz'))

arazi_mat <- arazi_exp[,-1] %>% as.matrix(); arazi_mat[is.na(arazi_mat)] <- 0
arazi_mat <- arazi_mat %>% as('sparseMatrix'); row.names(arazi_mat) <- arazi_exp[[1]]

seurat_arazi <- CreateSeuratObject(arazi_mat,
                                   meta.data = data.frame(disease_status = arazi_meta$disease, 
                                                          row.names = arazi_meta$cell_name),
                                   min.cells = min_cells_per_features, 
                                   min.features	= MinFeature)

seurat_arazi$disease_status %>% table()
seurat_arazi$paper <- 'arazi_2019'
seurat_arazi$channel <- 'arazi_2019'

all_list$arazi_2019 <- seurat_arazi

#2. Liao_2020
liao_patient_gsm <- c('GSM4145204', 'GSM4145205', 'GSM4145206')
for (i in 1:length(liao_patient_gsm)){
  gsm_i <- liao_patient_gsm[[i]]
  liao_dir <- str_c('./raw_data/wgcna/scRNA_data/Liao_2020', '/', gsm_i)
  dir.create(liao_dir, recursive = TRUE)
  getGEOSuppFiles(gsm_i, makeDirectory = FALSE, baseDir = liao_dir)
  for (j in list.files(liao_dir)){if(str_detect(j, ".*_")){
    message(str_c('executing ', j, ' renameing'))
    file.rename(str_c(liao_dir, '/', j), str_c(liao_dir, '/', str_replace(j, ".*_", "")))
  }}
  Matrix_i <- Read10X(liao_dir, gene.column = 2)
  seurat_i <- CreateSeuratObject(Matrix_i, min.cells = min_cells_per_features, 
                                 min.features	= MinFeature
  )
  seurat_i$disease_status <- 'normal'
  seurat_i$patient_name <- gsm_i
  seurat_i$paper <- "liao_2020"
  seurat_i$channel <- gsm_i
  all_list[[as.character(gsm_i)]] <- seurat_i
}

#3. Zheng_2020
zheng_gse <- 'GSE127136'
zheng_dir <- str_c('./raw_data/wgcna/scRNA_data/', zheng_gse); dir.create(zheng_dir, recursive = TRUE)
getGEOSuppFiles(zheng_gse, makeDirectory = FALSE, baseDir = zheng_dir)
zheng_exp <- fread(list.files(zheng_dir, full.names = TRUE))

zheng_mat <- zheng_exp[,-1] %>% as.matrix() %>% as(Class = 'sparseMatrix')
row.names(zheng_mat) <- zheng_exp[[1]]
zheng_mat <- zheng_mat[,!str_detect(colnames(zheng_mat), 'PBM')]

zheng_meta <- data.frame(disease_status = str_sub(colnames(zheng_mat), 1, -8),
                         patient_name = str_sub(colnames(zheng_mat), 1, -5),
                         channel = str_sub(colnames(zheng_mat), 1, -5),
                         row.names = colnames(zheng_mat)) %>% 
  mutate(disease_status = recode(disease_status, `IgAN` = 'igan', `NM` = 'normal'))

seurat_zheng <- CreateSeuratObject(zheng_mat,
                                   meta.data = zheng_meta,
                                   min.cells = min_cells_per_features, 
                                   min.features	= MinFeature)
seurat_zheng$paper <- 'zheng_2020'
all_list[[zheng_gse]] <- seurat_zheng

#4. young_2018
young_dir <- str_c('./raw_data/attempt_2/5_young_2018')
young_metadata <- readxl::read_xlsx('./raw_data/PTs/pt96_103/metadata.xlsx', sheet = 11)

young_matrix <- Seurat::ReadMtx(mtx = str_c(young_dir, '/matrix.mtx.gz'),
                                cells = str_c(young_dir, '/barcodes.tsv.gz'),
                                features = str_c(young_dir, '/features.tsv.gz'),
                                cell.column = 2,
                                feature.column = 3)

young_metadata <- young_metadata %>% mutate(ID = str_replace(Source, "_.*", ""),
                                            is_interest = ID %in% c('pRCC','RCC1','RCC2','RCC3','Trans','Wilms1','Wilms2','Wilms3'),
                                            is_normal = !str_detect(young_metadata$Source, "_T_"),
                                            is_pass = (QCpass & is_interest & is_normal)
)

seurat_young <- CreateSeuratObject(young_matrix)
seurat_young$is_pass <- young_metadata$is_pass
seurat_young$paper <- "Young_2018"
seurat_young$patient_name <- young_metadata$ID
seurat_young$channel <- young_metadata$SangerID
seurat_young$disease_status <- 'normal'
seurat_young <- subset(seurat_young, is_pass == TRUE)

all_list$seurat_young <- seurat_young


#5. Suryawanshi_2022

suryawanshi_gse <- 'GSE151671'
suryawanshi_dir <- str_c('./raw_data/wgcna/scRNA_data/', suryawanshi_gse); dir.create(suryawanshi_dir, recursive = TRUE)

suryawanshi_gsm <- c('GSM4587971', 'GSM4587972', 'GSM4587973')
suryawanshi_disease_status <- c('normal', 'allograft','allograft')

for (i in 1:length(suryawanshi_gsm)){
  gsm_i <- suryawanshi_gsm[[i]]; disease_i <- suryawanshi_disease_status[[i]]
  dir_i <- str_c(suryawanshi_dir, '/', gsm_i)
  dir.create(dir_i, recursive = TRUE)
  getGEOSuppFiles(gsm_i, makeDirectory = FALSE, baseDir = dir_i)
  
  suryawanshi_matrix_gz <- fread(list.files(dir_i, full.names = TRUE)[1])
  suryawanshi_matrix <- suryawanshi_matrix_gz[,-1] %>% as.matrix() %>% as(Class = 'sparseMatrix')
  row.names(suryawanshi_matrix) <- suryawanshi_matrix_gz[[1]]
  
  seurat_suryawanshi <- CreateSeuratObject(suryawanshi_matrix, min.cells = min_cells_per_features, 
                                           min.features	= MinFeature)
  
  seurat_suryawanshi <- CreateSeuratObject(suryawanshi_matrix, min.cells = min_cells_per_features, 
                                           min.features	= MinFeature)
  seurat_suryawanshi$disease_status <- disease_i
  seurat_suryawanshi$patient_name <- gsm_i
  seurat_suryawanshi$paper <- 'Suryawanshi_2022'
  seurat_suryawanshi$channel <- gsm_i
  all_list[[as.character(gsm_i)]] <- seurat_suryawanshi
}

#6. Tange_2021
tange_gse <- 'GSE171314'
tange_dir <- str_c('./raw_data/wgcna/scRNA_data/', tange_gse); dir.create(tange_dir, recursive = TRUE)

tange_gsm <- c('GSM5222730', 'GSM5222731', 'GSM5222732', 'GSM5222733', 'GSM5222734')
tange_disease_status <- c(rep('igan', 4), 'normal')

for (i in 1:length(tange_gsm)){
  gsm_i <- tange_gsm[[i]]; disease_i <- tange_disease_status[[i]]
  dir_i <- str_c(tange_dir, '/', gsm_i)
  dir.create(dir_i, recursive = TRUE)
  getGEOSuppFiles(gsm_i, makeDirectory = FALSE, baseDir = dir_i)
  
  tange_matrix_gz <- fread(list.files(dir_i, full.names = TRUE))
  tange_matrix <- tange_matrix_gz[,-1] %>% as.matrix() %>% as(Class = 'sparseMatrix')
  row.names(tange_matrix) <- tange_matrix_gz[[1]]
  
  seurat_tange <- CreateSeuratObject(tange_matrix, min.cells = min_cells_per_features, 
                                           min.features	= MinFeature)
  
  seurat_tange <- CreateSeuratObject(tange_matrix, min.cells = min_cells_per_features, 
                                           min.features	= MinFeature)
  seurat_tange$disease_status <- disease_i
  seurat_tange$patient_name <- gsm_i
  seurat_tange$paper <- 'tange_2021'
  seurat_tange$channel <- gsm_i
  all_list[[as.character(gsm_i)]] <- seurat_tange
}

#7. xu_2021
xu_gse <- 'GSE171458'
xu_dir <- str_c('./raw_data/wgcna/scRNA_data/', xu_gse); dir.create(xu_dir, recursive = TRUE)

xu_gsm <- c('GSM5225900', 'GSM5225901', 'GSM5225902', 'GSM5225903', 'GSM5225904', 'GSM5225905', 'GSM5225906', 'GSM5225907')
xu_disease_status <- c(rep('mn', 6), 'normal', 'normal')

for (i in 1:length(xu_gsm)){
  gsm_i <- xu_gsm[[i]]; disease_i <- xu_disease_status[[i]]
  dir_i <- str_c(xu_dir, '/', gsm_i)
  dir.create(dir_i, recursive = TRUE)
  getGEOSuppFiles(gsm_i, makeDirectory = FALSE, baseDir = dir_i)
  
  xu_matrix_gz <- fread(list.files(dir_i, full.names = TRUE))
  xu_matrix <- xu_matrix_gz[,-1] %>% as.matrix() %>% as(Class = 'sparseMatrix')
  row.names(xu_matrix) <- xu_matrix_gz[[1]]
  
  seurat_xu <- CreateSeuratObject(xu_matrix, min.cells = min_cells_per_features, 
                                     min.features	= MinFeature)
  
  seurat_xu <- CreateSeuratObject(xu_matrix, min.cells = min_cells_per_features, 
                                     min.features	= MinFeature)
  seurat_xu$disease_status <- disease_i
  seurat_xu$patient_name <- gsm_i
  seurat_xu$paper <- 'xu_2021'
  seurat_xu$channel <- gsm_i
  all_list[[as.character(gsm_i)]] <- seurat_xu
}

#7. Yu_2022
yu_gsm <- 'GSM5988612'
yu_dir <- str_c('./raw_data/wgcna/scRNA_data/', yu_gsm); dir.create(yu_dir, recursive = TRUE)
yu_disease_status <- 'mn'
getGEOSuppFiles(yu_gsm, makeDirectory = FALSE, baseDir = yu_dir)
yu_matrix_gz <- fread(list.files(yu_dir, full.names = TRUE))
yu_matrix <- yu_matrix_gz[,-1] %>% as.matrix() %>% as(Class = 'sparseMatrix')
row.names(yu_matrix) <- yu_matrix_gz[[1]]

seurat_yu <- CreateSeuratObject(yu_matrix, min.cells = min_cells_per_features, 
                                min.features	= MinFeature)

seurat_yu <- CreateSeuratObject(yu_matrix, min.cells = min_cells_per_features, 
                                min.features	= MinFeature)
seurat_yu$disease_status <- yu_disease_status
seurat_yu$patient_name <- yu_gsm
seurat_yu$paper <- 'yu_2022'
seurat_yu$channel <- yu_gsm
all_list[[as.character(yu_gsm)]] <- seurat_yu


#8. Zhong_2021
zhong_gse <- 'GSE174220'
zhong_dir <- str_c('./raw_data/wgcna/scRNA_data/', zhong_gse); dir.create(zhong_dir, recursive = TRUE)

zhong_gsm <- c('GSM5289544', 'GSM5289545')
zhong_disease_status <- c(rep('aki', 2))

for (i in 1:length(zhong_gsm)){
  gsm_i <- zhong_gsm[[i]]; disease_i <- zhong_disease_status[[i]]
  dir_i <- str_c(zhong_dir, '/', gsm_i)
  dir.create(dir_i, recursive = TRUE)
  getGEOSuppFiles(gsm_i, makeDirectory = FALSE, baseDir = dir_i)
  
  zhong_matrix_gz <- fread(list.files(dir_i, full.names = TRUE))
  zhong_matrix <- zhong_matrix_gz[,-1] %>% as.matrix() %>% as(Class = 'sparseMatrix')
  row.names(zhong_matrix) <- zhong_matrix_gz[[1]]
  
  seurat_zhong <- CreateSeuratObject(zhong_matrix, min.cells = min_cells_per_features, 
                                  min.features	= MinFeature)
  
  seurat_zhong <- CreateSeuratObject(zhong_matrix, min.cells = min_cells_per_features, 
                                  min.features	= MinFeature)
  seurat_zhong$disease_status <- disease_i
  seurat_zhong$patient_name <- gsm_i
  seurat_zhong$paper <- 'zhong_2021'
  seurat_zhong$channel <- gsm_i
  all_list[[as.character(gsm_i)]] <- seurat_zhong
}

#9. Tang_2021
tang_gse <- 'GSE174219'
tang_dir <- str_c('./raw_data/wgcna/scRNA_data/', tang_gse); dir.create(tang_dir, recursive = TRUE)

tang_gsm <- c('GSM5289541', 'GSM5289542', 'GSM5289543')
tang_disease_status <- c('htn','normal','htn')

for (i in 1:length(tang_gsm)){
  gsm_i <- tang_gsm[[i]]; disease_i <- tang_disease_status[[i]]
  dir_i <- str_c(tang_dir, '/', gsm_i)
  dir.create(dir_i, recursive = TRUE)
  getGEOSuppFiles(gsm_i, makeDirectory = FALSE, baseDir = dir_i)
  
  tang_matrix_gz <- fread(list.files(dir_i, full.names = TRUE))
  tang_matrix <- tang_matrix_gz[,-1] %>% as.matrix() %>% as(Class = 'sparseMatrix')
  row.names(tang_matrix) <- tang_matrix_gz[[1]]
  
  seurat_tang <- CreateSeuratObject(tang_matrix, min.cells = min_cells_per_features, 
                                     min.features	= MinFeature)
  
  seurat_tang <- CreateSeuratObject(tang_matrix, min.cells = min_cells_per_features, 
                                     min.features	= MinFeature)
  seurat_tang$disease_status <- disease_i
  seurat_tang$patient_name <- gsm_i
  seurat_tang$paper <- 'tang_2021'
  seurat_tang$channel <- gsm_i
  all_list[[as.character(gsm_i)]] <- seurat_tang
}


#10. Su_2022
su_gsm <- 'GSM4630031'
su_dir <- str_c('./raw_data/wgcna/scRNA_data/', su_gsm); dir.create(su_dir, recursive = TRUE)
su_disease_status <- 'normal'
getGEOSuppFiles(su_gsm, makeDirectory = FALSE, baseDir = su_dir)
su_matrix <- Read10X(su_dir) # tar -xvf
seurat_su <- CreateSeuratObject(su_matrix, min.cells = min_cells_per_features, 
                                min.features	= MinFeature)

seurat_su <- CreateSeuratObject(su_matrix, min.cells = min_cells_per_features, 
                                min.features	= MinFeature)
seurat_su$disease_status <- su_disease_status
seurat_su$patient_name <- su_gsm
seurat_su$paper <- 'su_2022'
seurat_su$channel <- su_gsm
all_list[[as.character(su_gsm)]] <- seurat_su

# 11. KPMP_2022

KPMP_path <- str_c('./raw_data/wgcna/scRNA_data/', 'KPMP_2022'); dir.create(KPMP_path, recursive = TRUE)
KPMP_names_file <- readxl::read_xlsx(str_c(KPMP_path, '/healthy_data_names.xlsx'), col_names = TRUE)  %>% 
  left_join(read_csv(str_c(KPMP_path, '/OpenAccessClinicalData.csv')),
            by = c("patient_name" = "Participant ID")) %>% 
  dplyr::select(2:15) %>% dplyr::rename ('age' = `Age (Years) (Binned)`) %>% 
  mutate (age = as.integer(str_sub(age, 1L,2L))+5)

KPMP_normal_path <- str_c(KPMP_path, '/normal')



for (i in list.files(KPMP_normal_path, pattern = '*.zip')){
  source('./r_code/220811-function.R')
  patient_identifier <- as.character(str_replace(as.character(i), "_.*",""))
  patient_name <- KPMP_names_file[match(i,KPMP_names_file$file_name), 'patient_name'] %>% pull()
  
  ex_dir <- str_c(KPMP_normal_path, '/',str_sub(i, 1, -5))
  unzip(str_c(KPMP_normal_path, '/', i), exdir = ex_dir)
  
  for (j in list.files(ex_dir)){
    if(!str_ends(j, '.gz')){
      message(str_c('executing ', j, ' gzip'))
      gzip(str_c(ex_dir, '/', j))
    }
  }
  for (j in list.files(ex_dir)){
    if(str_detect(j, ".*_")){
      message(str_c('executing ', j, ' renameing'))
      file.rename(str_c(ex_dir, '/', j), str_c(ex_dir, '/', str_replace(j, ".*_", "")))
    }
  }
  for (j in list.files(ex_dir)){
    if(str_detect(j, "genes")){
      message(str_c('genes ', j, ' renameing'))
      file.rename(str_c(ex_dir, '/', j), str_c(ex_dir, '/', str_replace(j, "genes", "features")))
    }
  }
  Matrix_i <- Read10X(ex_dir, gene.column = 2)
  rownames(Matrix_i) <- str_replace(rownames(Matrix_i), "hg38____", "")
  if(ncol(Matrix_i) > 300000){
    Matrix_i <- Matrix_i %>% RawMat_To_SoupXChannel_Using_DropletUtils %>% 
      SoupXChannel_To_Adjust_Matrix_Removing_Amb_mRNA_Using_SoupX
  }
  CSO <- function(Matrix_i){
    seurat_i <- CreateSeuratObject(Matrix_i, min.cells = 3, min.features	= 200)
    seurat_i$patient_name <- patient_name
    seurat_i$channel <- patient_identifier
    seurat_i$disease_status <- 'normal'
    seurat_i$paper <- "KPMP_2022"
    return(seurat_i)
  }
  all_list[[patient_identifier]] <- CSO(Matrix_i)
  if(length(list.files(ex_dir, pattern = '*.tsv$'))>1){print ('more than 2 tsv') }
}

KPMP_dn_path <- str_c(KPMP_path, '/dn')
dn_names_file <- readxl::read_xlsx(str_c(KPMP_path, '/dn_data_names.xlsx'), col_names = TRUE)  %>% 
  left_join(read_csv(str_c(KPMP_path, '/OpenAccessClinicalData.csv')),
            by = c("patient_name" = "Participant ID")) %>% 
  dplyr::select(2:15) %>% dplyr::rename ('age' = `Age (Years) (Binned)`) %>% 
  mutate (age = as.integer(str_sub(age, 1L,2L))+5)

for (i in list.files(KPMP_dn_path, pattern = '*.zip')){
  source('./r_code/220811-function.R')
  patient_name <- dn_names_file[match(i,dn_names_file$file_name), 'patient_name'] %>% pull()
  patient_identifier <- as.character(str_replace(as.character(i), "_.*",""))
  
  ex_dir <- str_c(KPMP_dn_path, '/',str_sub(i, 1, -5))
  unzip(str_c(KPMP_dn_path, '/', i), exdir = ex_dir)
  
  for (j in list.files(ex_dir)){
    if(!str_ends(j, '.gz')){
      message(str_c('executing ', j, ' gzip'))
      gzip(str_c(ex_dir, '/', j))
    }
  }
  for (j in list.files(ex_dir)){
    if(str_detect(j, ".*_")){
      message(str_c('executing ', j, ' renameing'))
      file.rename(str_c(ex_dir, '/', j), str_c(ex_dir, '/', str_replace(j, ".*_", "")))
    }
  }
  
  for (j in list.files(ex_dir)){
    if(str_detect(j, "genes")){
      message(str_c('genes ', j, ' renameing'))
      file.rename(str_c(ex_dir, '/', j), str_c(ex_dir, '/', str_replace(j, "genes", "features")))
    }
  }
  Matrix_i <- Read10X(ex_dir, gene.column = 2)
  rownames(Matrix_i) <- str_replace(rownames(Matrix_i), "hg38____", "")
  if(ncol(Matrix_i) > 300000){
    Matrix_i <- Matrix_i %>% RawMat_To_SoupXChannel_Using_DropletUtils %>% 
      SoupXChannel_To_Adjust_Matrix_Removing_Amb_mRNA_Using_SoupX
  }
  CSO <- function(Matrix_i){
    seurat_i <- CreateSeuratObject(Matrix_i, min.cells = 3, min.features	= 200)
    seurat_i$patient_name <- patient_name
    seurat_i$channel <- patient_identifier
    seurat_i$disease_status <- 'dn'
    seurat_i$paper <- "KPMP_2022"
    return(seurat_i)
  }
  all_list[[patient_identifier]] <- CSO(Matrix_i)
}

#12. McEvoy_2022
mcevoy_gse <- 'GSE202109'
mcevoy_gsm <- str_c('GSM60946', 52:70)
mcevoy_dir <- str_c('./raw_data/wgcna/scRNA_data/', mcevoy_gse); dir.create(mcevoy_dir, recursive = TRUE)
mcevoy_disease_status <- 'normal'

for (i in 1:length(mcevoy_gsm)){
  gsm_i <- mcevoy_gsm[[i]]; disease_i <- mcevoy_disease_status
  dir_i <- str_c(mcevoy_dir, '/', gsm_i)
  dir.create(dir_i, recursive = TRUE)
  getGEOSuppFiles(gsm_i, makeDirectory = FALSE, baseDir = dir_i)
  for (j in list.files(dir_i)){if(str_detect(j, ".*_")){
    message(str_c('executing ', j, ' renameing'))
    file.rename(str_c(dir_i, '/', j), str_c(dir_i, '/', str_replace(j, ".*_", "")))
  }}
  Matrix_i <- Read10X(dir_i, gene.column = 2)
  seurat_i <- CreateSeuratObject(Matrix_i, min.cells = min_cells_per_features, 
                                 min.features	= MinFeature
  )
  seurat_i$disease_status <- 'normal'
  seurat_i$patient_name <- gsm_i
  seurat_i$paper <- "mcevoy_2022"
  seurat_i$channel <- gsm_i
  all_list[[as.character(gsm_i)]] <- seurat_i
}


# merging
require(harmony)

First_Seurat <- all_list[[1]]
all_list[[1]]  <-  NULL
all_Seurat <- merge(First_Seurat, all_list)
all_Seurat$percent.mt <- PercentageFeatureSet(object = all_Seurat, pattern = "^MT-")

all_Seurat <- subset(all_Seurat, subset = nFeature_RNA > MinFeature & 
                       nFeature_RNA < MaxFeature & percent.mt < MaxMTpercent)

all_Seurat <- all_Seurat %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nFeatures = num_HVG) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = 'channel')

save(all_Seurat, file = './raw_data/wgcna/raw_v1.RData')

# immune extract
load(file = './raw_data/wgcna/raw_v1.RData')

all_Seurat <- all_Seurat %>% clustering(resol = 0.6)

quick(all_Seurat)
quickdot(all_Seurat)
save(all_Seurat, file = './raw_data/wgcna/raw_v1.RData')

all_Seurat_immune <- subset(all_Seurat, idents = c(3, 6, 10, 13, 19, 27, 29, 30))

save(all_Seurat_immune, file = './raw_data/wgcna/immune_v1.RData')

#13. KORNERSTONE

# snuh_1 <- Read10X('./raw_data/wgcna/scRNA_data/snuh_1') %>% CreateSeuratObject(min.cells = min_cells_per_features, 
#                                                                                min.features	= MinFeature)
# snuh_1$channel <- 'snuh_1'

#1. pre-immune cell sorting
snuh_2 <- Read10X('./raw_data/wgcna/scRNA_data/snuh_2') %>% CreateSeuratObject(min.cells = min_cells_per_features, 
                                                                              min.features	= MinFeature)
snuh_2$channel <- 'snuh_2'
snuh_2$disease_status <- 'normal'
snuh_2$paper <- 'KEYNOTE'
snuh_2$patient_name <- 'snuh_2'

load(file = './raw_data/wgcna/immune_v1.RData')
load(file = './raw_data/wgcna/raw_v1.RData')

all_Seurat <- all_Seurat %>% merge(snuh_2)
all_Seurat <- all_Seurat %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nFeatures = num_HVG) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = 'channel')

all_Seurat <- all_Seurat %>% clustering(resol = 0.6)
save(all_Seurat, file = './raw_data/wgcna/raw_v2.RData')

quickdot(all_Seurat)
all_Seurat_immune <- subset(all_Seurat, idents = c(2, 7, 9, 11, 12, 25, 27, 28))
all_Seurat_immune <- all_Seurat_immune %>% clustering(resol = 0.5)
quickdot(all_Seurat_immune, feat = markers_immune_park)

save(all_Seurat_immune, file = './raw_data/wgcna/immune_v2.RData')


# #2. post-immune cell sorting
# snuh_2 <- Read10X('./raw_data/wgcna/scRNA_data/snuh_2') %>% CreateSeuratObject(min.cells = min_cells_per_features, 
#                                                                                min.features	= MinFeature)
# snuh_2$channel <- 'snuh_2'
# 
# load(file = './raw_data/wgcna/immune_v2.RData')
# load(file = './raw_data/wgcna/raw_v1.RData')
# 
# all_Seurat_immune <- all_Seurat_immune %>% merge(snuh_2)
# all_Seurat_immune <- all_Seurat_immune %>% 
#   NormalizeData() %>% 
#   FindVariableFeatures(selection.method = "vst", nFeatures = num_HVG) %>% 
#   ScaleData() %>% 
#   RunPCA() %>% 
#   RunHarmony(group.by.vars = 'channel') 
# all_Seurat_immune <- all_Seurat_immune %>% clustering(resol = 0.6)
# 
# #
# d <- merge(snuh_1, snuh_2)
# d$percent.mt <- PercentageFeatureSet(object = d, pattern = "^MT-")
# d <- subset(d, subset = percent.mt < 10)
# 
# d <- d %>% 
#   NormalizeData() %>% 
#   FindVariableFeatures(selection.method = "vst", nFeatures = num_HVG) %>% 
#   ScaleData() %>% 
#   RunPCA() %>% 
#   RunHarmony(group.by.vars = 'channel')
# 
# d <- d %>% clustering()
# 
# quickdot(d)
# 
# quickdot(d, feat = markers_immune_park)
# quick(d)
# 
# table(d$channel, d$seurat_clusters)
#          

#3. only DN
load(file = './raw_data/wgcna/raw_v2.RData')

dn_all <- all_Seurat %>% subset (disease_status %in% c('dn','normal'))

dn_all <- dn_all %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nFeatures = num_HVG) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = 'channel')

dn_all <- dn_all %>% clustering(resol = 0.6)

quickdot(dn_all)
save(dn_all, file = './raw_data/wgcna/dn_v1.RData')

load(file = './raw_data/wgcna/dn_v1.RData')

dn_all$paper %>% table() 

dn_immune <- subset(dn_all, idents = c(2, 5, 8, 10, 11, 12, 22, 24, 25))
dn_immune <- dn_immune %>% clustering(resol = 0.5)

save(dn_immune, file = './raw_data/wgcna/dn_immune_v1.RData')
