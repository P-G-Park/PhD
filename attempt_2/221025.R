unlink('.Rdata'); rm(list=ls()) ; gc()
pacman::p_load(tidyverse, readxl, Seurat, data.table, ggsci, ggpubr, harmony)
source('./r_code/220811-function.R')
`%notin%` <- Negate(`%in%`)

load(file = './raw_data/processed_raw_v1.RData')

metadata <- All_Seurat_young_old_v1@meta.data %>% select(paper, channel, patient_name, patient_age, age_group)

Seurat_EC <- CreateSeuratObject(GetAssayData(All_Seurat_young_old_v1, 'count'),
                   min.cells = min_cells_per_features, min.features	= MinFeature
)
Seurat_EC$percent.mt <- PercentageFeatureSet(object = Seurat_EC, pattern = "^MT-")                   
Seurat_EC@meta.data <- cbind(Seurat_EC@meta.data, metadata)
Seurat_EC <- Seurat_EC %>% SCTransform(vars.to.regress = "percent.mt", conserve.memory = TRUE) %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = 'channel')

save(Seurat_EC, file = './raw_data/processed_raw_v1.RData')

load(file = './raw_data/processed_raw_v1.RData')

Seurat_EC <- Seurat_EC %>% clustering(resol = 0.6)
Seurat_EC <- Seurat_EC %>% subset(seurat_clusters != 13)

Seurat_EC1 <-  Seurat_EC
quick(Seurat_EC1)
quickdot(Seurat_EC1, feat = c(markers_park_221025))

Seurat_EC <-  Seurat_EC1 %>% RenameIdents('0' = 'PT', '1' = 'PT', '2' = 'PT', '3' = 'PT', '4' = 'PT', '10' = 'PT', '22' = 'PT',
                                         '20' = 'DL/tAL', 
                                         '5' = 'TAL',
                                         '16' = 'DCT',
                                         '14' = 'CD-P',
                                         '9' = 'CD-I',
                                         '19' = 'SMC/PERI',
                                         '12' = 'gEC',
                                         '7' = 'ptEC',
                                         '18' = 'PODO', 
                                         '6' = 'IMM-1', 
                                         '11' = 'IMM-2', 
                                         '17' = 'IMM-3')

quick(Seurat_EC)
quickdot(Seurat_EC, feat = markers_park_221025)

###########################################################################
Seurat_EC$HBB_Exp <- FetchData(Seurat_EC, 'HBB')
Seurat_EC_immune <- Seurat_EC %>% subset(idents = c('IMM-1','IMM-2','IMM-3')) %>% 
  subset(percent.mt <= 5) %>%
  subset(HBB_Exp <1) %>%
  clustering()
Seurat_EC_immune1 <- Seurat_EC_immune
quick(Seurat_EC_immune1)
quickdot(Seurat_EC_immune1, feat = markers_immune_park)

Seurat_EC_immune <- Seurat_EC_immune1 %>% RenameIdents(  `3` = 'KRM',
                                                        `1` = '?',
                                                        `2` = '??',
                                                        `4` = 'CD4+ T',
                                                        `5` = 'CD8+ TN',
                                                        `0` = 'CD8+ TEM',
                                                        `6` = 'B')
quick(Seurat_EC_immune)
quickdot(Seurat_EC_immune, feat = markers_immune_park)
KRM_cluster <- subset(Seurat_EC_immune, idents = 'KRM') %>% clustering()
DimPlot(KRM_cluster, group.by = 'age_group')

quick(KRM_cluster)
FeaturePlot(KRM_cluster, 'CEBPD')

Idents(KRM_cluster) <- KRM_cluster$age_group
old_DEG <- FindMarkers(KRM_cluster, ident.1 = 'old')
write.csv(Young_DEG, 'young.csv', row.names = T)


# volcano plot (221103)

EnhancedVolcano::EnhancedVolcano(old_DEG, x = 'avg_log2FC',
                                 y = 'p_val', 
                                 lab = row.names(old_DEG),
                                 pCutoff = 10e-4,
                                 col=c('black', 'black', 'black', 'red3'),
                                 drawConnectors = TRUE,
                                 title = NULL,
                                 subtitle = NULL)

save(Seurat_EC, Seurat_EC_immune, KRM_cluster, file = './raw_data/attempt_2/processed_EC_v1.RData')
###########################################################################
# EC
Seurat_EC_ptEC <- Seurat_EC %>% subset(idents = 'ptEC') %>% 
  clustering()
DimPlot(Seurat_EC_ptEC, group.by = 'age_group')
Idents(Seurat_EC_ptEC) <- Seurat_EC_ptEC$age_group
old_ptEC_DEG <- FindMarkers(Seurat_EC_ptEC, ident.1 = 'old')
old_ptEC_DEG %>% arrange(desc(avg_log2FC))
write.csv(old_ptEC_DEG, './output_csv_data/221103_ptEC_DEG.csv', row.names = T)

EnhancedVolcano::EnhancedVolcano(old_ptEC_DEG, x = 'avg_log2FC',
                                 y = 'p_val', 
                                 lab = row.names(old_ptEC_DEG),
                                 pCutoff = 10e-4,
                                 col=c('black', 'black', 'black', 'red3'),
                                 drawConnectors = TRUE,
                                 title = NULL,
                                 subtitle = NULL)

pacman::p_load(fgsea, msigdbr, ggsci)

old_gene <- old_ptEC_DEG$avg_log2FC
names(old_gene) <- rownames(old_ptEC_DEG)

h_pathway <- msigdbr(species = "human", category = "H")
h_list <- split(x = h_pathway$gene_symbol, f = h_pathway$gs_name)

gsea_H <- fgsea(pathways=h_list, stats=old_gene, nperm=10000) %>% 
  as_tibble() %>% 
  arrange(padj)
gsea_H <- gsea_H %>% 
  mutate(pathway_new = str_replace(pathway, 'HALLMARK_',''),
         pathway_new = str_replace(pathway_new, '_',' '),
         pathway_new = str_replace(pathway_new, '_',' '),
         pathway_new = str_replace(pathway_new, '_',' '),
         significant = if_else(pval < 0.05, 'p Value < 0.05', 'p value > 0.05'),
         Edge = unlist(lapply(gsea_H$leadingEdge, function(x)str_c(x, collapse=', ')))) %>% 
  arrange(desc(NES))
ggplot(gsea_H, aes(reorder(pathway_new, NES), NES)) +
  geom_col(aes(fill=significant)) +
  scale_fill_manual(values = c(pal_nejm(alpha = 0.8)(1),pal_nejm(alpha = 0.9)(1)),
                    breaks = c('p Value < 0.05')) +
  coord_flip() +
  labs(
    title="Hallmark pathways NES from GSEA") + 
  guides(fill=guide_legend(title=NULL)) +
  xlab('')+
  ylab('Normalized enrichment score')+
  theme_minimal()+
  theme(legend.position = c(0.85, 0.5))

NFKB <- gsea_H$leadingEdge[[1]]
p53 <- gsea_H$leadingEdge[[3]]

Seurat_EC_ptEC <- AddModuleScore(
  object = Seurat_EC_ptEC,
  features = p53,
  ctrl = 5,
  name = 'p53'
)

Seurat_EC_ptEC <- AddModuleScore(
  object = Seurat_EC_ptEC,
  features = NFKB,
  ctrl = 5,
  name = 'NFKB'
)

ms <- tibble(p53 = Seurat_EC_ptEC$p531,
             NFKB = Seurat_EC_ptEC$NFKB1,
             age = Seurat_EC_ptEC$age_group)

ggplot(ms, aes(x = p53, y = NFKB, col = age)) +
  geom_jitter(width = 0.1, height = 0.1) +
  theme_pubr(base_size = 12, legend = 'right') +
  guides(col=guide_legend(title= 'Module score'))

set.seed(123)
ggplot(ms, aes(x = p53, y = NFKB)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = age),
                  bins = 9) +
  theme_pubr(base_size = 12, legend = 'right') +
  guides(alpha = "none",
         fill = guide_legend(title= NULL))

###########################################################################
# for cpdb
load('./raw_data/attempt_2/processed_EC_v1.RData')

Idents_All <- as.character(Idents(Seurat_EC))
Imm_Order <- match(names(Idents(Seurat_EC_immune)), names(Idents(Seurat_EC)))
Idents_All[Imm_Order] <- as.character(Idents(Seurat_EC_immune))

Seurat_EC_cpdb <- Seurat_EC
Seurat_EC_cpdb$cell_type <- Idents_All

Seurat_EC_cpdb <- subset(Seurat_EC_cpdb, cell_type %notin% c('IMM-1','IMM-2','IMM-3','IMM-4'))
Idents(Seurat_EC_cpdb) <- Seurat_EC_cpdb$cell_type

Seurat_EC_cpdb_old <- subset(Seurat_EC_cpdb, age_group == "old")
Seurat_EC_cpdb_young <- subset(Seurat_EC_cpdb, age_group == "young")

# Seurat_to_10x(Seurat_EC_cpdb_old, path = './output_csv_data/221025_cpdb_old')
# Seurat_to_10x(Seurat_EC_cpdb_young, path = './output_csv_data/221025_cpdb_young')

Seurat_to_10x(Seurat_EC_cc_old, path = './output_csv_data/221025_cpdb_old')
Seurat_to_10x(Seurat_EC_cc_young, path = './output_csv_data/221025_cpdb_young')


###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################


###########################################################################
# for cellchat

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)


nojam <- c('MAC/Mono', 'Neutrophil', 'B', 'CD8+ TN', 'CD8+ TEM', 'CD4+ T')
Seurat_EC_cc <- Seurat_EC_cpdb %>% subset(cell_type %notin% nojam)


Seurat_EC_cc_old <- subset(Seurat_EC_cc, age_group == 'old')
Seurat_EC_cc_young <- subset(Seurat_EC_cc, age_group == 'young')

# old

data_old <- GetAssayData(Seurat_EC_cc_old, assay = "RNA", slot = "data")
meta_old <- data.frame(labels =  Idents(Seurat_EC_cc_old), 
                       row.names = names(Idents(Seurat_EC_cc_old))) 
cellchat_old <- createCellChat(object = data_old, meta = meta_old, group.by = "labels")
cellchat_old@DB <- CellChatDB.human

cellchat_old <- cellchat_old %>% subsetData() %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  projectData(PPI.human)
cellchat_old <- cellchat_old %>% computeCommunProb(raw.use = FALSE) %>%  
  filterCommunication(min.cells = 10) %>% 
  computeCommunProbPathway() %>%
  aggregateNet() %>% 
  netAnalysis_computeCentrality() %>% 
  identifyCommunicationPatterns(pattern = "outgoing", k = 3)

# young
data_young <- GetAssayData(Seurat_EC_cc_young, assay = "RNA", slot = "data")
meta_young <- data.frame(labels =  Idents(Seurat_EC_cc_young), 
                         row.names = names(Idents(Seurat_EC_cc_young))) 
cellchat_young <- createCellChat(object = data_young, meta = meta_young, group.by = "labels")
cellchat_young@DB <- CellChatDB.human

cellchat_young <- cellchat_young %>% subsetData() %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  projectData(PPI.human)
cellchat_young <- cellchat_young %>% 
  computeCommunProb(raw.use = FALSE) %>%  
  filterCommunication(min.cells = 10) %>% 
  computeCommunProbPathway() %>%
  aggregateNet() %>% 
  netAnalysis_computeCentrality() %>% 
  identifyCommunicationPatterns(pattern = "outgoing",k = 3)

# all

data_oldyoung <- GetAssayData(Seurat_EC_cc, assay = "RNA", slot = "data")
meta_oldyoung <- data.frame(labels =  Idents(Seurat_EC_cc), 
                       row.names = names(Idents(Seurat_EC_cc))) 
cellchat_oldyoung <- createCellChat(object = data_oldyoung, meta = meta_oldyoung, group.by = "labels")
cellchat_oldyoung@DB <- CellChatDB.human

cellchat_oldyoung <- cellchat_oldyoung %>% subsetData() %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  projectData(PPI.human)
cellchat_oldyoung <- cellchat_oldyoung %>% computeCommunProb(raw.use = FALSE) %>%  
  filterCommunication(min.cells = 10) %>% 
  computeCommunProbPathway() %>%
  aggregateNet() %>% 
  netAnalysis_computeCentrality() %>% 
  identifyCommunicationPatterns(pattern = "outgoing", k = 3)


netAnalysis_signalingRole_scatter(cellchat_oldyoung)

###########################################################################
save(Seurat_EC_cc_old, Seurat_EC_cc_young, 
     cellchat_old, cellchat_young, file = './raw_data/attempt_2/221025_cellchat.RData')



###########################################################################
unlink('.Rdata'); rm(list=ls()) ; gc()
load(file = './raw_data/attempt_2/221025_cellchat.RData')

object.list <- list(old = cellchat_old, young = cellchat_young)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

netVisual_diffInteraction(cellchat, weight.scale = T, measure = 'weight')
netVisual_diffInteraction(cellchat, sources.use = 'KRM', weight.scale = T, measure = "weight")

netAnalysis_signalingRole_scatter(cellchat_old)
netAnalysis_signalingRole_scatter(cellchat_young)

rankNet(cellchat, mode = "comparison", sources.use = 'KRM', targets.use = 'ptEC',
        stacked = T, do.stat = TRUE, measure = 'weight')
rankNet(cellchat, mode = "comparison", sources.use = 'KRM', targets.use = 'gEC',
        stacked = T, do.stat = TRUE, measure = 'weight')

netVisual_chord_gene(cellchat_old, sources.use = 'KRM', targets.use = 'ptEC', signaling = 'IL4')
netVisual_aggregate(cellchat_old, signaling = 'IL4', layout = "circle")

netVisual_aggregate(cellchat_old, signaling = 'VISTA', layout = "circle")
netVisual_aggregate(cellchat_young, signaling = 'VISTA', layout = "circle")

netVisual_aggregate(cellchat_old, signaling = 'TGFb', layout = "circle")
netVisual_aggregate(cellchat_young, signaling = 'TGFb', layout = "circle")


netVisual_chord_gene(cellchat_old, sources.use = 'KRM', targets.use = 'ptEC', signaling = 'VISTA')
netVisual_chord_gene(cellchat_young, sources.use = 'KRM', targets.use = 'ptEC', signaling = 'VISTA')


####################################################################
# cpdb
require(tidyverse)
require(stringr)

lr_old <- cpdb_dot(means_path = './output_csv_data/221025_cpdb_old/out/means.txt',
                   pvalues_path = './output_csv_data/221025_cpdb_old/out/pvalues.txt',
                   return_list = TRUE)
lr_young <- cpdb_dot(means_path = './output_csv_data/221025_cpdb_young/out/means.txt',
                     pvalues_path = './output_csv_data/221025_cpdb_young/out/pvalues.txt',
                     return_list = TRUE)

lr_old %>% filter(L == 'KRM', R == 'ptEC')
lr_young %>% filter(L == 'KRM', R == 'ptEC')

lr_old %>% filter(L == 'KRM', R == 'gEC')
lr_young %>% filter(L == 'KRM', R == 'gEC')

lr_old %>% filter(L == 'KRM') %>% select(R) %>% table()
lr_young %>% filter(L == 'KRM') %>% select(R) %>% table()

KRM_EC_old <- lr_old %>% filter(L=='KRM', R=='ptEC') 
KRM_EC_young <- lr_young %>% filter(L=='KRM', R=='ptEC') 


lr_old %>% filter(L == 'KRM') %>% arrange(sender) %>% View()
lr_young %>% filter(L == 'KRM') %>% arrange(sender) %>% View()

## circos old
require(circlize)
circos.clear()

Group <- c(structure(KRM_EC_old$L, names = KRM_EC_old$sender)[!duplicated(KRM_EC_old$sender)], 
           structure(KRM_EC_old$R, names = KRM_EC_old$receiver)[!duplicated(KRM_EC_old$receiver)]
)
grid_col <- dplyr::recode(Group, KRM = pal_nejm(alpha=0.8)(2)[1],
                   ptEC = pal_nejm(alpha=0.8)(2)[2])


chordDiagram(tibble(KRM_EC_old$sender, KRM_EC_old$receiver, KRM_EC_old$mean),
             #link.arr.type = 'big.arrow',
             col = 'white',
             directional = 1,
             direction.type = 'arrows',
             link.arr.col = 'black', link.arr.length = 0.2,
             annotationTrack = 'grid',
             grid.col = grid_col,
             preAllocateTracks = list(track.height = mm_h(4),
                                      track.margin = c(mm_h(24), 0)),
             group = Group)
circos.track(ylim = c(0,1), track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] - mm_y(24), CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 1))
}, bg.border = NA)

highlight.sector(KRM_EC_old$sender, track.index = 1, col = pal_nejm(alpha=0.8)(2)[1], 
                 text = "KRM", cex = 0.8, text.col = "white", niceFacing = TRUE)
highlight.sector(KRM_EC_old$receiver, track.index = 1, col = pal_nejm(alpha=0.8)(2)[2], 
                 text = "ptEC", cex = 0.8, text.col = "white", niceFacing = TRUE)

# 221114 
old_krm <-  lr_old %>% filter (L == 'KRM') 
young_krm <- lr_young %>% filter(L == 'KRM')
ord <- c('KRM', 'PT','DL/tAL','TAL','DCT/CD-P','CD-I','SMC/PERI','gEC','ptEC','PODO')

ggplot(old_krm %>% filter(pair %notin% young_krm$pair, R != 'KRM'),aes(x=factor(R, level = ord), y=pair)) +
  geom_point(aes(size=mean)) +
  # scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

####################################################################
# scenic
setwd('./python')
SeuratDisk::as.loom(Seurat_EC_cc_old, filename = 'old.loom')
SeuratDisk::as.loom(Seurat_EC_cc_young, filename = 'young.loom')

###########################################################################

FeaturePlot(Seurat_EC, 'IL13')
FeaturePlot(Seurat_EC, 'IL13RA1')

Seurat_EC_old <- subset(Seurat_EC, age_group == 'old')
Seurat_EC_young <- subset(Seurat_EC, age_group == 'young')

FeaturePlot(Seurat_EC_old, 'TGFB1') + FeaturePlot(Seurat_EC_young, 'TGFB1')
FeaturePlot(Seurat_EC_old, 'TGFBR3') + FeaturePlot(Seurat_EC_young, 'TGFBR3')


cellchat_old@net$prob

names(cellchat_old@net)
dim(cellchat_old@net$pval)
dim(cellchat_old@net$prob)

signal_old <- subsetCommunication(cellchat_old) %>% 
  filter(source == 'KRM', target == 'ptEC') %>% arrange(desc(prob))
signal_young <- subsetCommunication(cellchat_young) %>% 
  filter(source == 'KRM', target == 'ptEC') %>% arrange(desc(prob))

FeaturePlot(Seurat_EC, 'C1QC')


FeaturePlot(Seurat_EC, 'MIF')
FeaturePlot(Seurat_EC, 'CD74')
FeaturePlot(Seurat_EC, 'SPP1')
FeaturePlot(Seurat_EC, 'IL13')

FeaturePlot(Seurat_EC_cc_old, 'GAS6')
FeaturePlot(Seurat_EC_cc_young, 'GAS6')

FeaturePlot(Seurat_EC_immune, 'IL13RA1')
FeaturePlot(Seurat_EC_immune, 'TGFB1')
FeaturePlot(Seurat_EC_immune, 'TGFBR3')
FeaturePlot(Seurat_EC_immune, 'TNFRSF1A')

writexl::write_xlsx(list(old = signal_old, young = signal_young), '221102_cellchat.xlsx')


####################################################################
# 221123_ver3
FeaturePlot(Seurat_EC, 'CEBPD')
FeaturePlot(Seurat_EC_immune, 'CEBPD')

load(file = './raw_data/attempt_2/221025_cellchat.RData')

netVisual_chord_gene(cellchat_old, sources.use = 'KRM', targets.use = 'ptEC', signaling = 'COMPLEMENT')
netVisual_aggregate(cellchat_old, signaling = 'COMPLEMENT', layout = "circle")

netVisual_chord_gene(cellchat_young, sources.use = 'KRM', targets.use = 'ptEC', signaling = 'COMPLEMENT')
netVisual_aggregate(cellchat_young, signaling = 'COMPLEMENT', layout = "circle")

FeaturePlot(Seurat_EC, 'C3')
FeaturePlot(Seurat_EC_immune, 'C3')

FeaturePlot(Seurat_EC, 'CR2')
FeaturePlot(Seurat_EC_immune, 'CR2')

FeaturePlot(Seurat_EC, 'CR1')
FeaturePlot(Seurat_EC_immune, 'CR1')

FeaturePlot(Seurat_EC, 'ITGAX')
FeaturePlot(Seurat_EC_immune, 'ITGAX')

FeaturePlot(Seurat_EC, 'ITGB2')
FeaturePlot(Seurat_EC_immune, 'ITGB2')

cellchat_old@net$prob['KRM',,'C3_ITGAX_ITGB2']
cellchat_young@net$prob['KRM',,'C3_ITGAX_ITGB2']

cellchat_old@net$prob['KRM',,'C3_CR2']
cellchat_young@net$prob['KRM',,'C3_CR2']

old_chat <- subsetCommunication(cellchat_old, thresh = 0.5) %>% filter(source == 'KRM', target == 'ptEC')
young_chat <- subsetCommunication(cellchat_young, thresh = 0.5) %>% filter(source == 'KRM', target == 'ptEC')

####################################################################

# 221110 susztak lr analysis

# all_markers <- FindAllMarkers(Seurat_EC_cpdb)
write.csv(hello, 'all_markers.csv', row.names = TRUE)

LRdb <- read.csv('./raw_data/LRdb_122019.txt', sep = '\t') %>% tibble()

cell_exp_old <- tibble()
for (i in levels(Seurat_EC_cc_old)){
  seurat_i <- subset(Seurat_EC_cc_old, idents = i)
  MAT <- GetAssayData(object = seurat_i, slot = "counts")
  MAT <-  MAT > 0
  cell_exp_old <- bind_rows(cell_exp_old, tibble (cell = i, gene = names(rowMeans(MAT)[rowMeans(MAT) > 0.2])))
}

krm_gene_old <- cell_exp_old %>% filter(cell == 'KRM') %>% arrange(gene) %>% .$gene
ptec_gene_old <- cell_exp_old %>% filter(cell == 'ptEC') %>% arrange(gene) %>% .$gene

LR_KRM_old <- LRdb %>% filter(ligand %in% krm_gene_old, receptor %in% ptec_gene_old)

load(file = './raw_data/attempt_2/221025_cellchat.RData')

DotPlot(object = Seurat_EC_cpdb_old, features = unique(LR_KRM_old$ligand))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
DotPlot(object = Seurat_EC_cpdb_young, features = unique(LR_KRM_old$ligand))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

DotPlot(object = Seurat_EC_cpdb_old, features = unique(LR_KRM_old$receptor))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
DotPlot(object = Seurat_EC_cpdb_young, features = unique(LR_KRM_old$receptor))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))


# young

cell_exp_young <- tibble()
for (i in levels(Seurat_EC_cc_young)){
  seurat_i <- subset(Seurat_EC_cc_young, idents = i)
  MAT <- GetAssayData(object = seurat_i, slot = "counts")
  MAT <-  MAT > 0
  cell_exp_young <- bind_rows(cell_exp_young, tibble (cell = i, gene = names(rowMeans(MAT)[rowMeans(MAT) > 0.2])))
}

krm_gene_young <- cell_exp_young %>% filter(cell == 'KRM') %>% arrange(gene) %>% .$gene
ptec_gene_young <- cell_exp_young %>% filter(cell == 'ptEC') %>% arrange(gene) %>% .$gene

LR_KRM_young <- LRdb %>% filter(ligand %in% krm_gene_young, receptor %in% ptec_gene_young)


DotPlot(object = Seurat_EC_cpdb_old, features = unique(LR_KRM_young$ligand))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
DotPlot(object = Seurat_EC_cpdb_young, features = unique(LR_KRM_young$ligand))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

DotPlot(object = Seurat_EC_cpdb_old, features = unique(LR_KRM_young$receptor))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
DotPlot(object = Seurat_EC_cpdb_young, features = unique(LR_KRM_young$receptor))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

####################################################################
# 221110 continuous variable

load(file = './raw_data/attempt_2/processed_raw_v1.RData')

Seurat_EC_age <- All_Seurat_v1 %>% clustering(resol = 0.6); quickdot(Seurat_EC_age)
Seurat_immune_age <- Seurat_EC_age %>% subset(idents = c(4, 9, 10, 12)) %>% clustering(); 

seurat_endo <- Seurat_EC_age %>% subset(idents = c(3, 13)) %>% clustering(0.1)
seurat_ptec <- seurat_endo %>% subset(idents = 0)

quickdot(Seurat_immune_age, feat = markers_immune_park)
KRM_age <- subset(Seurat_immune_age, idents = 2)

# KRM_age <- KRM_cluster

require(limma)
lmMat <- data.frame(row.names = rownames(KRM_age))
Age <- vector()
for (i in unique(KRM_age$patient_name)){
  KRM_i <- subset(KRM_age, patient_name == i)
  MAT <- GetAssayData(object = KRM_i, slot = "data")
  df <- data.frame(rowMeans(MAT))
  colnames(df) <- i
  lmMat <- cbind(lmMat, df)
  Age <- c(Age, KRM_i$patient_age[1])
}

ggplot(tibble(age = Age, feature = unlist(lmMat['CEBPD',])), aes(x=age, y=feature)) +
  geom_jitter(width = 0.3) +
  theme_pubr(base_size = 12, legend = 'right') +
  stat_smooth(method = 'lm', se=T, color='black') + 
  stat_cor(method = "pearson",  label.x = 3, label.y = 2)+
  ylab('CEBPD expression')+
  xlab('Age')

feature_pval <- tibble()
for (i in 1:nrow(lmMat)){
  ps_test <- cor.test(unlist(lmMat[i,]), Age)
  tb <- tibble(feature = rownames(lmMat)[i], pval = ps_test$p.value, corr = ps_test$estimate)
  feature_pval <- bind_rows(feature_pval, tb)
}
sig <- feature_pval %>% filter(pval < 0.05, corr > 0) %>% 
  arrange(pval)

# ptec
# lmMat <- data.frame(row.names = rownames(seurat_ptec))
# Age <- vector()
# for (i in unique(seurat_ptec$patient_name)){
#   ptec_i <- subset(seurat_ptec, patient_name == i)
#   MAT <- GetAssayData(object = ptec_i, slot = "data")
#   df <- data.frame(rowMeans(MAT))
#   colnames(df) <- i
#   lmMat <- cbind(lmMat, df)
#   Age <- c(Age, ptec_i$patient_age[1])
# }
# 
# ggplot(tibble(age = Age, C3 = unlist(lmMat['ITGB2',])), aes(x=age, y=C3)) +
#   geom_jitter(width = 0.3) +
#   theme_pubr(base_size = 12, legend = 'right') +
#   stat_smooth(method = 'lm', se=F, color='black') + 
#   stat_cor(method = "pearson", label.x = 3, label.y = .3)+
#   ylab('ITGB2 expression')


# # 221116
# Sex <- tibble(patient_name = colnames(AverageExpression(Seurat_EC_age, group.by = 'patient_name', features =  'RPS4Y1')$RNA),
#               sex = if_else(AverageExpression(Seurat_EC_age, group.by = 'patient_name', features =  'RPS4Y1')$RNA > 0.1, 'M','F')
# )
# male_name <- unname(Sex$patient_name[Sex$sex == 'M'])
# female_name <- unname(Sex$patient_name[Sex$sex == 'F'])
# 
# KRM_age_M <- subset(KRM_age, patient_name %in% male_name)
# KRM_age_F <- subset(KRM_age, patient_name %in% female_name)
# 
# lmMat_M <- data.frame(row.names = rownames(KRM_age_M))
# lmMat_F <- data.frame(row.names = rownames(KRM_age_F))
# Age_M <- vector()
# Age_F <- vector()
# 
# for (i in unique(KRM_age_M$patient_name)){
#   KRM_i <- subset(KRM_age_M, patient_name == i)
#   MAT <- GetAssayData(object = KRM_i, slot = "data")
#   df <- data.frame(rowMeans(MAT))
#   colnames(df) <- i
#   lmMat_M <- cbind(lmMat_M, df)
#   Age_M <- c(Age_M, KRM_i$patient_age[1])
# }
# 
# for (i in unique(KRM_age_F$patient_name)){
#   KRM_i <- subset(KRM_age_F, patient_name == i)
#   MAT <- GetAssayData(object = KRM_i, slot = "data")
#   df <- data.frame(rowMeans(MAT))
#   colnames(df) <- i
#   lmMat_F <- cbind(lmMat_F, df)
#   Age_F <- c(Age_F, KRM_i$patient_age[1])
# }
# 
# ggplot(tibble(age = Age, feature = unlist(lmMat['CEBPD',])), aes(x=age, y=feature)) +
#   geom_jitter(width = 0.3) +
#   theme_pubr(base_size = 12, legend = 'right') +
#   stat_smooth(method = 'gam', se=T, color='black') + 
#   stat_cor(method = "pearson", label.x = 3, label.y = 2)+
#   ylab('CEBPD expression: ALL')
# 
# ggplot(tibble(age = Age_M, feature = unlist(lmMat_M['CEBPD',])), aes(x=age, y=feature)) +
#   geom_jitter(width = 0.3) +
#   theme_pubr(base_size = 12, legend = 'right') +
#   stat_smooth(method = 'gam', se=T, color='black') + 
#   stat_cor(method = "pearson", label.x = 3, label.y = 2)+
#   ylab('CEBPD expression: MALE')
# 
# ggplot(tibble(age = Age_F, feature = unlist(lmMat_F['CEBPD',])), aes(x=age, y=feature)) +
#   geom_jitter(width = 0.3) +
#   theme_pubr(base_size = 12, legend = 'right') +
#   stat_smooth(method = 'lm', se=T, color='black') + 
#   stat_cor(method = "pearson", label.x = 3, label.y = 2)+
#   ylab('CEBPD expression: FEMALE')



## 221128 scWGCNA
pacman::p_load(tidyverse, readxl, Seurat, data.table, ggsci, ggpubr,
               hdWGCNA, cowplot, patchwork, WGCNA, harmony)
theme_set(theme_cowplot()); set.seed(12345)
source('./r_code/220811-function.R')
`%notin%` <- Negate(`%in%`)

load(file = './raw_data/attempt_2/processed_raw_v1.RData')
Seurat_EC_age <- All_Seurat_v1 %>% clustering(resol = 0.6)
Seurat_immune_age <- Seurat_EC_age %>% subset(idents = c(4, 9, 10, 12)) %>% clustering() 
KRM_age <- subset(Seurat_immune_age, idents = 2)

KRM_cluster <- KRM_age

# load('./raw_data/attempt_2/processed_EC_v1.RData')
KRM_cluster <- KRM_cluster %>% SetupForWGCNA(
  gene_select = "fraction", fraction = 0.05, wgcna_name = "KRM") %>% 
  MetacellsByGroups(group.by = "patient_name", 
                    k = 3, min_cells = 2,  max_shared = 3, 
                    ident.group = 'patient_name') %>% 
  NormalizeMetacells() %>% 
  SetDatExpr(group_name = names(table(KRM_cluster$patient_name)),
             use_metacells = TRUE, group.by='patient_name') %>% 
  TestSoftPowers(networkType = 'signed', corFnc = 'cor')

plot_list <- PlotSoftPowers(KRM_cluster); wrap_plots(plot_list, ncol=2)

KRM_cluster <- ConstructNetwork(KRM_cluster, setDatExpr=FALSE, tom_name = 'KRM')

PlotDendrogram(KRM_cluster, main='Dendrogram')

KRM_cluster <- ScaleData(KRM_cluster, features=VariableFeatures(KRM_cluster)) %>% 
  ModuleEigengenes(group.by.vars="patient_name") %>% 
  ModuleConnectivity() %>% 
  ResetModuleNames(new_name = "KRM") %>% 
  ModuleExprScore(n_genes = 25, method='Seurat')

PlotKMEs(KRM_cluster)

modules <- GetModules(KRM_cluster); head(modules[,1:6])
hub_df <- GetHubGenes(KRM_cluster, n_hubs = 10);hub_df
wrap_plots(ModuleFeaturePlot(KRM_cluster, features='hMEs', order=TRUE))

MEs <- GetMEs(KRM_cluster, harmonized=TRUE); 
mods <- colnames(MEs); mods <- mods[mods != 'grey']
KRM_cluster@meta.data <- cbind(KRM_cluster@meta.data, MEs)
DotPlot(KRM_cluster, features=mods, group.by = 'age_group', scale = F)


# 221128 nichenet
pacman::p_load(nichenetr)
Seurat_EC_cpdb

ligand_target_matrix = readRDS('./ext_source/ligand_target_matrix.rds')
lr_network = readRDS('./ext_source/lr_network.rds')
weighted_networks = readRDS('./ext_source/weighted_networks.rds')
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

## receiver
receiver = "ptEC"
expressed_genes_receiver = get_expressed_genes(receiver, Seurat_EC_cpdb, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("KRM")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, Seurat_EC_cpdb, 0.10)
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#2
seurat_obj_receiver = subset(Seurat_EC_cpdb, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["age_group"]])

condition_oi = "old"
condition_reference = "young" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#3
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#4
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

DotPlot(Seurat_EC_cpdb, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

#5
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands in KRM","Predicted target genes in ptEC", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network + NoLegend()

#6
DE_table_all = Idents(Seurat_EC_cpdb) %>% levels() %>% intersect(sender_celltypes) %>% 
  lapply(get_lfc_celltype, seurat_obj = Seurat_EC_cpdb, condition_colname = "age_group", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join)
DE_table_all[is.na(DE_table_all)] = 0

ligand_activities_de = ligand_activities %>% 
  select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% 
  left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,] %>% as.matrix()

colnames(vis_ligand_lfc) = 'KRM'

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("","Log FC in KRM", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

# ligand activities
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))

p_ligand_pearson 

FeaturePlot(Seurat_EC, 'HMGB2')
FeaturePlot(Seurat_EC, 'C1QC')

FeaturePlot(Seurat_EC_cpdb_young, 'TIPIN')
FeaturePlot(Seurat_EC_cpdb_old, 'TIPIN')
