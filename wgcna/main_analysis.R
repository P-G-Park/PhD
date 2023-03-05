rm(list=ls()); gc()
setwd('..')

pacman::p_load(tidyverse, readxl, Seurat, data.table, gridExtra, flextable, harmony)
source('./r_code/220811-function.R')

#########################################################################
load(file = './raw_data/wgcna/dn_v1.RData')
metadata <- dn_all@meta.data
dn_all <- CreateSeuratObject(
  GetAssayData(dn_all, 'counts'),
  min.cells = min_cells_per_features, min.features	= MinFeature
)
dn_all@meta.data <- cbind(dn_all@meta.data, metadata %>% 
                            select(patient_name, disease_status, paper, channel))
dn_all$percent.mt <- PercentageFeatureSet(object = dn_all, pattern = "^MT-")

dn_all <- dn_all %>% subset(is.na(patient_name) | patient_name %notin% c('Wilms1', 'Wilms2','Wilms3','Trans'))

dn_all$paper %>% table(useNA = 'always')

dn_all <- dn_all %>% 
  NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nFeatures = num_HVG) %>% 
  ScaleData(vars.to.regress = "percent.mt") %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = 'channel')

save(dn_all, file = './raw_data/wgcna/dn_v2.RData')

#########################################################################

load(file = './raw_data/wgcna/dn_v2.RData')
dn_all <- dn_all %>% clustering(resol = 0.3)
dn_all1 <- dn_all

quick(dn_all)
quickdot(dn_all)

dn_all1 <- dn_all %>% RenameIdents(
  '0' = 'PT',
  '1' = 'PT',
  '13' = 'PT',
  '14' = 'PT',
  '15' = 'PT',
  '18' = 'PT',
  '22' = 'PT',
  '16' = 'DL/ tAL',
  '6' = 'DL/ tAL',
  '3' = 'TAL',
  '5' = 'DCT/ CD-P',
  '8' = 'CD-I',
  '12' = 'SMC/ PERI',
  '4' = 'EC',
  '17' = 'PODO',
  '2' = 'IMM-1',
  '7' = 'IMM-2',
  '9' = 'IMM-3',
  '10' = 'IMM-4',
  '11' = 'IMM-5',
  '19' = 'IMM-6',
  '20' = 'IMM-7',
  '21' = 'IMM-8'
)

#DEG_ident(dn_all1)

quick(dn_all1)
quickdot(dn_all1)

dn_immune <- subset(dn_all1, idents = c('IMM-1',
                                        'IMM-2',
                                        'IMM-3',
                                        'IMM-4',
                                        'IMM-5',
                                        'IMM-6',
                                        'IMM-7',
                                        'IMM-8'))

save(dn_immune, file = './raw_data/wgcna/dn_immune_v2.RData')


# prepare for 14, 22 sorting (230210)

Idents_All <- as.character(Idents(dn_all1))
Imm_Order <- match(names(Idents(dn_immune)), names(Idents(dn_all1)))
Idents_All[Imm_Order] <- as.character(Idents(dn_immune1))

dn_all2 <- dn_all1
dn_all2$cell_type <- Idents_All

dn_all2 <- subset(dn_all2, cell_type %notin% c('IMM-1','IMM-2','IMM-3','IMM-4','IMM-5','IMM-6','IMM-7','IMM-8'))
Idents(dn_all2) <- dn_all2$cell_type


#########################################################################
load(file = './raw_data/wgcna/dn_immune_v2.RData')

dn_immune$percent.mt <- PercentageFeatureSet(object = dn_immune, pattern = "^MT-")
dn_immune <- dn_immune %>% subset(percent.mt <= 10) %>% clustering(resol = 1.1) 

mark <- c('HLA-DRA',  'LYZ', 'CD14', 'CD68',   'C1QC', 'MRC1','LYVE1',
          'FCN1','PLAC8',  'THBS1', 'VCAN',  'FCER1A', 'CD1C', 'S100A8','IL1R2',  
          'CD3E', 'TRAC', 'CD4', 'CD40LG','CD8A', 'TNF', 'GZMB', 'IFNG', 'IGKC', 'JCHAIN', 'CD79A','NKG7', 'KLRD1'
          )

quick(dn_immune)
quickdot(dn_immune, feat = mark)


dn_immune1 <- dn_immune %>% RenameIdents(
  '11' = 'KRM',
  '10' = 'Infiltrating Mac',
  '16' = 'Infiltrating Mac',
  '17' = 'Monocyte',
  '7' = 'Neutrophil',
  '20' = 'Neutrophil',
  '12' = 'cDC',
  '3' = 'CD4+ T',
  '2' = 'CD4+ T',
  '6' = 'CD4+ T',
  '1' = 'CD8+ T',
  '4' = 'CD8+ T, effector',
  '8' = 'CD8+ T, effector',
  '0' = 'NK',
  '9' = 'NK',
  '5' = 'B',
  '13' = 'B',
  '15' = 'B',
  '21' = 'B',
  '18' = 'Cytokine-secreting B',
  '19' = 'Proliferating cell'
)

dn_immune1$cell_type <- as.character(Idents(dn_immune1))
dn_immune1 <- dn_immune1 %>% subset(cell_type %notin% c(14, 22))

quick(dn_immune1)
quickdot(dn_immune1, feat = rev(mark))

a <- list('krm_mono' = (FindMarkers(dn_immune1, ident.1 = 'KRM', ident.2 = 'Monocyte') %>% rownames_to_column('gene')),
     'krm_infilt' = (FindMarkers(dn_immune1, ident.1 = 'KRM', ident.2 = 'Infiltrating Mac') %>% rownames_to_column('gene'))
)

writexl::write_xlsx(a, 'krm_deg.xlsx')

#DEG_ident(dn_all1)
#DEG_ident(dn_immune1)

non_immune <- c('PT','DL, tAL','TAL','DCT, CD-P','CD-I','SMC, PERI','EC', 'PODO')
immune <- c('KRM', 'Infiltrating Mac', 'Monocyte','Neutrophil', 'cDC',
            'CD4+ T', 'CD8+ T','CD8+ T, effector','NK','B')


# DEG_dn <- list()
# for (i in immune){
#   subset_seurat <- dn_immune1 %>% subset(idents = i)
#   Idents(subset_seurat) <- subset_seurat$disease_status
#   dn_deg <- FindMarkers(subset_seurat, ident.1 = 'dn', only.pos = TRUE)
#   DEG_dn[[i]] <- dn_deg %>% rownames_to_column('Gene_name')
# }
# for (i in non_immune){
#   subset_seurat <- dn_all1 %>% subset(idents = i)
#   Idents(subset_seurat) <- subset_seurat$disease_status
#   dn_deg <- FindMarkers(subset_seurat, ident.1 = 'dn')
#   DEG_dn[[i]] <- dn_deg %>% rownames_to_column('Gene_name')
# }
# 
# writexl::write_xlsx(DEG_dn, 'DEG_dn.xlsx')

################################# pathway analysis, DEGs
KRM <- subset(dn_immune1, idents = 'KRM') %>% clustering(resol = 0.3)
Idents(KRM) <- KRM$disease_status
KRM_dn_DEG <- FindMarkers(KRM, ident.1 = 'dn')
KRM_dn_DEG %>% arrange(desc(avg_log2FC))

pacman::p_load(fgsea, msigdbr, ggsci)

dn_gene <- KRM_dn_DEG$avg_log2FC
names(dn_gene) <- rownames(KRM_dn_DEG)

h_pathway <- msigdbr(species = "human", category = "H")
h_list <- split(x = h_pathway$gene_symbol, f = h_pathway$gs_name)

gsea_H <- fgsea(pathways=h_list, stats=dn_gene, nperm=100000) %>% 
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

# gsea_H %>% write_xlsx('gsea_hallmark_230213.xlsx')

ggplot(gsea_H %>% filter(significant == 'p Value < 0.05'), aes(reorder(pathway_new, NES), NES)) +
  geom_col(aes(fill=significant)) +
  scale_fill_manual(values = c(pal_nejm(alpha = 0.8)(1),pal_nejm(alpha = 0.9)(1)),
                    breaks = c('p Value < 0.05')) +
  coord_flip() +
  labs(title="Hallmark pathways NES from GSEA") + 
  guides(fill=guide_legend(title=NULL)) +
  xlab('')+
  ylab('Normalized enrichment score')+
  ggpubr::theme_pubr()+
 # theme_classic2()+
  theme(legend.position = c(0.85, 0.5))+
  NoLegend()

subset_gsea <- function(subset_seurat){
  Idents(subset_seurat) <- subset_seurat$disease_status
  subset_DEG <- FindMarkers(subset_seurat, ident.1 = 'dn')
  subset_DEG %>% arrange(desc(avg_log2FC))
  
  pacman::p_load(fgsea, msigdbr, ggsci)
  
  dn_gene <- subset_DEG$avg_log2FC
  names(dn_gene) <- rownames(subset_DEG)
  
  h_pathway <- msigdbr(species = "human", category = "H")
  h_list <- split(x = h_pathway$gene_symbol, f = h_pathway$gs_name)
  
  gsea_H <- fgsea(pathways=h_list, stats=dn_gene, nperm=100000) %>% 
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
  
  gsea_H %>% writexl::write_xlsx('gsea_hallmark_230213.xlsx')
  
  p <- ggplot(gsea_H, aes(reorder(pathway_new, NES), NES)) +
    geom_col(aes(fill=significant)) +
    scale_fill_manual(values = c(pal_nejm(alpha = 0.8)(1),pal_nejm(alpha = 0.9)(1)),
                      breaks = c('p Value < 0.05')) +
    coord_flip() +
    labs(title="Hallmark pathways NES from GSEA") + 
    guides(fill=guide_legend(title=NULL)) +
    xlab('')+
    ylab('Normalized enrichment score')+
    ggpubr::theme_pubr()+
    # theme_classic2()+
    theme(legend.position = c(0.85, 0.5))
  return(p)
}
KRM <- S4Vectors::subset(dn_immune1, idents = 'KRM')
Infilt <- S4Vectors::subset(dn_immune1, idents = 'Infiltrating Mac')
Mono <- S4Vectors::subset(dn_immune1, idents = 'Monocyte')

subset_gsea(Infilt)
subset_gsea(Mono)

Idents(KRM) <- KRM$disease_status
Idents(Infilt) <- Infilt$disease_status
Idents(Mono) <- Mono$disease_status


cytokines <- c('IL1A', 'IL1B', 'IL2', 'IL3', 'IL4', 'IL5', 'IL6', 'IL7', 'IL8', 'IL10',
               'IL11', 'IL12A', 'IL13', 'IL16', 'IL18', 'CXCL1', 'CXCL2', 'CXCL3', 
               'TNF')

phagocytosis <- c('AXL', 'MERTK', 'TYRO3', 'CD34', 'CD36', 'OLR1', 'STAB2', 'AGER', 'TIMD4')

fibrosis <- c('COL1A1', 'COL3A1', 'AREG',  'MMP9', 'TIMP1', 'VIM',
              'TGFB1', 'FN1', 'ACTA2', 'CTGF')

Heatmap_DEGs <- function(x){
  bind_cols(
  FoldChange(KRM, ident.1 = 'dn', features = x) %>% select(1),
  FoldChange(Infilt, ident.1 = 'dn', features = x) %>% select(1),
  FoldChange(Mono, ident.1 = 'dn', features = x) %>% select(1)
  ) %>% 
  `colnames<-`(c('KRM', 'Infiltrating Mac', 'Monocyte')) %>% 
  ComplexHeatmap::Heatmap(name = ' ',
                          cluster_rows = FALSE, cluster_columns = FALSE,   
                          show_column_names = TRUE,
                          column_names_rot = 30)
}

Heatmap_DEGs(cytokines)
Heatmap_DEGs(phagocytosis)
Heatmap_DEGs(fibrosis)


OxPhos <- gsea_H$leadingEdge[[1]]
DNA_repair <- gsea_H$leadingEdge[[2]]
Coag <- gsea_H$leadingEdge[[3]]

OxPhos_select <- KRM_dn_DEG %>% rownames_to_column('gene') %>% 
  filter(gene %in% OxPhos, avg_log2FC > .5, p_val_adj < 10^-5) %>% 
  pull(gene)

my_lab <- c(expression(list(Log['2']*FC > 0.5,Adj.Pval < 10^-5)))
p <- EnhancedVolcano::EnhancedVolcano(KRM_dn_DEG, x = 'avg_log2FC',
                                      y = 'p_val_adj', 
                                      lab = row.names(KRM_dn_DEG),
                                      selectLab = OxPhos_select,
                                      pCutoff = 10e-5,
                                      FCcutoff = .5,
                                      col=c('black', 'black', 'black', 'red3'),
                                      drawConnectors = TRUE,
                                      gridlines.major = F,
                                      gridlines.minor = F,
                                      title = NULL,
                                      subtitle = NULL,
                                      max.overlaps = Inf)
p + scale_colour_manual(values = c(NS = 'black', 
                                   FC = 'black', 
                                   P = 'black', 
                                   FC_P = 'red'),
                        labels = my_lab,
                        breaks = c(FC_P = 'FC_P')) +
  xlim(c(-4, 3))


KRM <- AddModuleScore(
  object = KRM,
  features = OxPhos,
  nbin = 12,
  ctrl = 100,
  name = 'OxPhos'
)
require(ggpubr)

rain_height <- .1
tibble(disease = KRM$disease_status, OxPhos = KRM$OxPhos1) %>%
  mutate(disease = fct_relevel(disease, c('normal','dn'))) %>% 
  ggplot(aes(x = "", y = OxPhos, fill = disease)) +
  introdataviz::geom_flat_violin(trim=FALSE, alpha = 0.4,
                                 position = position_nudge(x = rain_height+.05)) +
  # rain
  #geom_point(aes(colour = disease), size = 2, alpha = .5, show.legend = FALSE, 
  #           position = position_jitter(width = rain_height * 0.5, height = 0)) +
  # boxplots
  #geom_boxplot(width = rain_height, alpha = 0.4, show.legend = FALSE, 
  #             outlier.shape = NA,
  #             position = position_nudge(x = -rain_height*2)) +
  # mean and SE point in the cloud
  stat_summary(fun.data = mean_cl_normal, mapping = aes(color = disease), show.legend = FALSE,
               position = position_nudge(x = rain_height * 2)) +
  stat_compare_means(group = 'disease', label = 'p.signif',
                     position = position_nudge(x = rain_height), label.y = 1)+
  coord_flip() +
  # custom colours and theme
  scale_fill_brewer(palette = "Dark2", name = "",
                    labels = c('Normal', 'DN')) +
  scale_x_discrete(name = "", expand = c(0, 0, 0, 0)) +
  scale_colour_brewer(palette = "Dark2") +
  ylab('OxPhos')+
  theme_pubclean()

ms <- tibble(umap_1 = KRM@reductions$umap@cell.embeddings[,1],
             umap_2 = KRM@reductions$umap@cell.embeddings[,2],
             disease_status = KRM$disease_status) %>% 
  mutate(disease = fct_relevel(disease_status, 'normal', 'dn'))

set.seed(123)

require(ggpubr)
ggplot(ms, aes(x = umap_1, y = umap_2)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = disease),
                  bins = 10) +
  theme_pubr(base_size = 12, legend = 'right') +
  ylim(c(-5, 5)) +
  guides(alpha = "none",
         fill = guide_legend(title= NULL)) +
  scale_fill_manual(values = c(gray(0.6, alpha = 1), 'firebrick4'),
                    labels = c('Normal', 'DN')) +
  theme(legend.position = c(.85, .9))

save(KRM, file = './raw_data/KRM.Rdata')



# wgcna
load(file = './raw_data/KRM.Rdata')

pacman::p_load(tidyverse, readxl, Seurat, data.table, ggsci, ggpubr,
               hdWGCNA, cowplot, patchwork, WGCNA, harmony)
theme_set(theme_cowplot()); set.seed(12345)

KRM <- subset(KRM, patient_name %in% names(table(KRM$patient_name))[table(KRM$patient_name)>=30])

KRM <- KRM %>% SetupForWGCNA(
  gene_select = "fraction", fraction = 0.07, wgcna_name = "KRM") %>% 
  MetacellsByGroups(group.by = 'patient_name', 
                    k = 10, min_cells = 30,  max_shared = 4, 
                    ident.group = 'patient_name') %>% 
  NormalizeMetacells() %>% 
  SetDatExpr(group_name = names(table(KRM$patient_name)),
             use_metacells = TRUE, group.by='patient_name') %>% 
  TestSoftPowers(networkType = 'signed', corFnc = 'cor')

plot_list <- PlotSoftPowers(KRM); wrap_plots(plot_list, ncol=2)

KRM <- ConstructNetwork(KRM, #soft_power,
                        minModuleSize = 60, deepSplit = 0,
                        setDatExpr=FALSE, tom_name = 'KRM', overwrite_tom = TRUE)

Module_num <- GetModules(KRM) %>% pull(module) %>% unique() %>% length() - 1

KRM <- KRM %>% 
  ModuleEigengenes(group.by.vars="patient_name") %>% 
  ModuleConnectivity() %>% 
  ResetModuleNames(new_name = "KRM")  %>% 
  ModuleExprScore(n_genes = 25, method='Seurat')

Col <- ggsci::pal_jco(alpha=0.8)(Module_num+1)
modules <- GetModules(KRM) 
modules$module_col <- Col[as.factor(modules$color) %>% as.integer()]
module_col <- modules %>% mutate(module_col = if_else(module == 'grey', 'grey', module_col)) %>% 
  pull(module_col)

KRM@misc$KRM$wgcna_modules$color <- module_col

PlotDendrogram(KRM, main='Dendrogram')
PlotKMEs(KRM, text_size = 3)

hub_df <- GetHubGenes(KRM, n_hubs = 10);hub_df
hub_df %>% arrange(module, -kME) %>% 
  select(-kME) %>% mutate(order = rep(1:10, Module_num)) %>% spread(key = module, value = gene_name) %>% 
  flextable::flextable()

wrap_plots(ModuleFeaturePlot(KRM, features='hMEs', order=TRUE))

MEs <- GetMEs(KRM, harmonized=TRUE); 
mods <- colnames(MEs); mods <- mods[mods != 'grey']
KRM@meta.data <- cbind(KRM@meta.data, MEs)

DotPlot(KRM, features=mods, group.by = 'disease_status', scale = T)

# GSEA
pacman::p_load('enrichR')
setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()
dbs <- c('MSigDB_Hallmark_2020')

gsea1 <- list()

hub_df_100 <- GetHubGenes(KRM, n_hubs = 100) %>% arrange(module, -kME)
# hub_df_100 %>% group_split(module) %>% writexl::write_xlsx('hub_df_230210.xlsx')

for (i in 1:Module_num){
  hub_df_i <- hub_df_100 %>% filter(module == str_c('KRM', i))
  enriched <- enrichr(hub_df_i$gene_name, dbs)
  
  gsea <- bind_rows(enriched)  %>% mutate(Term = str_replace(Term, ' \\(GO.*', ''))
  
  gsea <- gsea %>% select(Term, Adjusted.P.value) %>% filter(Adjusted.P.value < 0.05) %>% 
    mutate(module = str_c('KRM ', i))
  gsea1[[i]] <- gsea
}

gsea1 <- bind_rows(gsea1) 
gsea1 <- gsea1 %>% mutate(number = factor(letters[1:nrow(gsea1)]) , value = -log10(Adjusted.P.value)) 
ggplot(gsea1, aes(x = number, y = value)) +
  theme_pubr() +
  geom_col(aes(color = module), fill = NA, size = .75) +
  scale_x_discrete(limits = letters[nrow(gsea1):1], labels = rev(gsea1$Term)) + xlab('') + ylab('-log(p value)')+
  coord_flip() +
  scale_fill_manual(values = 'white', labels = '')+
  scale_color_manual(values = Col[c(1,2,4,5)]) +
  theme(legend.title = element_blank())

hub_2 <-  hub_df_100 %>% filter(module == 'KRM2')
enriched_2 <- enrichr(hub_2$gene_name, dbs)

compare_KRM <- KRM@meta.data %>% select(disease_status, starts_with('KRM')) 
compare_KRM <- compare_KRM[,colnames(compare_KRM) %>% sort()] 

compare_KRM %>% 
  pivot_longer(cols = 2:(Module_num+1)) %>%
  mutate(disease_status = fct_relevel(disease_status, 'normal','dn')) %>% 
  ggplot(aes(x = name, y = value, fill = disease_status)) +
  introdataviz::geom_split_violin(alpha = 1, trim = FALSE, width = 2) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
  scale_fill_brewer(palette = "Dark2", name = "",
                    labels = c('Normal', 'DN')) + xlab('') + ylab('Module Score')+
  stat_compare_means(aes(group = disease_status), label = 'p.signif',
                     label.y = 27) + 
  theme_pubr()

save(KRM, Module_num,  file = './raw_data/wgcna/wgcna_v1.RData')
load(file = './raw_data/wgcna/wgcna_v1.RData')

# 230208 nichenet
Idents_All <- as.character(Idents(dn_all1))
Imm_Order <- match(names(Idents(dn_immune1)), names(Idents(dn_all1)))
Idents_All[Imm_Order] <- as.character(Idents(dn_immune1))


dn_all2 <- dn_all1
dn_all2$cell_type <- Idents_All

dn_all2 <- subset(dn_all2, cell_type %notin% c('IMM-1','IMM-2','IMM-3','IMM-4','IMM-5','IMM-6','IMM-7','IMM-8'))
Idents(dn_all2) <- dn_all2$cell_type

save(dn_all2, file = './raw_data/wgcna/dn_nichenetr.RData')

pacman::p_load(nichenetr)
load(file = './raw_data/wgcna/dn_nichenetr.RData')

ligand_target_matrix = readRDS('./ext_source/ligand_target_matrix.rds')
lr_network = readRDS('./ext_source/lr_network.rds')
weighted_networks = readRDS('./ext_source/weighted_networks.rds')
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))


# Receiver KRM -- sender other
sender_celltypes <- c('PT','DL/ tAL','TAL','DCT/ CD-P','CD-I','SMC/ PERI','EC', 'PODO', 
                      'Infiltrating Mac', 'Monocyte','cDC','Neutrophil','CD4+ T', 
                      'CD8+ T','CD8+ T, effector','NK','B')

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = dn_all2, 
  receiver = "KRM", 
  condition_colname = "disease_status", condition_oi = "dn", condition_reference = "normal", 
  sender = sender_celltypes, 
  geneset = 'up',
  top_n_ligands = 14,
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "human")

DotPlot(dn_all2, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
DotPlot(dn_all2, features = nichenet_output$top_ligands %>% rev(), split.by = "disease_status") + RotatedAxis()

nichenet_output$ligand_differential_expression_heatmap
nichenet_output$ligand_target_heatmap

nichenet_output$top_targets %in% OxPhos
nichenet_output$top_receptors

?nichenet_seuratobj_aggregate

# Receiver other -- sender KRM

# scenic
pacman::p_load(Seurat, SCopeLoomR, SCENIC)
# 
# exprMat <- as.matrix(GetAssayData(KRM, slot = 'count'))
# exprMat <-  exprMat[rowSums(exprMat) >= 3,]
# SCopeLoomR::build_loom('wgcna_v1.loom', dgem = exprMat)
load(file = './raw_data/wgcna/wgcna_v1.RData')
pyscenic_wgcna <- SCopeLoomR::open_loom('./pyscenic/wgcna_v1_output.loom')

regulons_incidMat <- get_regulons(pyscenic_wgcna, column.attr.name = 'Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(pyscenic_wgcna, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(pyscenic_wgcna)
embeddings <- get_embeddings(pyscenic_wgcna)



regulonActivity_byCellType <- sapply(split(colnames(KRM), KRM$disease_status),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

Top_10_regulon <- regulonActivity_byCellType %>% as_tibble(rownames = 'gene') %>% 
  mutate(relative = dn / normal) %>% 
  arrange(-relative) %>% 
  slice(c(1:20, (n()-20):n())) %>% 
  mutate(gene = str_replace(gene, '\\(\\+\\)',''))

col_fun = circlize::colorRamp2(c(0, 1, 3), c("yellow", "white", "green"))

ComplexHeatmap::Heatmap(Top_10_regulon %>% column_to_rownames('gene') %>% select(relative) %>% 
                          as.matrix(),
                        col = col_fun,
                        cluster_rows = FALSE,
                        show_column_names = FALSE,
                        heatmap_legend_param = list(title = ''))

?Heatmap

?ComplexHeatmap::Heatmap

?make_heatmap_ggplot



Top_10_regulon

Top_10_regulon %>% flextable::flextable() %>% flextable::colformat_double(digits = 5)

ggplot(Top_10_regulon %>% pivot_longer(2:3), aes(x = name, y = gene)) +
  geom_point(aes(size = value, color = ggsci::pal_aaas()(1))) +
  scale_y_discrete(limits = rev(Top_10_regulon$gene)) +
  scale_radius(range = c(1,5)) +
  ggpubr::theme_pubr() +
  NoLegend() + xlab('') + ylab('')


make_heatmap_ggplot

#################################################################################
# pathway_analysis_inferring regulator
pacman::p_load(CIE, rcytoscapejs, clusterProfiler, ReactomePA)
?runCIE
KRM_dn_DEG
CIE::runCIE()

asdf <- KRM_dn_DEG %>% bind_cols(select(org.Hs.eg.db, keys = rownames(KRM_dn_DEG), columns = 'ENTREZID', keytype = 'SYMBOL')) %>% 
  dplyr::select('entrez' = ENTREZID, 'pval' = p_val_adj, 'fc' = avg_log2FC) %>% 
  filter(!is.na(entrez)) 
geneList <- 2^asdf$fc
names(geneList) <- asdf$entrez
geneList <- sort(geneList, decreasing  = T)
y <- gsePathway(geneList, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)
head(y)

viewPathway("Cristae formation", 
            readable = TRUE, 
            foldChange = geneList)

viewPathway("Mitochondrial biogenesis", 
            readable = TRUE, 
            foldChange = geneList)

viewPathway("Formation of ATP by chemiosmotic coupling", 
            readable = TRUE, 
            foldChange = geneList)

viewPathway("Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.", 
            readable = TRUE, 
            foldChange = geneList)

viewPathway("The citric acid (TCA) cycle and respiratory electron transport", 
            readable = TRUE, 
            foldChange = geneList)

# trajectory analysis
pacman::p_load(monocle3)

myeloid_cells <- dn_immune1 %>% subset(idents = c('cDC', 'KRM', 'Monocyte', 'Neutrophil', 'Infiltrating Mac'))
myeloid_cells <- myeloid_cells %>% RunUMAP(reduction = 'harmony', dims = 1:30)

myeloid_cds <- new_cell_data_set(
  myeloid_cells@assays$RNA@data,
  cell_metadata  = myeloid_cells@meta.data
)

myeloid_cds <- preprocess_cds(myeloid_cds, num_dim = 50)
myeloid_cds <- align_cds(myeloid_cds, alignment_group = "channel")
myeloid_cds <- reduce_dimension(myeloid_cds)

myeloid_cds@reduce_dim_aux$UMAP$model$umap_model$embedding <- myeloid_cells@reductions$umap@cell.embeddings
myeloid_cds@int_colData@listData$reducedDims@listData$UMAP <- myeloid_cells@reductions$umap@cell.embeddings

plot_cells(myeloid_cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell_type", reduction_method = 'UMAP')

myeloid_cds <- cluster_cells(myeloid_cds)
plot_cells(myeloid_cds, color_cells_by = "partition")
myeloid_cds <- learn_graph(myeloid_cds)
plot_cells(myeloid_cds,
           color_cells_by = "cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size = 4,
           group_label_size = 4)
myeloid_cds <- order_cells(myeloid_cds)
plot_cells(myeloid_cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

?plot_cells

save(myeloid_cds, myeloid_cells, file = './raw_data/monocle_v1.Rdata')
load(file = './raw_data/monocle_v1.Rdata')
