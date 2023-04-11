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


mito <- c('COX4I1', 'COX4I2', 'COX5A', 'NDUFS1', 'NDUFS2', 'MT-ND3', 'SDHA')

quickdot(dn_all1, feat = mito)

# prepare for 14, 22 sorting (230210)

Idents_All <- as.character(Idents(dn_all1))
Imm_Order <- match(names(Idents(dn_immune)), names(Idents(dn_all1)))
Idents_All[Imm_Order] <- as.character(Idents(dn_immune1))

dn_all2 <- dn_all1
dn_all2$cell_type <- Idents_All

dn_all2 <- subset(dn_all2, cell_type %notin% c('IMM-1','IMM-2','IMM-3','IMM-4','IMM-5','IMM-6','IMM-7','IMM-8'))
Idents(dn_all2) <- dn_all2$cell_type

quickdot(dn_all2, feat = mito)

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

# a <- list('krm_mono' = (FindMarkers(dn_immune1, ident.1 = 'KRM', ident.2 = 'Monocyte') %>% rownames_to_column('gene')),
#      'krm_infilt' = (FindMarkers(dn_immune1, ident.1 = 'KRM', ident.2 = 'Infiltrating Mac') %>% rownames_to_column('gene'))
# )
# 
# writexl::write_xlsx(a, 'krm_deg.xlsx')

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
# KRM <- S4Vectors::subset(dn_immune1, idents = 'KRM')
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

tlr <- c('TLR1', 'TLR2', 'TLR3', 'TLR4',
         'TNFRSF1A', 'TNFRSF1B', 'LTBR', 'CD40', 'FAS', 'TNFRSF4', 'CD27')
tlr1 <- c('TLR1', 'TLR2', 'TLR3','TLR4','TLR5','TLR6','TLR7','TLR8','TLR9', 'TLR10')

FoldChange(KRM, ident.1 = 'dn', features = tlr1) %>% select(1) %>% 
  `colnames<-`('KRM') %>% 
  ComplexHeatmap::Heatmap(name = ' ',
                          cluster_rows = FALSE, cluster_columns = FALSE,   
                          show_column_names = TRUE,
                          column_names_rot = 30)
quickdot(dn_immune1, feat = rev(tlr))
quickdot(dn_immune1, feat = rev(tlr1))

quickdot(dn_immune1, feat = 'IL10RA')


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
KRM$disease_status %>% table()
KRM
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

Col0 <- ggsci::pal_jco(alpha=0.8)(Module_num+1)
Col <- Col0[c(1,2,4,3,5,6)]
modules <- GetModules(KRM) 
modules$module_col <- Col[as.factor(modules$color) %>% as.integer()]
module_col <- modules %>% mutate(module_col = if_else(module == 'grey', 'grey', module_col)) %>% 
  pull(module_col)

KRM@misc$KRM$wgcna_modules$color <- module_col



PlotDendrogram(KRM, main='Dendrogram')

ggplot(tibble(clarity = 1, cut = letters[1:6]), aes(clarity, fill = cut)) +
  geom_bar() +
  scale_fill_manual(values = c(Col[c(5,2,6,3,1)], 'grey'), labels = c('KRM1', 'KRM2','KRM3','KRM4','KRM5', 'Not defined'))


PlotKMEs(KRM, text_size = 3)
ModuleNetworkPlot(KRM)

hub_df <- GetHubGenes(KRM, n_hubs = 10);hub_df
hub_df %>% arrange(module, -kME) %>% 
  select(-kME) %>% mutate(order = rep(1:10, Module_num)) %>% spread(key = module, value = gene_name) %>% 
  flextable::flextable()

wrap_plots(ModuleFeaturePlot(KRM, features='hMEs', order=TRUE))

MEs <- GetMEs(KRM, harmonized=TRUE); 
mods <- colnames(MEs); mods <- mods[mods != 'grey']
KRM@meta.data <- cbind(KRM@meta.data, MEs)

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
  
  gsea <- gsea %>% select(Term, Adjusted.P.value) %>% 
    filter(Adjusted.P.value < 0.05) %>% 
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

for (i in 1:Module_num){
  hub_df_i <- hub_df_100 %>% filter(module == str_c('KRM', i))
  enriched <- enrichr(hub_df_i$gene_name, dbs)
  
  gsea <- bind_rows(enriched)  %>% mutate(Term = str_replace(Term, ' \\(GO.*', ''))
  
  gsea <- gsea %>% select(Term, Adjusted.P.value) %>% 
    #filter(Adjusted.P.value < 0.05) %>% 
    mutate(module = str_c('KRM ', i))
  gsea1[[i]] <- gsea
}
gsea1 <- bind_rows(gsea1) 
gsea1 <- gsea1 %>% mutate(number = factor(letters[1:nrow(gsea1)]) , value = -log10(Adjusted.P.value)) 

gsea2 <- gsea1 %>% filter(str_detect(Term, 'xidat')) 
gsea2 %>% 
  ggplot(aes(x = module, y = value)) +
  geom_col(aes(color = module), fill = NA, linewidth = .75) + xlab('') +
  ylab('OXPHOS adjusted p-value')+
  coord_flip() +
  theme_pubr() +
  theme(legend.title = element_blank())

hub_2 <-  hub_df_100 %>% filter(module == 'KRM2')
enriched_2 <- enrichr(hub_2$gene_name, dbs)
hub_5 <-  hub_df_100 %>% filter(module == 'KRM5')
enriched_5 <- enrichr(hub_5$gene_name, dbs)

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



# DMEs
pacman::p_load(ggrepel, igraph)
group1 <- KRM@meta.data %>% subset(disease_status == 'dn') %>% rownames
group2 <- KRM@meta.data %>% subset(disease_status == 'normal') %>% rownames

PlotDMEsVolcano(KRM, FindDMEs(KRM, group1, group2))

HubGeneNetworkPlot(
  KRM,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)

c1 <- c('NDUFA1', 'NDUFA2', 'NDUFA3', 'NDUFB1', 'NDUFB2')
c2 <- c('SDHA', 'SDHB', 'SDHC')
c3 <- c('UQCRB', 'UQCRQ', 'CYC1')
c4 <- c('COX5A', 'COX5B', 'COX6C', 'COX8A', 'COX7B')
c5 <- c('ATP5PD','ATP5PF','ATP5MG', 'ATP5MC2','ATP5MC3', 'ATP5ME', 'ATP5F1E')

library(circlize)
col_fun = colorRamp2(c(0, 1), c("white", "red"))

FoldChange(KRM, ident.1 = 'dn', features = c(c1, c2, c3, c4, c5)) %>% select(1) %>% 
  `colnames<-`('KRM') %>% 
  ComplexHeatmap::Heatmap(name = ' ',col = col_fun,
                          cluster_rows = FALSE, cluster_columns = FALSE,   
                          show_column_names = TRUE,
                          column_names_rot = 30)


save(KRM, Module_num,  file = './raw_data/wgcna/wgcna_v1.RData')
load(file = './raw_data/wgcna/wgcna_v1.RData')

# Network Plot
library(igraph)
ModuleNetworkPlot(KRM, mods = 'KRM2')
?ModuleNetworkPlot
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

p <- nichenet_output$ligand_differential_expression_heatmap

Immune_cluster <- c('SMC/PERI', 'PT', 'PODO', 'DL/tAL', 'CD4+ T', 'DCT/CD-P', 'TAL', 'CD-I',
                    'CD8+ T, effector', 'NK', 'Neutrophil', 'EC', 'CD8+ T', 'B', 'Infiltrating Mac', 'Monocyte')
p + scale_x_discrete(labels = Immune_cluster, position = 'top')

?scale_x_discrete

nichenet_output$ligand_target_heatmap

nichenet_output$top_targets %in% OxPhos
nichenet_output$top_receptors

nichenet_output$ligand_activities
  
nichenet_output$ligand_activity_target_heatmap

ligand_activity <- nichenet_output$ligand_activities %>% pull(pearson)

ligand_activity <- ligand_activity[1:11] %>% as.matrix(ncol = 1) 

ligand_pearson_matrix = (nichenet_output$ligand_activities) %>% select(pearson) %>% as.matrix() %>% 
  magrittr::set_rownames(nichenet_output$ligand_activities$test_ligand)
vis_ligand_pearson = ligand_pearson_matrix[14:1, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>%
  make_heatmap_ggplot("Prioritized ligands","Ligand activity", 
                      color = "darkorange",legend_position = "top", x_axis_position = "top", 
                      legend_title = "Pearson correlation coefficient\ntarget gene prediction ability")

ligand_target_matrix <- (nichenet_output$ligand_target_matrix)[,c(4,5, 1:3, 6:39)]

ligand_target_matrix %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", 
                                             color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple") + theme(axis.text.x = element_text(face = "italic"))


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
  dplyr::slice(c(1:20, (n()-20):n())) %>% 
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
# trajectory analysis
pacman::p_load(monocle3)

myeloid_cells <- dn_immune1 %>% subset(idents = c('KRM', 'Monocyte', 'Infiltrating Mac'))
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
myeloid_cds <- order_cells(myeloid_cds)
plot_cells(myeloid_cds,
           color_cells_by = "cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size = 4,
           group_label_size = 4)
plot_cells(myeloid_cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

save(myeloid_cds, myeloid_cells, file = './raw_data/monocle_v1.Rdata')
load(file = './raw_data/monocle_v1.Rdata')


# cellchat
pacman::p_load(CellChat, ComplexHeatmap)


dn_immune1_dn <- subset(dn_immune1, disease_status == 'dn')
dn_immune1_normal <- subset(dn_immune1, disease_status == 'normal')

#dn

data_dn <- GetAssayData(dn_immune1_dn, assay = "RNA", slot = "data")
meta_dn <- data.frame(labels =  Idents(dn_immune1_dn), 
                       row.names = names(Idents(dn_immune1_dn))) 
cellchat_dn <- createCellChat(object = data_dn, meta = meta_dn, group.by = "labels")
cellchat_dn@DB <- CellChatDB.human

cellchat_dn <- cellchat_dn %>% subsetData() %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  projectData(PPI.human)
cellchat_dn <- cellchat_dn %>% computeCommunProb(raw.use = FALSE) %>%   
  filterCommunication(min.cells = 10) %>% 
  computeCommunProbPathway() %>%
  aggregateNet() %>% 
  netAnalysis_computeCentrality() %>% 
  identifyCommunicationPatterns(pattern = "outgoing", k = 3)

netAnalysis_signalingRole_heatmap(cellchat_dn, pattern = "outgoing",
                                  height = 18)
netAnalysis_dot(cellchat_dn, pattern = "outgoing")

# normal
data_normal <- GetAssayData(dn_immune1_normal, assay = "RNA", slot = "data")
meta_normal <- data.frame(labels =  Idents(dn_immune1_normal), 
                         row.names = names(Idents(dn_immune1_normal))) 
cellchat_normal <- createCellChat(object = data_normal, meta = meta_normal, group.by = "labels")
cellchat_normal@DB <- CellChatDB.human

cellchat_normal <- cellchat_normal %>% subsetData() %>% 
  identifyOverExpressedGenes() %>% 
  identifyOverExpressedInteractions() %>% 
  projectData(PPI.human)
cellchat_normal <- cellchat_normal %>% 
  computeCommunProb(raw.use = FALSE) %>%  
  filterCommunication(min.cells = 10) %>% 
  computeCommunProbPathway() %>%
  aggregateNet() %>% 
  netAnalysis_computeCentrality() %>% 
  identifyCommunicationPatterns(pattern = "outgoing",k = 3)


object.list <- list(old = cellchat_old, young = cellchat_young)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, sources.use = 'KRM', weight.scale = T, measure = "weight")
netVisual_diffInteraction(cellchat, targets.use = 'KRM', weight.scale = T, measure = "weight")
netVisual_heatmap(cellchat,measure = "weight" )
rankNet(cellchat, mode = "comparison", sources.use = 'KRM', targets.use = 'EC',
        stacked = T, do.stat = TRUE, measure = 'weight')


pacman::p_load(singlecellsignalR)
LRdb <- SingleCellSignalR::LRdb
require(SingleCellSignalR)

dn_immune1_dn$class <- as.integer(factor(Idents(dn_immune1_dn)))
dn_immune1_normal$class <- as.integer(factor(Idents(dn_immune1_normal)))

signal_dn <- SingleCellSignalR::cell_signaling(data = dn_immune1_dn@assays$RNA@data,
                             genes = rownames(dn_immune1_dn@assays$RNA@data),
                             cluster = dn_immune1_dn$class,
                             c.names = levels(factor(Idents(dn_immune1_dn))),
                             write = FALSE)
signal_normal <- SingleCellSignalR::cell_signaling(data = dn_immune1_normal@assays$RNA@data,
                               genes = rownames(dn_immune1_normal@assays$RNA@data),
                               cluster = dn_immune1_normal$class,
                               c.names = levels(factor(Idents(dn_immune1_normal))),
                               write = FALSE)

inter_dn <- inter_network(data = dn_immune1_dn@assays$RNA@data,
                           signal = signal_dn,
                           genes = rownames(dn_immune1_dn@assays$RNA@data),
                          c.names = levels(factor(Idents(dn_immune1_dn))),
                           cluster = dn_immune1_dn$class, write = FALSE)

inter_normal <- inter_network(data = dn_immune1_normal@assays$RNA@data,
                           signal = signal_normal,
                           genes = rownames(dn_immune1_normal@assays$RNA@data),
                           c.names = levels(factor(Idents(dn_immune1_normal))),
                           cluster = dn_immune1_normal$class, write = FALSE)

signal_dn

namechange <- function(sig){
  signal1 <- lapply(sig, function(x){
    x$cell_L <- colnames(x)[1]
    x$cell_R <- colnames(x)[2]
    colnames(x)[1] <- 'gene_L'
    colnames(x)[2] <- 'gene_R'
    return(x)})
  signal1 <- do.call(rbind, signal1)
  return(signal1)
}

interest_dn <- namechange(signal_dn) %>% filter(cell_R == 'KRM', cell_L %in% c('CD4+ T', 'CD8+ T', 'CD8+ T, effector'))



two_chord <- function(x, cellname){
  circos.par(canvas.ylim=c(-1.2,1.2), # edit  canvas size 
             #track.margin = c(0.01, 0.05), # adjust bottom and top margin
             track.height = 0.3)
  chordDiagram(chord(x),
               group = structure(c(x$cell_L, x$cell_R),
                                 names = c(x$gene_L, x$gene_R)),
               annotationTrack = 'none', col = col_7[1],
               directional = 1, direction.type = c("diffHeight"), diffHeight = -mm_h(2), 
               preAllocateTracks = list(track.height = mm_h(4)),
               link.arr.type = 'big.arrow')
  circos.track(ylim = c(0,1), track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] - mm_y(-6), CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 1))
  }, bg.border = NA)
  highlight.sector(x$gene_L,  
                   text = "KRM",col = col_7[2], 
                   cex = 0.8, text.col = "white", niceFacing = TRUE,
                   track.index = 1,  lwd = 10)
  highlight.sector(x$gene_R, 
                   text = cellname, col = col_7[3], 
                   cex = 0.8, text.col = "white", niceFacing = TRUE,
                   track.index = 1,  lwd = 10)
}

require(circlize)
col_7 <- ggsci::pal_nejm(alpha=0.6)(7)
chord <- function(x){tibble(from = x$gene_L, to = x$gene_R, value = x$LRscore)}


interest_dn1 <- interest_dn %>% filter(LRscore>=0.85)
two_chord(interest_dn, 'KRM')
interest_dn

chordDiagram(chord(interest_dn),
             group = structure(c(interest_dn$cell_L, interest_dn$cell_R),
                               names = c(interest_dn$gene_L, interest_dn$gene_R)),
             annotationTrack = 'none', col = col_7[1],
             directional = 1, direction.type = c("diffHeight"), diffHeight = -mm_h(2), 
             preAllocateTracks = list(track.height = mm_h(4)),
             link.arr.type = 'big.arrow')

visualize_interactions(signal_dn, show.in = c(50))
visualize_interactions(signal_dn, show.in = c(60))
visualize_interactions(signal_dn, show.in = c(50, 60, 70))

visualize_interactions(signal_normal, show.in = c(53))
visualize_interactions(signal_normal, show.in = 60)
visualize_interactions(signal_normal, show.in = c(53, 60, 71))

circos.clear()
chordDiagram(chord(interest_dn1), annotationTrack = 'none')
circos.track(ylim = c(0,1), track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] - mm_y(-6), CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 1))
}, bg.border = NA)

highlight.sector(interest_dn1 %>% filter(cell_L == 'CD4+ T') %>% pull(gene_L),  
                 text = "CD4+ T",col = col_7[2], 
                 cex = 0.8, text.col = "white", niceFacing = TRUE,
                 track.index = 1,  lwd = 10)
highlight.sector(interest_dn1 %>% filter(cell_L == 'CD8+ T') %>% pull(gene_L),  
                 text = "CD8+ T",col = col_7[2], 
                 cex = 0.8, text.col = "white", niceFacing = TRUE,
                 track.index = 1,  lwd = 10)
highlight.sector(interest_dn1 %>% filter(cell_L == 'CD8+ T, effector') %>% pull(gene_L),  
                 text = "CD8+ T, effector",col = col_7[2], 
                 cex = 0.8, text.col = "white", niceFacing = TRUE,
                 track.index = 1,  lwd = 10)
highlight.sector(interest_dn1$gene_R, 
                 text = 'KRM', col = col_7[3], 
                 cex = 0.8, text.col = "white", niceFacing = TRUE,
                 track.index = 1,  lwd = 10)


interest_dn1 %>% transmute(from = str_c(cell_L, ' ' ,gene_L), to = gene_R, value = LRscore) %>% 
  chordDiagram(annotationTrack = 'none')
circos.info()
