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
          'CD3E', 'TRAC', 'CD4', 'CD40LG','CD8A', 'GZMB', 'IFNG', 'IGKC', 'JCHAIN', 'CD79A','NKG7', 'KLRD1'
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

#DEG_ident(dn_all1)
#DEG_ident(dn_immune1)

non_immune <- c('PT','DL, tAL','TAL','DCT, CD-P','CD-I','SMC, PERI','EC', 'PODO')
immune <- c('KRM', 'Infiltrating Mac', 'Monocyte','Neutrophil', 'cDC',
            'CD4+ T', 'CD8+ T','CD8+ T, effector','NK','B')


DEG_dn <- list()
for (i in immune){
  subset_seurat <- dn_immune1 %>% subset(idents = i)
  Idents(subset_seurat) <- subset_seurat$disease_status
  dn_deg <- FindMarkers(subset_seurat, ident.1 = 'dn', only.pos = TRUE)
  DEG_dn[[i]] <- dn_deg %>% rownames_to_column('Gene_name')
}
for (i in non_immune){
  subset_seurat <- dn_all1 %>% subset(idents = i)
  Idents(subset_seurat) <- subset_seurat$disease_status
  dn_deg <- FindMarkers(subset_seurat, ident.1 = 'dn')
  DEG_dn[[i]] <- dn_deg %>% rownames_to_column('Gene_name')
}

writexl::write_xlsx(DEG_dn, 'DEG_dn.xlsx')

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
  ggpubr::theme_pubclean()+
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
  
  gsea_H %>% write_xlsx('gsea_hallmark_230213.xlsx')
  
  p <- ggplot(gsea_H, aes(reorder(pathway_new, NES), NES)) +
    geom_col(aes(fill=significant)) +
    scale_fill_manual(values = c(pal_nejm(alpha = 0.8)(1),pal_nejm(alpha = 0.9)(1)),
                      breaks = c('p Value < 0.05')) +
    coord_flip() +
    labs(title="Hallmark pathways NES from GSEA") + 
    guides(fill=guide_legend(title=NULL)) +
    xlab('')+
    ylab('Normalized enrichment score')+
    ggpubr::theme_pubclean()+
    # theme_classic2()+
    theme(legend.position = c(0.85, 0.5))
  return(p)
}
Infilt <- S4Vectors::subset(dn_immune1, idents = 'Infiltrating Mac')
subset_gsea(Infilt)
Mono <- S4Vectors::subset(dn_immune1, idents = 'Monocyte')
subset_gsea(Mono)

Idents(dn_immune1)

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

tibble(disease = KRM$disease_status, OxPhos = KRM$OxPhos1) %>%
  mutate(disease = fct_relevel(disease, c('normal','dn'))) %>% 
  ggviolin(x = 'disease',y = 'OxPhos', fill = 'disease', add = 'mean_sd') +
  xlab('') + 
  theme_pubr()  +
  stat_compare_means(group = 'disease', label = 'p.signif',
                     label.x = 1.5, label.y = 5) +
  theme(legend.title = element_blank()) +
  scale_fill_brewer(palette = "Dark2", labels = c('Normal', 'DN'))+
  scale_x_discrete(labels = c('', '')) +
  ylab('OxPhos')
  

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

?scale_x_discrete
quick(KRM)

ms <- tibble(umap_1 = KRM@reductions$umap@cell.embeddings[,1],
             umap_2 = KRM@reductions$umap@cell.embeddings[,2],
             disease_status = KRM$disease_status) %>% 
  mutate(disease = fct_relevel(disease_status, 'normal', 'dn'))

set.seed(123)

require(ggpubr)
ggplot(ms, aes(x = umap_1, y = umap_2)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = disease),
                  bins = 7) +
  theme_pubr(base_size = 12, legend = 'right') +
  ylim(c(-5, 5)) +
  guides(alpha = "none",
         fill = guide_legend(title= NULL)) +
  scale_fill_manual(values = c(gray(0.6, alpha = 1), 'firebrick1'),
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
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "human")


DotPlot(dn_all2, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
DotPlot(dn_all2, features = nichenet_output$top_ligands %>% rev(), split.by = "disease_status") + RotatedAxis()

nichenet_output$ligand_differential_expression_heatmap
nichenet_output$ligand_target_heatmap

nichenet_output$top_targets %in% OxPhos
nichenet_output$top_receptors

# Receiver other -- sender KRM
sender = 'KRM'

expressed_genes_sender = get_expressed_genes(sender, dn_all2, pct = 0.10)
background_expressed_genes = expressed_genes_sender %>% .[. %in% colnames(ligand_target_matrix)]

seurat_obj_sender = subset(dn_all2, idents = sender)
seurat_obj_sender = SetIdent(seurat_obj_sender, value = seurat_obj_sender[["disease_status"]])

condition_oi = "dn"
condition_reference = "normal" 

DE_table_sender = FindMarkers(object = seurat_obj_sender, ident.1 = condition_oi, 
                                ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

ligand_oi = DE_table_sender %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
ligand_oi = ligand_oi %>% .[. %in% colnames(ligand_target_matrix)]


## receiver
receiver_celltypes <- c('PT','DL/ tAL','TAL','DCT/ CD-P','CD-I','SMC/ PERI','EC', 'PODO', 
                      'Infiltrating Mac', 'Monocyte','cDC','Neutrophil','CD4+ T', 
                      'CD8+ T','CD8+ T, effector','NK','B')
list_expressed_genes_receiver = receiver_celltypes %>% unique() %>% lapply(get_expressed_genes, dn_all2, 0.10)
expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()

# list_DEGs_receiver <- list()
# for (i in receiver_celltypes){
#   list_DEGs_receiver[[i]] <- FindMarkers(dn_all2, ident.1 = i, only.pos = TRUE)
# }
# DEGs_receiver <- list_DEGs_receiver %>% bind_rows()
# saveRDS(DEGs_receiver, file = './ext_source/DEGs_receiver.rds')
DEGs_receiver <- readRDS(file = './ext_source/DEGs_receiver.rds')
geneset_oi = DEGs_receiver %>% rownames_to_column('gene') %>% 
  filter(p_val_adj <= 0.05 & avg_log2FC >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#3
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#4
ligand_activities = ligand_oi[ligand_oi %in% potential_ligands]


#5
active_ligand_target_links_df = ligand_activities %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(ligand_activities, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands in other cells",
                                                                    "Predicted target genes in KRM", color = "purple",
                                                                    legend_position = "top", x_axis_position = "top",
                                                                    legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network + NoLegend()

#6
DE_table_all = Idents(dn_all2) %>% levels() %>% intersect(sender_celltypes) %>% 
  lapply(get_lfc_celltype, seurat_obj = dn_all2, condition_colname = "disease_status", 
         condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% 
  reduce(full_join)
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
