rm(list=ls())
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

`%notin%` <- Negate(`%in%`)

dn_all <- dn_all %>% subset(is.na(patient_name) | patient_name %notin% c('Wilms1', 'Wilms2','Wilms3','Trans'))
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
  '6' = 'DL, tAL',
  '3' = 'TAL',
  '5' = 'DCT, CD-P',
  '8' = 'CD-I',
  '12' = 'SMC, PERI',
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

#########################################################################


load(file = './raw_data/wgcna/dn_immune_v2.RData')

dn_immune$percent.mt <- PercentageFeatureSet(object = dn_immune, pattern = "^MT-")
dn_immune <- dn_immune %>% subset(percent.mt <= 10) %>% clustering(resol = .6) 

dn_immune <- dn_immune %>% FindClusters(resolution=1.1)

mark <- c('HLA-DRA',  'LYZ', 'CD14', 'CD68',   'C1QC', 'MRC1','LYVE1',
          'FCN1','PLAC8',  'THBS1', 'VCAN',  'FCER1A', 'CD1C', 'S100A8','IL1R2',  
          'CD3E', 'TRAC', 'CD4','CD8A', 'GZMB', 'IFNG', 'IGKC', 'JCHAIN', 'CD79A','NKG7', 'KLRD1',
          'CD40LG')

quick(dn_immune)
quickdot(dn_immune, feat = mark)


dn_immune1 <- dn_immune %>% RenameIdents(
  '11' = 'KRM',
  '7' = 'Mac/Mono',
  '12' = 'Neutrophil',
  '5' = 'cDC',
  '3' = 'CD4+ T',
  '2' = 'CD4+ T',
  '15' = 'CD4+ T',
  '17' = 'CD8+ T',
  '1' = 'CD8+ T, effector',
  '3' = 'NK',
  '4' = 'NK',
  '10' = 'NK(?)',
  '7' = 'B',
  '8' = 'B',
  '9' = 'B',
  '13' = 'B'
)

quick(dn_immune1)

# DEG_ident(dn_immune1)

non_immune <- c('PT','DL, tAL','TAL','DCT, CD-P','CD-I','SMC, PERI','EC-1', 'EC-2', 'PODO', 20)
immune <- c('KRM', 'Infiltrating Mac', 'Monocyte','cDC','Neutrophil','CD4+ T',
            'CD8+ T','CD8+ T, effector','NK','B')

DEG_ident(dn_all1)
DEG_ident(dn_immune)

DEG_dn <- list()
for (i in 0:18){
  subset_seurat <- dn_immune %>% subset(idents = i)
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

asdf <- FindMarkers(dn_all1, ident.1 = 'EC-1',ident.2 = 'EC-2')

writexl::write_xlsx(DEG_dn, 'DEG_dn.xlsx')

################################# clustering
KRM <- subset(dn_immune1, idents = 'KRM') %>% clustering(resol = 0.3)

Idents(KRM) <- KRM$disease_status
KRM_dn_DEG <- FindMarkers(KRM, ident.1 = 'dn')
KRM_dn_DEG %>% arrange(desc(avg_log2FC))
KRM_dn_DEG$p_val_adj
EnhancedVolcano::EnhancedVolcano(KRM_dn_DEG, x = 'avg_log2FC',
                                 y = 'p_val_adj', 
                                 lab = row.names(KRM_dn_DEG),
                                 pCutoff = 10e-10,
                                 col=c('black', 'black', 'black', 'red3'),
                                 drawConnectors = TRUE,
                                 title = NULL,
                                 subtitle = NULL)

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

OxPos <- gsea_H$leadingEdge[[1]]
KRM <- AddModuleScore(
  object = KRM,
  features = OxPos,
  nbin = 12,
  ctrl = 100,
  name = 'OxPos'
)
require(ggpubr)
tibble(disease = KRM$disease_status, OxPos = KRM$OxPos1) %>% 
  ggboxplot(x = 'disease',y = 'OxPos', color = 'disease') +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = .8), aes(color = disease))+
  xlab('') + 
  theme_pubr()  +
  stat_compare_means(group = 'disease', label = 'p.signif') +
  guides(fill = guide_legend(title = NULL)) 


quick(KRM)

ms <- tibble(umap_1 = KRM@reductions$umap@cell.embeddings[,1],
             umap_2 = KRM@reductions$umap@cell.embeddings[,2],
             disease_status = KRM$disease_status)

set.seed(123)
require(ggpubr)
ggplot(ms, aes(x = umap_1, y = umap_2)) +
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level.., fill = disease_status),
                  bins = 7) +
  theme_pubr(base_size = 12, legend = 'right') +
  ylim(c(-5, 5)) +
  guides(alpha = "none",
         fill = guide_legend(title= NULL))

# wgcna
pacman::p_load(tidyverse, readxl, Seurat, data.table, ggsci, ggpubr,
               hdWGCNA, cowplot, patchwork, WGCNA, harmony)
theme_set(theme_cowplot()); set.seed(12345)

KRM <- subset(KRM, patient_name %in% names(table(KRM$patient_name))[table(KRM$patient_name)>=20])

KRM <- KRM %>% SetupForWGCNA(
  gene_select = "fraction", fraction = 0.05, wgcna_name = "KRM") %>% 
  MetacellsByGroups(group.by = 'patient_name', 
                    k = 7, min_cells = 20,  max_shared = 5, 
                    ident.group = 'patient_name') %>% 
  NormalizeMetacells() %>% 
  SetDatExpr(group_name = names(table(KRM$patient_name)),
             use_metacells = TRUE, group.by='patient_name') %>% 
  TestSoftPowers(networkType = 'signed', corFnc = 'cor')

plot_list <- PlotSoftPowers(KRM); wrap_plots(plot_list, ncol=2)

KRM <- ConstructNetwork(KRM, #soft_power,
                        minModuleSize = 50, deepSplit = 0,
                        setDatExpr=FALSE, tom_name = 'KRM', overwrite_tom = TRUE)

PlotDendrogram(KRM, main='Dendrogram')

KRM <- NormalizeData(KRM) %>%  
  FindVariableFeatures(selection.method = "vst", nFeatures = num_HVG) %>% 
  ScaleData() %>% 
  ModuleEigengenes(group.by.vars="patient_name")
KRM <- KRM %>% ModuleConnectivity() %>% 
  ResetModuleNames(new_name = "KRM")  
KRM <- KRM %>% ModuleExprScore(n_genes = 25, method='Seurat')



modules <- GetModules(KRM); head(modules[,1:6])
hub_df <- GetHubGenes(KRM, n_hubs = 10);hub_df
hub_df %>% arrange(module, -kME) %>% 
  select(-kME) %>% mutate(order = rep(1:10, 7)) %>% spread(key = module, value = gene_name) %>% 
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
for (i in 1:7){
  hub_df_i <- hub_df %>% filter(module == str_c('KRM', i))
  enriched <- enrichr(hub_df_i$gene_name, dbs)
  
  gsea <- bind_rows(enriched)  %>% mutate(Term = str_replace(Term, ' \\(GO.*', ''))
  
  gsea <- gsea %>% select(Term, Adjusted.P.value) %>% filter(Adjusted.P.value < 0.05) %>% 
    mutate(module = str_c('KRM ', i))
  gsea1[[i]] <- gsea
}
gsea1 <- bind_rows(gsea1) %>% mutate(number = letters[1:11] , value = -log10(Adjusted.P.value)) 

ggplot(gsea1, aes(x = number, y = value)) +
  theme_pubr() +
  geom_col(aes(fill = module)) +
  scale_x_discrete(limits = letters[11:1], labels = rev(gsea1$Term)) + xlab('') + ylab('-log(p value)')+
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 3) +theme(legend.title = element_blank())


tibble(disease = KRM$disease_status, 
       KRM1 = KRM$KRM1,
       KRM2 = KRM$KRM2,
       KRM3 = KRM$KRM3,
       KRM4 = KRM$KRM4,
       KRM5 = KRM$KRM5,
       KRM6 = KRM$KRM6,
       KRM7 = KRM$KRM7) %>% 
  pivot_longer(cols = 2:8) %>% 
  ggboxplot(x = 'name',y = 'value', color = 'disease') +
  geom_jitter(position=position_jitterdodge(jitter.width = .1, dodge.width = .8), aes(color = disease))+
  xlab('') + 
  theme_pubr()  +
  stat_compare_means(aes(group = disease), label = 'p.signif') +
  guides(fill = guide_legend(title = NULL)) 

save(KRM, dn_immune1, file = './raw_data/wgcna/wgcna_v1.RData')
