##########TRANSCRIPTIONAL FACTORS ACTIVITY INFERENCE######################
#BiocManager::install('decoupler')
library('decoupleR')
library(tidyverse)
library(cowplot)
#remotes::install_github("cvarrichio/Matrix.utils")
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)

Idents(Tubo_tumor_annotated)= Tubo_tumor_annotated$neigb_type
FCan_Data<- subset(Tubo_tumor_annotated, idents=c('Fibroblast/Stromal cells', 'Cancer cells', 'Cancer cells_Fibroblast/Stromal cells'))
FCan_CR<-subset(FCan_Data, subset= survival=='CR')
table(FCan_CR$neigb_class)
FCan_PD<-subset(FCan_Data, subset= survival=='PD')
table(FCan_PD$neigb_class)
FCan_PR<- subset(FCan_Data, subset=survival=='PR')
table(FCan_PR$neigb_class)

######pseudobulk approach#####################################################
Idents(Tubo_tumor_annotated)= Tubo_tumor_annotated$neigb_type

Olb_Data<- subset(Tubo_tumor_annotated, idents=c('Fibroblast/Stromal cells', 'Cancer cells', 
                                                 'Cancer cells_Fibroblast/Stromal cells', 'Fibroblast/Stromal cells_Cancer cells' ))
Olb_Data$neigb_type<- str_replace(Olb_Data$neigb_type,
                                  'Fibroblast/Stromal cells_Cancer cells',
                                  'Cancer cells_Fibroblast/Stromal cells')
Idents(Olb_Data) <- Olb_Data$neigb_type
cluster_names <- levels(Idents(Olb_Data)) 
groups<- Olb_Data@meta.data[, c("neigb_type", "survival")]
head(Olb_Data@assays$RNA@counts)
head(groups)
# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(Olb_Data@assays$RNA@data), 
                                groupings = groups, fun = "sum") 
agg_data<- as.data.frame(t(aggr_counts)) #%>% head()
acts <- run_ulm(mat=agg_data, net=collectri, .source='source', .target='target',
                .mor='mor', minsize = 5)

saveRDS(acts, 'acts.rds')

###viewing the top 25
n_tfs <- 25

my_TF <- acts %>% filter(p_value <0.05) %>% group_by(condition) %>% 
  slice_max(abs(score), n=25) %>% pull(source)

# Transform to wide matrix
sample_acts_mat <- acts %>% filter(source %in% my_TF) %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

#saveRDS(sample_acts_mat, 'tf_25_cluster.rds')

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))
library(pheatmap)
# Plot
p1 <- pheatmap(sample_acts_mat, border_color = NA, color=my_color, 
               breaks = my_breaks, main = 'Transcriptional factor activity', treeheight_col = 0)

pdf("/mnt/Data/shameed/Olbretch/figures/TF_pseudobulk.pdf", width = 10.0, height = 3.5)
p1
dev.off()

####Differentially expressed TF activity (Doublets Vs Doublets)###################################
Idents(FCan_Data)<- FCan_Data$TF_class
##CR Vs PD
deg <- FindMarkers(FCan_Data, ident.1 = 'CR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PD_Cancer cells_Fibroblast/Stromal cells',min.pct = 0.25)
##CR Vs PR
deg <- FindMarkers(FCan_Data, ident.1 = 'CR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PR_Cancer cells_Fibroblast/Stromal cells',min.pct = 0.25)
##PR Vs PD
deg <- FindMarkers(FCan_Data, ident.1 = 'PR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PD_Cancer cells_Fibroblast/Stromal cells',min.pct = 0.25)

deg<- deg[deg$p_val_adj <= 0.05,]
# Run ulm
contrast_acts <- run_ulm(mat=deg[, 'avg_log2FC', drop=FALSE], net=collectri, .source='source', .target='target',
                         .mor='mor', minsize = 5)

#saveRDS(contrast_acts, 'contrast_acts.CR_PD.rds')
#saveRDS(contrast_acts, 'contrast_acts.CR_PR.rds')
#saveRDS(contrast_acts, 'contrast_acts.PR_PD.rds')

###viewing the top 30
n_tfs <- 35
# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)

# Plot
p1<-ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors") + ggtitle('Partial remission Vs Progressive diseases')

pdf("/mnt/Data/shameed/Olbretch/figures/TF_PR_PD.pdf", width = 7.5, height = 3.0)
p1
dev.off()


#############Differentially expressed TF activity (Doublets Vs Singlets)###################################

Idents(FCan_CR)<- FCan_CR$neigb_class
deg <- FindMarkers(FCan_CR, ident.1 = 'Doublet',ident.2 = 'Singlet',min.pct = 0.25)
Idents(FCan_PD)<- FCan_PD$neigb_class
deg <- FindMarkers(FCan_PD, ident.1 = 'Doublet',ident.2 = 'Singlet',min.pct = 0.25)
Idents(FCan_PR)<- FCan_PR$neigb_class
deg <- FindMarkers(FCan_PR, ident.1 = 'Doublet',ident.2 = 'Singlet',min.pct = 0.25)

deg<- deg[deg$p_val_adj <= 0.05,]
# Run ulm
contrast_acts <- run_ulm(mat=deg[, 'avg_log2FC', drop=FALSE], net=collectri, .source='source', .target='target',
                         .mor='mor', minsize = 5)
contrast_act
s
#saveRDS(contrast_acts, 'contrast_acts_Doub_Sing_PR.rds')

###viewing the top 30
n_tfs <- 30
# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)

# Plot
p9<-ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors") + ggtitle('Progressive diseases: Doublet Vs Singlet')

pdf("/mnt/Data/shameed/Olbretch/figures/TF_doubSing_PD.pdf", width = 7.5, height = 3.0)
p9
dev.off()
####Repeat plots for CR and PR groups!!

########PATHWAY ANALYSIS############################################################################################
FCan_Data$TF_class<- paste0(FCan_Data$survival, '_', FCan_Data$neigb_type)
Idents(FCan_Data)<- FCan_Data$TF_class
table(FCan_Data$TF_class)
###Doublets Vs Doublets####
deg_C_PD <- FindMarkers(FCan_Data, ident.1 = 'CR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PD_Cancer cells_Fibroblast/Stromal cells',min.pct = 0.25)
deg_C_PR <- FindMarkers(FCan_Data, ident.1 = 'CR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PR_Cancer cells_Fibroblast/Stromal cells',min.pct = 0.25)
deg_PR_PD <- FindMarkers(FCan_Data, ident.1 = 'PR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PD_Cancer cells_Fibroblast/Stromal cells',min.pct = 0.25)

deg_C_PD<- deg_C_PD[deg_C_PD$p_val_adj <= 0.05,]###3016
deg_C_PR<- deg_C_PR[deg_C_PR$p_val_adj <= 0.05,]###1081
deg_PR_PD<- deg_PR_PD[deg_PR_PD$p_val_adj <= 0.05,]###2486

#H<- msigdbr::msigdbr(species = 'Homo sapiens', category = 'H')
#C2_KEGG <-msigdbr:: msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
C2_REACTOME <-msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
#C5 <-msigdbr:: msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
#C7 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7") #immunologic signature gene sets
#C8 <- msigdbr(species = "Homo sapiens", category = "C8") #cell type signature gene sets

C.symbol<-C2_REACTOME%>% select(c(gs_name, gene_symbol)) %>% group_by(gs_name) %>%
  summarise(all.genes= list(gene_symbol)) %>% deframe()
gsea_gene<- deg_PR_PD %>% #### repeat by replacing with DEGs of the other groups
  rownames_to_column(var = 'genes') %>% 
  dplyr::select(c( genes, avg_log2FC))
gsea_gene<- gsea_gene[order(-gsea_gene$avg_log2FC),] #arrange(desc(gsea_gene$log2FoldChange))
gene_list<- gsea_gene$avg_log2FC
#gene_list<- jitter(gsea_gene$avg_log2FC, factor = 0.01) ##if there is error due to many matches in avg_log2FC
names(gene_list)<- gsea_gene$genes
gsea_path<- fgseaSimple(pathways = C.symbol, stats = gene_list, nperm = 1000)
#gsea_path$pathway<- str_replace(gsea_path$pathway, 'GOBP_', '')
gsea_path <- gsea_path%>% filter(pval <0.05)

ggplot(gsea_path[1:40], aes(x= reorder(pathway, NES), y= NES, fill= pval)) + 
  geom_col() + coord_flip() + labs(y='normalised enrichment scores',
                                   x= 'Biological processes') + ggtitle('Doublet Partial remission Vs Progressive disease')
CR_PD_reac <- gsea_path
CR_PR_reac <- gsea_path
PR_PD_reac <- gsea_path

CR_PD_reac$group <- 'CR Vs PD'
PR_PD_reac$group <- 'PR Vs PD'
CR_PR_reac$group <- 'CR Vs PR'

#saveRDS(PR_PD_reac, 'PR_PD_reactome.rds')

Pathways_Reac <- rbind(CR_PD_reac, CR_PR_reac, PR_PD_reac)
Pathways_filt<- Pathways_Reac %>% group_by(group) %>% slice_max(n=20, order_by = abs(NES)) %>% pull(pathway)
Pathways_Reac <- Pathways_Reac[Pathways_Reac$pathway %in% unique(Pathways_filt),]

Pathways_Reac$pathway <- gsub('REACTOME_', '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub('_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND', '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub("_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS", '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub('AND_EIFS_AND_SUBSEQUENT_BINDING_TO_4', '', Pathways_Reac$pathway)

#saveRDS(Pathways_Reac, 'DEG_combined_reactome.rds')

p1<- ggplot(Pathways_Reac, aes(x = factor(group, levels = c("CR Vs PD", "PR Vs PD","CR Vs PR")), y = pathway)) + 
  geom_point(aes(fill = NES, size = -log10(padj)), shape = 21) + 
  labs( x= "", y = "", fill = "NES", size = "-log10(padj)")  + 
  scale_fill_gradient2(low="#3B4992FF", mid="white", high="#EE0000FF", oob = scales::squish) +
  theme(legend.key=element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth  = 1.2), 
        legend.position = "right", axis.text = element_text(face = "bold", size = 9), 
        legend.text = element_text(face = "bold", size = 10), 
        legend.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 90),
        title = element_text(face = "bold", size = 12)) #+ 
ggtitle('Reactome pathway analysis: Cancer_Fibroblasts Doublets DEGs')

pdf("/mnt/Data/shameed/Olbretch/figures/doublet_DEG_pathways.pdf", width = 8.5, height = 6.2)
p1
dev.off()

#######Doublet vs singlet##########

#H<- msigdbr::msigdbr(species = 'Homo sapiens', category = 'H')
#C2_KEGG <-msigdbr:: msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
C2_REACTOME <-msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
#C5 <-msigdbr:: msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
#C7 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7") #immunologic signature gene sets
#C8 <- msigdbr(species = "Homo sapiens", category = "C8") #cell type signature gene sets

C.symbol<-C2_REACTOME%>% select(c(gs_name, gene_symbol)) %>% group_by(gs_name) %>%
  summarise(all.genes= list(gene_symbol)) %>% deframe()
gsea_gene<- deg_PD_DS %>% ####repeat by replacing with DEGs of the other groups
  rownames_to_column(var = 'genes') %>% 
  dplyr::select(c( genes, avg_log2FC))
gsea_gene<- gsea_gene[order(-gsea_gene$avg_log2FC),] #arrange(desc(gsea_gene$log2FoldChange))
gene_list<- gsea_gene$avg_log2FC
#gene_list<- jitter(gsea_gene$avg_log2FC, factor = 0.01) ##if there is error due to many matches in avg_log2FC
names(gene_list)<- gsea_gene$genes
gsea_path<- fgseaSimple(pathways = C.symbol, stats = gene_list, nperm = 1000)
#gsea_path$pathway<- str_replace(gsea_path$pathway, 'GOBP_', '')
gsea_path <- gsea_path%>% filter(pval <0.05)

ggplot(gsea_path, aes(x= reorder(pathway, NES), y= NES, fill= pval)) + 
  geom_col() + coord_flip() + labs(y='normalised enrichment scores',
                                   x= 'Biological processes') + ggtitle('Doublet Partial remission Vs Progressive disease')

PD_reac <- gsea_path
CR_reac <- gsea_path
PR_reac <- gsea_path  ####only 6 processes

PD_reac$group <- 'PD'
PR_reac$group <- 'PR'
CR_reac$group <- 'CR'

#saveRDS(CR_reac, 'CR_Doub_Sing_reactome.rds')

Pathways_Reac <- rbind(CR_reac,PR_reac, PD_reac)
Pathways_filt <- Pathways_Reac %>% group_by(group) %>% slice_max(n=30, order_by = abs(NES)) %>% pull(pathway)
Pathways_Reac <- Pathways_Reac[Pathways_Reac$pathway %in% Pathways_filt,]
#Pathways_Reac <- Pathways_Reac %>% arrange(desc(abs(NES)))
#Pathways_Reac <- rbind(Pathways_Reac[1:30], Qian_Reac)
Pathways_Reac$pathway <- gsub('REACTOME_', '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub('_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND', '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub("_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS", '', Pathways_Reac$pathway)
Pathways_Reac$pathway <- gsub('AND_EIFS_AND_SUBSEQUENT_BINDING_TO_4', '', Pathways_Reac$pathway)

#saveRDS(Pathways_Reac, 'Doub_Sing_combined_reactome.rds')

ggplot(Pathways_Reac, aes(x = group, y = pathway)) + 
  geom_point(aes(fill = NES, size = -log10(padj)), shape = 21) + 
  labs( x= "", y = "", fill = "NES", size = "-log10(padj)")  + 
  scale_fill_gradient2(low="#3B4992FF", mid="white", high="#EE0000FF", oob = scales::squish) +
  theme(legend.key=element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth  = 1.2), 
        legend.position = "right", axis.text = element_text(face = "bold", size = 9), 
        legend.text = element_text(face = "bold", size = 10), 
        legend.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 90),
        title = element_text(face = "bold", size = 12)) + 
  ggtitle('Reactome: Cancer_Fibroblasts Doublets Vs Singlet')


##################TRAJECTORY INFERENCES################################################################################

#install.packages("devtools")
#devtools::install_github('cole-trapnell-lab/monocle3')
#devtools::install_github('satijalab/seurat-wrappers')
install.packages("igraph")
library(igraph)
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)

Idents(Tubo_tumor_annotated)= Tubo_tumor_annotated$neigb_type
FCan_Data<- subset(Tubo_tumor_annotated, idents=c('Fibroblast/Stromal cells', 'Cancer cells', 'Cancer cells_Fibroblast/Stromal cells'))
Idents(FCan_Data)= FCan_Data$neigb_class
FCan_doub <- subset(FCan_Data, idents='Doublet')
FCan_doub$TF_class <- paste0(FCan_doub$survival, '_', FCan_doub$neigb_type)
Idents(FCan_doub) <- FCan_doub$TF_class

#####reclustering doublets
FCan_doub<-ScaleData(FCan_doub, features = rownames(FCan_doub))
FCan_doub<-FindNeighbors(FCan_doub,dims = 1:10)
FCan_doub<-FindClusters(FCan_doub,resolution = 0.5)
FCan_doub<-RunUMAP(FCan_doub, dims = 1:10)
DimPlot(FCan_doub,reduction = "umap")
plot1 <-DimPlot(FCan_doub, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'TF_class') 
plot2<-DimPlot(FCan_doub, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'seurat_clusters') + NoLegend()

plot1 | plot2

####trajectory: create cell dataset object and cluster cells
FCan_doub@meta.data <- FCan_doub@meta.data %>% select(-sample_name)
cds <-as.cell_data_set(FCan_doub)
cds <-cluster_cells(cds, resolution=1e-3)
fData(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

#####if more than one partitions, subset the partition of interest
integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1) ## all cells belong to partition 1
cds <- as.cell_data_set(integrated.sub)

######learn pseudotime trajectory pattern
cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)
pseudotime(cds)

plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "TF_class",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size = 4)

#####order by pseudotime by choosing the root cluster based on biological knowledge
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 3]))
plot3<-plot_cells(cds,
                  color_cells_by = "pseudotime",
                  group_cells_by = "cluster",
                  label_cell_groups = F,
                  label_groups_by_cluster=F,
                  label_leaves=F,
                  label_branch_points=F,
                  label_roots = F,
                  trajectory_graph_color = "grey60")

plot4 <-DimPlot(FCan_doub, reduction = "umap", label = F, pt.size = 0.5, group.by = 'TF_class') 

plot4 | plot3

integrated.sub <- as.Seurat(cds, assay = NULL)
Idents(integrated.sub) <- integrated.sub$TF_class
FeaturePlot(integrated.sub, "monocle3_pseudotime", label = T)

#########visualise pseudotime distribution by cell clusters/labels
cds$pseudotimr_values <- pseudotime(cds)
pseudo_data <- as.data.frame(colData(cds))
ggplot(pseudo_data, aes(pseudotimr_values, reorder(TF_class,pseudotimr_values, median),
                        fill= TF_class)) + geom_boxplot() +
  labs(x = 'pseudotime',
       y= '') + coord_flip()

######perform a differential expression to know gene trajectory patterns
cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)

#saveRDS(cds_graph_test_results, 'cds_graph_test_results.rds')

rowData(cds)$gene_short_name <- row.names(rowData(cds))
deg_pseudotime <- subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, 
                                                      decreasing = TRUE),], q_value < 0.05)
deg_ids <- rownames(deg_pseudotime)

####view some of the genes along the trajectory
plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)

FeaturePlot(integrated.sub, features = "RPL34")

########identify modules of coexpressed genes
cds <- preprocess_cds(cds, num_dim = 100)
gene_module_df <- find_gene_modules(cds[deg_ids,], resolution=1e-2)

#saveRDS(gene_module_df, 'gene_module_df.rds')

plot_cells(cds, 
           genes=gene_module_df %>% filter(module %in% c(8, 28, 33, 37)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)

cell_groups <- data.frame(cell = row.names(colData(cds)),
                          cell_group = colData(cds)$TF_class)
agg_mat <- aggregate_gene_expression(cds,
                                     gene_group_df = gene_module_df,
                                     cell_group_df = cell_groups)
dim(agg_mat)

#saveRDS(agg_mat, 'agg_mat.rds')

###plotting modules by survival groups
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- str_replace(colnames(agg_mat),'_Cancer cells_Fibroblast/Stromal cells', '') 
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize_col = 15,
                   fontsize_row = 12,
                   main = 'Enrichment patterns of co-expressed gene modules'
)

#############subset genes and modules for CR and PD groups
row.names(agg_mat) <- str_replace(row.names(agg_mat), "Module ", '')
CR_modules <- agg_mat %>% as.data.frame() %>% filter(CR > PD & PR > PD) %>% rownames()
CR_trajectory <- gene_module_df %>% filter(module %in% CR_modules)
gene_module_df$module %in% CR_modules %>% sum()
PD_modules <- agg_mat %>% as.data.frame() %>% filter(PD >CR & PD > PR) %>% rownames()
PD_trajectory <- gene_module_df %>% filter(module %in% PD_modules)

#saveRDS(CR_trajectory, 'CR_trajectory.rds')
#saveRDS(PD_trajectory, 'PD_trajectory.rds')

saveRDS(cds, 'cds.rds')


######functional annotation of module genes#################
databases<- c('GO_Molecular_Function_2021', 'GO_Biological_Process_2021', 'KEGG_2019_Human',
              'KEA_2013', 'Reactome_2022','MSigDB_Oncogenic_Signatures', 'MSigDB_Hallmark_2020')

upregulated<-CR_trajectory$id
upregulated<-PD_trajectory$id

enriched_up<- enrichr(upregulated, databases = databases)
head(enriched_up$KEGG_2019_Human, 20)
head(enriched_up$GO_Molecular_Function_2021,20)
head(enriched_up$GO_Biological_Process_2021, 20)
head(enriched_up$MSigDB_Oncogenic_Signatures,20)
head(enriched_up$Reactome_2022,20)

plotEnrich(enriched_up$Reactome_2022, showTerms = 30, numChar = 60, xlab = 'Recactome Enriched terms', title = 'CR: Trajectory module genes')
plotEnrich(enriched_up$KEGG_2019_Human, showTerms = 30, numChar = 40, xlab = 'KEGG Enriched pathways', title = 'CR: Trajectory module genes')
plotEnrich(enriched_up$GO_Molecular_Function_2021, showTerms = 30, xlab = 'Molecular function', title =  'Complete remission')
plotEnrich(enriched_up$GO_Biological_Process_2021, showTerms = 30, xlab = 'Biological processes', numChar = 50,title = 'CR: Trajectory module genes')


