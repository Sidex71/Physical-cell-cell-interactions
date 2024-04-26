######PACKAGES#######################################
library(Seurat)
library(patchwork)
library(hdf5r)
library(tidyverse)
library(Neighborseq)
library(RColorBrewer)
library(Neighborseq)
library(RColorBrewer)
library('celldex')
library('SingleR')
library(ggpubr)
library(igraph)

message("Loading matrix...")
tuboData<- readRDS('/mnt/Data/shameed/Olbretch/2095-Olbrecht_counts_matrix.rds')
meta <- read_csv("/mnt/Data/shameed/Olbretch/2093-Olbrecht_metadata.csv")
colnames(tuboData)
message("Metadata per cell...")
my_meta<- data.frame(cell_label=colnames(tuboData))
my_meta <- my_meta %>% separate(., col = cell_label, into = c('cell', 'sample_name'), sep = '_', remove = F)
my_meta <- my_meta %>% select(-cell)
meta_seurat <- left_join(x = my_meta,
                         y = as.data.frame(meta),
                         by = "sample_name")
# add cell labels as rownames
rownames(meta_seurat) <- meta_seurat$cell_label
rownames(tuboData) <- gsub("_", "-", rownames(tuboData))

message("Creating Seurat object...")
tubo_ovarian<- CreateSeuratObject(tuboData,
                                  project              = "OV",
                                  min.cells            = 10,     # only genes > 10 cells
                                  min.genes            = 200,    # only cells with > 200 genes
                                  is.expr              = 0,      # expression threshold for 'detected' gene
                                  normalization.method = NULL,   # no normalization for now
                                  scale.factor         = 10000,  # scale factor in normalization process
                                  # (not used given norm.method == NULL)
                                  do.scale             = FALSE,  # gene-based z-score
                                  do.center            = FALSE,  # gene-based centering
                                  names.field          = 1,
                                  names.delim          = "_",
                                  meta.data            = meta_seurat,
                                  display.progress     = TRUE)

tubo_ovarian[["percent.mt"]] <-PercentageFeatureSet(tubo_ovarian, pattern = "^MT-")
VlnPlot(tubo_ovarian,features = c("nFeature_RNA",
                                  "nCount_RNA",
                                  "percent.mt"),ncol = 3)
plot1<- FeatureScatter(tubo_ovarian,feature1 ="nCount_RNA",
                       feature2 = "percent.mt")
plot2<- FeatureScatter(tubo_ovarian,feature1 ="nCount_RNA", 
                       feature2 = "nFeature_RNA")
#view(tubo_ovarian@meta.data)
plot1 + plot2  
tubo_ovarian$tissue<- tubo_ovarian$sample_type
VlnPlot(tubo_ovarian, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = 'tissue')
ggplot(tubo_ovarian@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, col = percent.mt)) +
  geom_point(size = 0.8) +
  scale_colour_viridis_c(option = 'F') + 
  lims(x = c(0, NA), y = c(0, NA)) +
  facet_wrap(~patient_id, nrow = 5, scales = 'free') +
  theme_minimal()
tubo_ovarian<-subset(tubo_ovarian,subset = nFeature_RNA >300 & nFeature_RNA <5000
                     & percent.mt < 25)
tubo_ovarian<-NormalizeData(tubo_ovarian, normalization.method = "LogNormalize",
                            scale.factor = 10000)
tubo_ovarian<-FindVariableFeatures(tubo_ovarian,selection.method = "vst",
                                   nfeatures = 5000)
#top20<-head(VariableFeatures(tubo_ovarian),20)
#plot3<- VariableFeaturePlot(tubo_ovarian)
#plot4<- LabelPoints(plot = plot3,points = top20)
#plot3 + plot4
tubo_ovarian<-ScaleData(tubo_ovarian, features = rownames(tubo_ovarian))
tubo_ovarian<-RunPCA(tubo_ovarian, features = VariableFeatures(object = tubo_ovarian))
#print(tubo_ovarian[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(tubo_ovarian, dims = 1:2, reduction = "pca")
DimPlot(tubo_ovarian, reduction = "pca")
ElbowPlot(tubo_ovarian)
tubo_ovarian<-FindNeighbors(tubo_ovarian,dims = 1:20)
tubo_ovarian<-FindClusters(tubo_ovarian,resolution= 0.35)
tubo_ovarian<-RunUMAP(tubo_ovarian, dims = 1:20)
DimPlot(tubo_ovarian,reduction = "umap", label=T) + NoLegend()
#DimPlot(tubo_ovarian,group.by = 'RNA_snn_res.0.3')
DimPlot(tubo_ovarian, reduction = 'umap', group.by = c('sample_name', 'patient_id', 'tissue', 'sample_site'), ncol = 2)
#saveRDS(tubo_ovarian, file = "./tubo_ovarian_tutorial.rds")
library('celldex')
library('SingleR')
sc_data<- tubo_ovarian@assays$RNA@data
ref1 <-celldex::HumanPrimaryCellAtlasData()

celltype_hpca_fine<- SingleR(test = sc_data, ref = ref1, assay.type.test=1,
                             labels = ref1$label.fine)
celltype_hpcamain<-SingleR(test = sc_data, ref = ref1, assay.type.test=1,
                           labels = ref1$label.main)
ref_encode<-celldex:: BlueprintEncodeData() 
celltype_encode_main<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.main)
celltype_encode_fine<-SingleR(test = sc_data, ref = ref_encode, assay.type.test=1,
                              labels = ref_encode$label.fine)
tubo_ovarian$celltype_hpca<- celltype_hpcamain$pruned.labels
tubo_ovarian$celltype_main<- celltype_encode_main$pruned.labels
tubo_ovarian$celltype_fine<- celltype_encode_fine$pruned.labels
#p2<-DimPlot(tubo_ovarian, reduction = 'umap', group.by = 'celltype_fine', label = T, repel = T)
p3<-DimPlot(tubo_ovarian, reduction = 'umap', group.by = 'celltype_main', label = T, repel = T)
#p4<-DimPlot(tubo_ovarian, reduction = 'umap', group.by = 'celltype_hpca', label = T, repel = T)
#p3 +p4
#x11()
#view(table(tubo_ovarian$celltype_main))
#view(table(tubo_ovarian$celltype_fine))

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
## DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = c("Immune system") # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = tubo_ovarian[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(tubo_ovarian@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(tubo_ovarian@meta.data[tubo_ovarian@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(tubo_ovarian@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
#sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
#print(sctype_scores[,1:3])
tubo_ovarian@meta.data$ScType_cell = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  tubo_ovarian@meta.data$ScType_cell[tubo_ovarian@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
#########cell reannotation#####
VlnPlot(tubo_ovarian, features = c('STAR','FOXL2','DLK1', 'ARX', 'LUM'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('PDGFRA','COL1A1','COL1A2', 'BGN', 'DCN','POSTN','COL6A1'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('FABP4','ADIPOQ','ACSL1', 'PDGFRA', 'PDGFRB'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('PDGFRA','COL1A1','COL1A2', 'BGN', 'DCN','COL6A1'), log = T, group.by = 'ScType_cell')
VlnPlot(tubo_ovarian, features = c('CD3D','CD3E', 'TRAC', 'CD19', 'CD79A','IGKC'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('CD1C','CD207','CLEC9A', 'CD1A', 'HLA-DRA','CLEC10A','LILRA4'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('CD1C','ITGAX', 'CD1A', 'CD1B','CD207','ITGAM','CD209'), log = T, group.by ='seurat_clusters')
VlnPlot(tubo_ovarian, features = c('CD68','CD14','FCGR3A', 'LYZ', 'MARCO','CD86','HLA-DRB1'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('CD8A','PRF1','NCR1', 'KLRC1', 'KLRC3'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('PECAM1','CLDN5','VWF', 'MCAM', 'ADGRF5', 'CD93'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('EPCAM','KRT18', 'KRT6A','TNNT2','BAMBI'), log = T, group.by = 'celltype_main')
VlnPlot(tubo_ovarian, features = c('ERMAP','HBA-A2', 'HBA-X','HBA-BT','BAMBI'), log = T, group.by = 'ScType_cell')

#############cell refinement#############
tubo_ovarian$refined_celltype<- str_replace(tubo_ovarian$celltype_main,'HSC','Cancer cells') %>%
  str_replace('Mesangial cells','Cancer cells') %>% str_replace('Monocytes', 'DC') %>%
  str_replace('Adipocytes', 'Fibroblast/Stromal cells') %>% str_replace('NA','Cancer cells') %>%
  str_replace('Astrocytes', 'Fibroblast/Stromal cells') %>% str_replace('Neurons', 'Fibroblast/Stromal cells') %>%
  str_replace('Epithelial cells', 'Cancer cells') %>% str_replace('Chondrocytes','Fibroblast/Stromal cells') %>%
  str_replace('Keratinocytes', 'Cancer cells') %>% str_replace('Myocytes','Fibroblast/Stromal cells') %>%
  str_replace('Skeletal muscle', 'Fibroblast/Stromal cells') %>% str_replace('Smooth muscle', 'Fibroblast/Stromal cells') %>%
  str_replace('Erythrocytes', 'Fibroblast/Stromal cells') %>%str_replace('Fibroblasts', 'Fibroblast/Stromal cells') %>%
  str_replace('Chondrocytes', 'Fibroblast/Stromal cells') %>% str_replace('Melanocytes', 'Fibroblast/Stromal cells')%>%
  replace_na('Cancer cells')
DimPlot(tubo_ovarian, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'refined_celltype')

markers <- c('STAR','FOXL2','DLK1', 'ARX', 'LUM','PDGFRA','COL1A1','COL1A2', 'BGN', 'DCN','POSTN','COL6A1',
             'FABP4','ADIPOQ','ACSL1', 'PDGFRB','CD3D','CD3E','CD4', 'CD8A', 'TRAC', 'CD19', 'CD79A','IGKC',
             'CD1C','CD207','CLEC9A', 'CD1A', 'HLA-DRA','CLEC10A','LILRA4',
             'CD68','CD14','FCGR3A', 'LYZ', 'MARCO','CD86','HLA-DRB1',
             'PRF1','NCR1', 'KLRC1','PECAM1','CLDN5','VWF', 'MCAM', 'ADGRF5', 'CD93',
             'EPCAM','KRT18', 'KRT6A','TNNT2','BAMBI')
p5<- DotPlot(tubo_ovarian_annotated, features = markers, group.by = 'refined_celltype') +
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = 'bottom') #+ coord_flip()
pdf("/mnt/Data/shameed/Olbretch/figures/marker_expression_2.pdf", width = 8.5, height = 3.5)
p5
dev.off()


####create new cluster labels for neigborseq analysis
tubo_ovarian_annotated <- tubo_ovarian
tubo_ovarian_annotated@meta.data<- rownames_to_column(tubo_ovarian_annotated@meta.data, var = 'my_rows')
tubo_ovarian_annotated$my_newrows<- tubo_ovarian_annotated$my_rows
tubo_ovarian_annotated@meta.data<- column_to_rownames(tubo_ovarian_annotated@meta.data, var = 'my_rows')

tubo_ovarian_annotated@meta.data<- tubo_ovarian_annotated@meta.data %>% 
  mutate(new.clusters= paste0(seurat_clusters, '_', refined_celltype))

test.cell<-as.data.frame(table(tubo_ovarian_annotated$new.clusters)) %>% 
  filter(Freq >40) %>% rename('new.clusters'= Var1) %>% select(-Freq) %>% 
  mutate(final_clusters=0:27) #%>% left_join(tubo_ovarian_annotated@meta.data, by='new.clusters')
#test.cell$final_clusters<- as.factor(test.cell$final_clusters)
tubo_ovarian_annotated<- subset(tubo_ovarian_annotated, subset= new.clusters %in% 
                                  test.cell$new.clusters)
tubo_ovarian_annotated@meta.data<- left_join(tubo_ovarian_annotated@meta.data,test.cell, by='new.clusters')

tubo_ovarian_annotated@meta.data<- column_to_rownames(tubo_ovarian_annotated@meta.data, var = 'my_newrows')
view(tubo_ovarian_annotated@meta.data)

###join survival information
response_info<- data.frame(patient_id=c('P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7'),
                           survival= c('CR','CR','PD', 'CR', 'PR','CR','PR'),
                           stage= c('IIIC', 'IVB','IVB', 'IC1', 'IVB', 'IIIC', 'IVB'))

meta.data <- left_join(tubo_ovarian_annotated@meta.data, response_info,
                       by = 'patient_id')
rownames(meta.data)<- meta.data$cell_label
tubo_ovarian_annotated@meta.data <- meta.data

#saveRDS(tubo_ovarian_annotated, file = "tubo_ovarian_annotated.rds")

#########plot the proportion of cell types##############################################
####tissue & cell types
tubo_bar<-tubo_ovarian_annotated@meta.data %>% group_by(tissue, refined_celltype) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n))
#saveRDS(tubo_bar, 'tubo_bar_1.rds')
p6<-ggplot(tubo_bar, aes(x=tissue, y=cell_percentage, fill= refined_celltype)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'Cell distribution by tissue type',
       y= 'Cell proportion (%)', x=NULL)+
  theme_bw() #+
theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
      axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
      legend.text = element_text(size = 17, face = 'bold'),
      legend.key.size = unit(1.2, "cm"),
      axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
      plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)

#write_csv(tubo_bar, 'cell_proportion_TvsN.csv')  
pdf("/mnt/Data/shameed/Olbretch/figures/cell_distribution_tissue.pdf", width = 4.5, height = 3.0)
p6
dev.off()
#####survival + tissue + cell types
tubo_bar<-tubo_ovarian_annotated@meta.data %>% group_by(survival, tissue, refined_celltype) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n))
tubo_bar$survival<- factor(tubo_bar$survival, levels = c('CR', 'PR', 'PD'))
ggplot(tubo_bar, aes(x=survival, y=cell_percentage, fill= refined_celltype)) + 
  geom_bar(stat = 'identity', width = 0.4) + 
  labs(title = 'Cell proportion by tissue and patient outcomes',
       y= 'Cell proportions', x=NULL)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(labels=c("Complete remission", "Partial remission",'Progressive disease'))+
  theme(axis.text.x = element_text(size=12, hjust = 0.5)) + facet_wrap(~tissue, ncol = 2)


###tissue + survival
library(ggplot2)
library(tidyverse)
tubo_bar<-tubo_ovarian_annotated@meta.data %>% group_by(survival, tissue) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n))
p7 <-ggplot(tubo_bar, aes(x=survival, y=cell_percentage, fill= tissue)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'Tissue distribution',
       y= 'Cell proportion (%)', x=NULL)+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(labels=c("CR", "PR",'PD'))+ 
  theme_bw() #+
theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
      axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
      legend.text = element_text(size = 17, hjust = 0.5, face = 'bold'),
      legend.key.size = unit(1.5, "cm"),
      axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
      plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)


#saveRDS(tubo_bar, 'tubo_bar.rds')

pdf("/mnt/Data/shameed/Olbretch/figures/tissue_distriution_survival.pdf", width = 3.5, height = 2.5)
p7
dev.off()  

########survival alone
Tubo_tumor_annotated <-subset(tubo_ovarian_annotated, subset= tissue=='Tumor') ##subsetting only cells from tumor tissues, to be used for all downstream analysis 

tubo_bar<-Tubo_tumor_annotated@meta.data %>% group_by(survival, refined_celltype) %>%
  tally() %>% mutate(cell_percentage= n*100/sum(n))
tubo_bar$survival<- factor(tubo_bar$survival, levels = c('CR', 'PR', 'PD'))
p8<-ggplot(tubo_bar, aes(x=survival, y=cell_percentage, fill= refined_celltype)) + 
  geom_bar(stat = 'identity', width = 0.2) + 
  labs(title = 'Cell proportion by patient outcomes in tumour tissue',
       y= 'Cell proportions', x=NULL)+ 
  scale_x_discrete(labels=c("CR", "PR",'PD'))+ 
  theme_bw() #+
theme(axis.text.x = element_text(size=18, hjust = 0.5, face = 'bold'),
      axis.text.y = element_text(size = 18, hjust = 0.5, face = 'bold'),
      legend.text = element_text(size = 17, face = 'bold'),
      legend.key.size = unit(1.2, "cm"),
      axis.title.y = element_text(size = 20, hjust = 0.5, face = 'bold'),
      plot.title = element_text(size = 25, hjust = 0.5, face = 'bold')) +
  scale_fill_discrete(name=NULL)

tubo_bar <- tubo_bar %>% pivot_wider(names_from = survival, values_from = c(n, cell_percentage))
write_csv(tubo_bar, 'cell_proportion_survival.csv')  

pdf("/mnt/Data/shameed/Olbretch/figures/cell_distriution_survival.pdf", width = 5.5, height = 3.0)
p8
dev.off()  
#saveRDS(Tubo_tumor_annotated, 'tubo_tumor_annotated.rds')

##########################batch correction############################################
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)
library(SeuratObject)
# install dataset
#InstallData("ifnb")

#LoadData("ifnb")
#ov_list <- SplitObject(ifnb, split.by = "stim")

ov_list<-SplitObject(tubo_ovarian_annotated, split.by = "patient_id")

cancerData.list <-lapply(X = ov_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = cancerData.list)
anchors <- FindIntegrationAnchors(object.list = cancerData.list, anchor.features = features, verbose = F)
cancerData.combined<- IntegrateData(anchorset = anchors)
DefaultAssay(cancerData.combined) <- 'integrated'
cancerData.combined <- ScaleData(cancerData.combined, verbose = FALSE)
cancerData.combined <- RunPCA(cancerData.combined, npcs = 50, verbose = FALSE)
cancerData.combined <- FindNeighbors(cancerData.combined, dims = 1:50)
cancerData.combined <- RunUMAP(cancerData.combined, dims = 1:50)
#cancerData.combined@meta.data %>% view()
cancerData.combined <- FindClusters(cancerData.combined, resolution = 0.35) 
DimPlot(cancerData.combined, reduction = 'umap', raster = F)
Idents(cancerData.combined) <- cancerData.combined$survival
#saveRDS(cancerData.combined, 'cancerData.combined.rds')
#cancerData.combined@assays$integrated@counts %>% dim()
cancerData.combined <- tubo_ovarian_corrected
p1<-DimPlot(cancerData.combined, reduction = 'umap', group.by = 'refined_celltype', label = T, repel = T)
p2<-DimPlot(cancerData.combined, reduction = 'umap', group.by = 'patient_id', label = T, repel = T)
p3<-DimPlot(cancerData.combined, reduction = 'umap', group.by = 'survival', label = T, repel = T)
p4<-DimPlot(cancerData.combined, reduction = 'umap', group.by = 'RNA_snn_res.0.35', repel = T)

DefaultAssay(cancerData.combined) <- 'integrated'
DefaultAssay(cancerData.combined) <- 'RNA'
saveRDS(cancerData.combined, 'tubo_ovarian_corrected.rds')

pdf("/mnt/Data/shameed/Olbretch/figures/seurat_clusters.pdf", width = 3.5, height = 3.5)
p4
dev.off()

pdf("/mnt/Data/shameed/Olbretch/figures/annotated_clusters.pdf", width = 6.5, height = 4.5)
p1
dev.off()

####Neighborseq analysis for doublet prediction######################################################################
set.seed(0113)
Idents(tubo_tumor_annotated)<- tubo_tumor_annotated$final_clusters
ns.data_tubo = prep_cell_mat(ge = Tubo_tumor_annotated@assays$RNA, logfc.threshold = 0.5, celltypes =Tubo_tumor_annotated$final_clusters)
#ns_tubo_test = neighborseq(ns.data_tubo$cell.mat, ns.data_tubo$celltypes, do.mroc = T, sample = Tubo_tumor_annotated_tumor$survival)
ns_tubo_test = neighborseq(ns.data_tubo$cell.mat, ns.data_tubo$celltypes, do.mroc = T)
class(tubo_ovarian_annotated$new.clusters)

plot_interactions(ns_tubo_test$result, 
                  legend.position = 'right', 
                  width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('combined tissue')

####adding Neigborseq predictions to the seurat object#########
pred<- ns_tubo_test$pred
neigb_type<- colnames(pred)[apply(pred,1, which.max)]
pred$neigb_type<- neigb_type
pred<- separate(pred, col = 'neigb_type', into = c('a','b'),sep = '_' )
colnames(cell_cluster)<- c('cellnames', 'a')
pred <- pred %>% left_join(cell_cluster, by='a')
pred<- pred[-435] 
pred<-pred %>% rename('a'=cellnames)
cell_cluster<- cell_cluster %>% rename(b='a')
pred <- pred %>% left_join(cell_cluster, by='b') 
pred<- pred[-435] %>% rename('b'=cellnames)

pred<- pred %>% mutate(neigb_type= if_else(!is.na(b), paste0(a,'_',b),a))
pred<- pred %>% mutate(neigb_class= if_else(!is.na(b), 'Doublet','Singlet'))
#colnames(Tubo_tumor_annotated@meta.data)

meta_test<- Tubo_tumor_annotated@meta.data
#meta_test<- meta_test[-grep('neigb', colnames(meta_test))]
neigb_test<- pred %>% select(c(neigb_class, neigb_type))
neigb_test<- rownames_to_column(neigb_test, var = 'new_row')
meta_test<- rownames_to_column(meta_test, var = 'my_row')
meta_test<- rownames_to_column(meta_test, var = 'new_row')
meta_test<- left_join(meta_test, neigb_test, by='new_row') %>% 
  column_to_rownames(var = 'my_row') %>% select(-new_row)
Tubo_tumor_annotated@meta.data <- meta_test
#saveRDS(Tubo_tumor_annotated, file = 'Tubo_tumor_annotated.rds')
#saveRDS(pred, file = 'predictions.rds')
#table(Tubo_tumor_annotated$refined_celltype, Tubo_tumor_annotated$seurat_clusters)
Tubo_tumor_annotated@meta.data<- rownames_to_column(Tubo_tumor_annotated@meta.data, var = 'my_rows')
Tubo_tumor_annotated$my_newrows<- Tubo_tumor_annotated$my_rows
Tubo_tumor_annotated@meta.data<- column_to_rownames(Tubo_tumor_annotated@meta.data, var = 'my_rows')
cell_cluster<- as.data.frame(table(Tubo_tumor_annotated@meta.data$refined_celltype))
cell_cluster$final_clusters2<- 0:8
cell_cluster<- cell_cluster %>% select(-Freq)
colnames(cell_cluster)[1]<- 'refined_celltype'
test.cell<- left_join(Tubo_tumor_annotated@meta.data, cell_cluster, by='refined_celltype') %>%
  column_to_rownames(var = 'my_newrows')
Tubo_tumor_annotated@meta.data<-test.cell

####completing Neigborseq: statistics and network plots#####
result = multiplet_pop_sig(pred = ns_tubo_test$pred, sample = Tubo_tumor_annotated$survival )
ns_tubo_test$result<- result
cell_cluster<- data.frame(Tubo_tumor_annotated$new.clusters, Tubo_tumor_annotated$final_clusters) %>% unique()
cell_cluster<- cell_cluster %>% separate(col=Tubo_tumor_annotated.new.clusters,
                                         into = c('number', 'cellnames'), sep = '_') 
cell_cluster<- cell_cluster[-1]
cell_cluster$Tubo_tumor_annotated.final_clusters<- as.factor(cell_cluster$Tubo_tumor_annotated.final_clusters)
names(cell_cluster)[2]<- 'Cell_1'
#view(cell_cluster)
result_new<- left_join(ns_tubo_test$result, cell_cluster, by= 'Cell_1')
result_new<-result_new[-2]
colnames(result_new)[7]<- 'Cell_1'
#result_new<- result_new %>% rename(Cell_1= Tubo_tumor_annotated.ScType_cell)

names(cell_cluster)[2]<- 'Cell_2'
result_new<- left_join(result_new, cell_cluster, by= 'Cell_2')
result_new<-result_new[-2]
colnames(result_new)[7]<- 'Cell_2'
#result_new<- result_new %>% rename(Cell_2= Tubo_tumor_annotated.ScType_cell)

plot_interactions(result_new,
                  legend.position = 'left', 
                  width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('Olbretch et al; High-grade serous tubo-ovarian cancer refined with single-cell RNA sequencing, Omentum, Ovary and Peritoneum')

result_PR<-result_new %>% filter(sample== 'PR')
result_CR<-result_new %>% filter(sample== 'CR')
result_PD<-result_new %>% filter(sample== 'PD')

#table(Tubo_tumor_annotated$sample_site)
p1<-plot_interactions(result_PR,
                      legend.position = 'left', 
                      width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('Partial Remission')

p2<-plot_interactions(result_CR,
                      legend.position = 'right',
                      width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('Complete Remission')
p3<-plot_interactions(result_PD,
                      legend.position = 'right', 
                      width_range = c(0.5,1), min_color = 'red', max_color = 'green') + ggtitle('Progressive Disease')

p1
p2
p3

pdf("/mnt/Data/shameed/Olbretch/figures/neib_plot_PR.pdf", width = 5.5, height = 2.5)
p1
dev.off()  

pdf("/mnt/Data/shameed/Olbretch/figures/neib_plot_CR.pdf", width = 7.5, height = 3.5)
p2
dev.off()  

pdf("/mnt/Data/shameed/Olbretch/figures/neib_plot_PD.pdf", width = 3.5, height = 2.5)
p3
dev.off()  

#saveRDS(result_CR, file = 'complete_response.rds')
#saveRDS(result_PR, file = 'partial_remission.rds')
#saveRDS(result_PD, file = 'progressive_disease.rds')
#saveRDS(result_new, 'result_survival.rds')


###merging the networks and plotting a heatmap#############
PR_adjusted<- result_PR %>% 
  mutate(cell_merged= paste0(Cell_1, '_', Cell_2)) %>%
  filter(!Cell_1==Cell_2 & Counts>10 & padj<0.05) %>% select(c(Counts, cell_merged))

PD_adjusted<-result_PD %>% 
  mutate(cell_merged= paste0(Cell_1, '_', Cell_2)) %>%
  filter(!Cell_1==Cell_2 & Counts>10 & padj<0.05) %>% select(c(Counts, cell_merged))

CR_adjusted<- result_CR %>% 
  mutate(cell_merged= paste0(Cell_1, '_', Cell_2)) %>%
  filter(!Cell_1==Cell_2 & Counts>10 & padj<0.05) %>% select(c(Counts, cell_merged))

colnames(PR_adjusted)[1]<- 'Count_PR'
colnames(PD_adjusted)[1]<- 'Count_PD'
colnames(CR_adjusted)[1]<- 'Count_CR'

PD_adjusted<- aggregate(PD_adjusted$Count_PD ~ PD_adjusted$cell_merged, PD_adjusted, sum)
PR_adjusted<- aggregate(PR_adjusted$Count_PR ~ PR_adjusted$cell_merged, PR_adjusted, sum)
CR_adjusted<- aggregate(CR_adjusted$Count_CR ~ CR_adjusted$cell_merged, CR_adjusted, sum)

colnames(PR_adjusted)<- c('cell_merged', 'Partial')
colnames(PD_adjusted)<- c('cell_merged', 'Progressive')
colnames(CR_adjusted)<- c('cell_merged', 'Complete')

response_data<- full_join(PR_adjusted, PD_adjusted, by='cell_merged') %>% 
  full_join(CR_adjusted, by='cell_merged')
response_data<- response_data %>% mutate(Complete= coalesce(Complete, 0),
                                         Partial= coalesce(Partial,0),
                                         Progressive= coalesce(Progressive,0))
response_data<- response_data %>% column_to_rownames(var = 'cell_merged')
colnames(response_data)<- c('Partial_remission', 'Progressive disease', 'Complete_remission')
#response_data<- column_to_rownames(response_data, var = 'cell_merged')
#response_norm <- normalised_patient_response
response_norm<- response_data %>% mutate(Partial_remission= Partial_remission*100/2381,
                                         Complete_remission= Complete_remission*100/7389,
                                         `Progressive disease`= `Progressive disease`*100/3184)
pheatmap::pheatmap(as.matrix(response_data), scale = 'row', fontsize_row = 8,
                   main = 'Response comparison',) #treeheight_col = 0, )

colnames(response_norm) <- c('PR', 'PD', 'CR')
p1<-pheatmap::pheatmap(as.matrix(response_norm), scale = 'row', fontsize_row = 9.5,
                       main = 'Normalised Response comparison') #treeheight_col = 0, )

pdf("/mnt/Data/shameed/Olbretch/figures/neigb_heatmap.pdf", width = 6.0, height = 4.5)
p1
dev.off()

p1<-pheatmap::pheatmap(as.matrix(response_norm), scale = 'row', fontsize_row = 25.5,
                       main = 'Normalised Response comparison',
                       fontsize_col = 16)


png("/mnt/Data/shameed/Olbretch/figures/neigb_heatmap_2.png", width = 15, height = 10.5, units = 'in', res = 600)
p1
dev.off()

#saveRDS(response_data, file = 'combined_patient_response.rds')
#saveRDS(response_norm, file = 'normalised_patient_response.rds')
#saveRDS(ns_tubo_test, 'neigb_result_0112.rds')
#saveRDS(ns.data_tubo, 'ns.data_0112.rds')
