library(Seurat)
library(tidyverse)

### Preparing for LDA analysis##########################
table(Tubo_tumor_annotated$survival, Tubo_tumor_annotated$neigb_type)
Idents(Tubo_tumor_annotated)= Tubo_tumor_annotated$neigb_type

###subsetting Cancer cells, stromal cells and Cancer_stromal doublets
FCan_Data<- subset(Tubo_tumor_annotated, idents=c('Fibroblast/Stromal cells', 'Cancer cells', 'Cancer cells_Fibroblast/Stromal cells'))
FCan_CR<-subset(FCan_Data, subset= survival=='CR') ##complete remission
table(FCan_CR$neigb_class)
FCan_PD<-subset(FCan_Data, subset= survival=='PD') ##progressive disease
table(FCan_PD$neigb_class)
FCan_PR<- subset(FCan_Data, subset=survival=='PR') ##partial remission
table(FCan_PR$neigb_class)
###exporting expression matrix and metadata for LDA analysis
write.csv(FCan_PD@meta.data, file = 'PD_meta.csv')
write.csv(FCan_PD@assays$RNA@counts, 'PD.csv')
write.csv(FCan_CR@meta.data, file = 'CR_meta.csv')
write.csv(FCan_CR@assays$RNA@counts, 'CR.csv')
write.csv(FCan_PR@meta.data, file = 'PR_meta.csv')
write.csv(FCan_PR@assays$RNA@counts, 'PR.csv')

#####Downstream analysis from LDA Output###############################
##importing and filtering LDA genes
CR_genes = read.csv('CR_FibCan_genes.csv') ##complete remission
CR_30<-CR_genes %>% filter(X1>5) %>% group_by(X3) %>% slice_max(n=30, order_by = X1) %>% ungroup()
colnames(CR_30)<- c('genes', 'cells', 'probability', 'topics')
PD_genes = read.csv('PD_FibCan_genes.csv') ##progressive disease
PD_30<-PD_genes %>% filter(X1>5) %>% group_by(X3) %>% slice_max(n=30, order_by = X1) %>%ungroup() 
colnames(PD_30)<- c('genes', 'cells', 'probability', 'topics')
PR_genes = read.csv('PR_FibCan_genes.csv') ##partial remission
PR_30<-PR_genes %>% filter(X1 >5) %>% group_by(X3) %>% slice_max(n=30, order_by = X1) %>%ungroup()
colnames(PR_30)<- c('genes', 'cells', 'probability', 'topics')

##normalizing by doublet counts
CR_filtered<- CR_30 %>% select(c(genes, cells)) %>% mutate(cells= cells*100/246)
PR_filtered<- PR_30 %>% select(c(genes, cells)) %>% mutate(cells= cells*100/82)
PD_filtered<- PD_30 %>% select(c(genes, cells)) %>% mutate(cells= cells*100/323)

###aggregating duplicated genes
CR_filtered<- aggregate(CR_filtered$cells ~ CR_filtered$genes, CR_filtered, mean) 
PR_filtered<- aggregate(PR_filtered$cells ~ PR_filtered$genes, PR_filtered, mean)
PD_filtered<- aggregate(PD_filtered$cells ~ PD_filtered$genes, PD_filtered, mean)
colnames(CR_filtered)<- c('genes', 'cells')
colnames(PR_filtered)<- c('genes', 'cells')
colnames(PD_filtered)<- c('genes', 'cells')

#saveRDS(CR_filtered, 'CR_filtered.rds')
#saveRDS(PR_filtered, 'PR_filtered.rds')
#saveRDS(PD_filtered, 'PD_filtered.rds')

######Enrichment analysis -enrichR#################
library(enrichR)
#BiocManager::install('pathview')
library(pathview)
library(org.Mm.eg.db)
#BiocManager::install('org.Hs.eg.db')
library(org.Hs.eg.db)
library(EBImage)
#BiocManager::install('EBImage')
library(EBImage)
library('clusterProfiler')
#install.packages("msigdbr")
library(msigdbr)
#BiocManager::install('fgsea')
library(fgsea)


databases<- c('GO_Molecular_Function_2021', 'GO_Biological_Process_2021', 'KEGG_2019_Human',
              'KEA_2013', 'Reactome_2022','MSigDB_Oncogenic_Signatures', 'MSigDB_Hallmark_2020')
my_CR<-CR_filtered %>% pull(genes)
my_PR<-PR_filtered %>% pull(genes)
my_PD<-PD_filtered %>% pull(genes)

upregulated<-my_PD
enriched_up<- enrichr(upregulated, databases = databases)
upregulated<- my_CR
enriched_up<- enrichr(upregulated, databases = databases)
upregulated<- my_PR
enriched_up<- enrichr(upregulated, databases = databases)
#head(enriched_up$KEGG_2019_Human, 20)
#head(enriched_up$GO_Molecular_Function_2021,20)
#head(enriched_up$GO_Biological_Process_2021, 20)
#head(enriched_up$MSigDB_Oncogenic_Signatures,20)
#head(enriched_up$Reactome_2022,20)

p1<- plotEnrich(enriched_up$Reactome_2022, showTerms = 30, numChar = 60, xlab = 'Recactome Enriched terms', title = 'Partial remission')
#plotEnrich(enriched_up$KEGG_2019_Human, showTerms = 30, numChar = 40, xlab = 'KEGG Enriched pathways', title = 'Complete remission')
#plotEnrich(enriched_up$GO_Molecular_Function_2021, showTerms = 30, xlab = 'Molecular function', title =  'Complete remission')
#plotEnrich(enriched_up$GO_Biological_Process_2021, showTerms = 30, xlab = 'Biological processes', numChar = 50,title = 'common LDA genes ')
pdf("/mnt/Data/shameed/Olbretch/figures/pathways_PR.pdf", width = 7.5, height = 5.5)
p1
dev.off()

####differentially expressed LDA genes###################
my_CR<-CR_filtered %>% pull(genes)
my_PR<-PR_filtered %>% pull(genes)
my_PD<-PD_filtered %>% pull(genes)

Idents(Tubo_tumor_annotated) <- Tubo_tumor_annotated$neigb_type
FCan_Data<- subset(Tubo_tumor_annotated, idents=c('Fibroblast/Stromal cells', 'Cancer cells', 'Cancer cells_Fibroblast/Stromal cells'))
FCan_Data$TF_class <- paste0(FCan_Data$survival, '_', FCan_Data$neigb_type)

Idents(FCan_Data) <- FCan_Data$TF_class

set.seed(040324)
###Progressive disease
DB_can<- FindMarkers(FCan_Data, ident.1 = 'PD_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PD_Cancer cells',min.pct = 0.25)
DB_fib<- FindMarkers(FCan_Data, ident.1 = 'PD_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PD_Fibroblast/Stromal cells',min.pct = 0.25)
DB_can_up <- DB_can %>% filter(p_val_adj <0.05 & avg_log2FC >0)
DB_fib_up <- DB_fib %>% filter(p_val_adj <0.05 & avg_log2FC >0)
PD_comm_up<- intersect(rownames(DB_fib_up), rownames(DB_can_up))
DB_can_down <- DB_can %>% filter(p_val_adj <0.05 & avg_log2FC <0)
DB_fib_down <- DB_fib %>% filter(p_val_adj <0.05 & avg_log2FC <0)
PD_comm_down<- intersect(rownames(DB_fib_down), rownames(DB_can_down))
PD_DEG <- c(PD_comm_down, PD_comm_up) ####Differentially expressed doublet genes against both cancer and stromal cells
PD_DEG_LDA <- PD_DEG[PD_DEG %in% my_PD] ##subsetting for overlaps with LDA genes
PD_DEG <- data_frame(genes = PD_DEG)
PD_DEG_LDA <- data.frame(genes = PD_DEG_LDA)
write_csv(PD_DEG_LDA, 'PD_DEG_LDA.csv')
write_csv(PD_DEG, 'PD_DEG.csv')

###checking directionality
PD_comm_down[PD_comm_down %in% my_PD] ###none
PD_comm_up[PD_comm_up %in% my_PD] ###116

#####complete remission
DB_can<- FindMarkers(FCan_Data, ident.1 = 'CR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'CR_Cancer cells',min.pct = 0.25)
DB_fib<- FindMarkers(FCan_Data, ident.1 = 'CR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'CR_Fibroblast/Stromal cells',min.pct = 0.25)
DB_can_up <- DB_can %>% filter(p_val_adj <0.05 & avg_log2FC >0)
DB_fib_up <- DB_fib %>% filter(p_val_adj <0.05 & avg_log2FC >0)
CR_comm_up<- intersect(rownames(DB_fib_up), rownames(DB_can_up))
DB_can_down <- DB_can %>% filter(p_val_adj <0.05 & avg_log2FC <0)
DB_fib_down <- DB_fib %>% filter(p_val_adj <0.05 & avg_log2FC <0)
CR_comm_down<- intersect(rownames(DB_fib_down), rownames(DB_can_down))
CR_DEG <- c(CR_comm_down, CR_comm_up)####Differentially expressed doublet genes against both cancer and stromal cells
CR_DEG_LDA <- CR_DEG[CR_DEG %in% my_CR] ##subsetting overlapping LDA genes
CR_DEG <- data_frame(genes = CR_DEG)
CR_DEG_LDA <- data.frame(genes = CR_DEG_LDA)

write_csv(CR_DEG_LDA, 'CR_DEG_LDA.csv')
write_csv(CR_DEG, 'CR_DEG.csv')

CR_comm_down[CR_comm_down %in% my_CR] ###10
CR_comm_up[CR_comm_up %in% my_CR] ###

###partial remission
DB_can<- FindMarkers(FCan_Data, ident.1 = 'PR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PR_Cancer cells',min.pct = 0.25)
DB_fib<- FindMarkers(FCan_Data, ident.1 = 'PR_Cancer cells_Fibroblast/Stromal cells',ident.2 = 'PR_Fibroblast/Stromal cells',min.pct = 0.25)
DB_can_up <- DB_can %>% filter(p_val_adj <0.05 & avg_log2FC >0)
DB_fib_up <- DB_fib %>% filter(p_val_adj <0.05 & avg_log2FC >0)
PR_comm_up<- intersect(rownames(DB_fib_up), rownames(DB_can_up))
DB_can_down <- DB_can %>% filter(p_val_adj <0.05 & avg_log2FC <0)
DB_fib_down <- DB_fib %>% filter(p_val_adj <0.05 & avg_log2FC <0)
PR_comm_down<- intersect(rownames(DB_fib_down), rownames(DB_can_down))
PR_DEG <- c(PR_comm_down, PR_comm_up) ####Differentially expressed doublet genes against both cancer and stromal cells
PR_DEG_LDA <- PR_DEG[PR_DEG %in% my_PR] ###subsetting overlapping LDA genes

PR_DEG <- data_frame(genes = PR_DEG)
PR_DEG_LDA <- data.frame(genes = PR_DEG_LDA)
#write_csv(PR_DEG_LDA, 'PR_DEG_LDA.csv')
#write_csv(PR_DEG, 'PR_DEG.csv')

###pulling all differentially expressed LDA genes
LDEG <- c(CR_DEG_LDA, PD_DEG_LDA, PR_DEG_LDA)
LDEG <- unique(LDEG)


######################survival analysis#################################
###works with R version 4.3
# Install the main RTCGA package
# Install the clinical and mRNA gene expression data packages
BiocManager::install('RTCGA.clinical')
BiocManager::install("RTCGA.mRNA")
BiocManager::install('TCGAbiolinks')
library(TCGAbiolinks)
library(SummarizedExperiment)

#####get the list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-OV')

#####Building query
query_TCGA<- GDCquery(project = 'TCGA-OV', data.category = 'Transcriptome Profiling')
output_query <- getResults(query_TCGA)

####refining the query
query_TCGA<- GDCquery(project = 'TCGA-OV', data.category = 'Transcriptome Profiling',
                      experimental.strategy = 'RNA-Seq',
                      workflow.type = 'STAR - Counts',
                      access = 'open')
output_query <- getResults(query_TCGA)

####downloading files
GDCdownload(query_TCGA)

library(RTCGA)

###prepare counts data from the downloaded files
TCGA_OV_data <-GDCprepare(query_TCGA, summarizedExperiment = T)
TCGA_OV_matrix <- assay(TCGA_OV_data, 'unstranded')
#TCGA_OV_matrix <- assay(TCGA_OV_data, 'tpm_unstrand')
#TCGA_OV_matrix <- assay(TCGA_OV_data, 'fpkm_unstrand')
#TCGA_OV_matrix <- assay(TCGA_OV_data, 'fpkm_uq_unstrand')

#######getting clinical data
library(tidyverse)
clinic_OV <- GDCquery_clinic('TCGA-OV')
clinic_OV %>% select('vital_status', 'days_to_last_follow_up', 'days_to_death')

########survival
#install.packages("survival")
#install.packages("survminer")
library('survival')
library(survminer)
library(DESeq2)

clinic_OV$deceased <- ifelse(clinic_OV$vital_status=='Alive', 'FALSE', 'TRUE')
clinic_OV$overall_survival <- ifelse(clinic_OV$vital_status =='Alive',
                                     clinic_OV$days_to_last_follow_up,
                                     clinic_OV$days_to_death)

#####extract row data
gene_metadata <- as.data.frame(rowData(TCGA_OV_data))
coldata <- as.data.frame(colData(TCGA_OV_data ))

####vst transform using DESEQ

dds <- DESeqDataSetFromMatrix(countData = round(TCGA_OV_matrix),
                              colData =coldata,
                              design = ~1)
dds <- dds[rowSums(counts(dds))>=10,] ###remove genes with <10 counts across all samples
vsd <- vst(dds, blind = F)
TCGA_OV_vst <- assay(vsd)

#####extract the gene of interest
TCGA_OV_vst <- TCGA_OV_vst %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(gene_metadata, by= 'gene_id')
#saveRDS(TCGA_OV_data, 'TCGA_OV_data.rds')
#saveRDS(TCGA_OV_matrix, 'TCGA_OV_matrix.rds')
#saveRDS(dds, 'dds.rds')
#saveRDS(vsd, 'vsd.rds')
#saveRDS(TCGA_OV_vst, 'TCGA_OV_vst.rds')
#saveRDS(clinic_OV, 'clinic_OV.rds')

##multiple variables cox regression###############
LDEG<- LDEG[!grepl('-', LDEG)]
TCGA_LDA <- TCGA_OV_vst[TCGA_OV_vst$gene_name %in% LDEG, ]
TCGA_LDA <- TCGA_LDA %>% select(c(gene_name, case_id, counts)) %>% pivot_wider(names_from = gene_name, values_from = counts)
#nchar(TCGA_LDA$case_id)
#nchar(clinic_OV$submitter_id)
TCGA_LDA$case_id <- substr(TCGA_LDA$case_id, 1, 12)
TCGA_LDA <- inner_join(TCGA_LDA, clinic_OV, by= c('case_id' = 'submitter_id'))
TCGA_LDA$deceased <- as.logical(TCGA_LDA$deceased)
#formula<- paste(colnames(TCGA_LDA)[2:11], collapse = '+')
time = TCGA_LDA$overall_survival
event <- TCGA_LDA$deceased
colnames(TCGA_LDA) #check to specify the range of columns containing gene names

set.seed(12022024)
formula<-as.formula(paste("survival::Surv(time, event) ~",paste(colnames(TCGA_LDA)[2:221], collapse = '+')))
cox_model <-survival:: coxph(formula, data = TCGA_LDA)
# Get the summary of the Cox model
cox_summary <- summary(cox_model)
#saveRDS(cox_summary, 'cox_coeff_all.rds')
# Filter the coefficients based on p-values
significant_coeffs <- cox_summary$coefficients[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05, ]
significant_coeffs_LDA<- significant_coeffs %>% as.data.frame() %>% rownames_to_column('genes')
#write_csv(significant_coeffs_LDA, 'cox_coeff_signficant.csv')
cox_summary <- as.data.frame(cox_summary$coefficients)
#write_csv(cox_summary, 'cox_coeff_all.csv')

####adjust pvalues
library(stats)
significant_coeffs_LDA$p.adj <- p.adjust(significant_coeffs_LDA$`Pr(>|z|)`, method = 'fdr')
significant_coeffs_LDA<- significant_coeffs_LDA[significant_coeffs_LDA[, "p.adj"] < 0.05, ]
#saveRDS(significant_coeffs_LDA,'/mnt/Data/shameed/Olbretch/cox_coeff_significant.rds' )


##############################selecting optimum cut-off for KM plots#################
significant_coeffs_LDA <- column_to_rownames(significant_coeffs_LDA, 'genes')
res.cut <- surv_cutpoint(TCGA_LDA, time ='overall_survival' , event = 'deceased',
                         variables = rownames(significant_coeffs_LDA))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
head(res.cat)

######multiple Log Rank and KM plots##########
to_check <- rownames(significant_coeffs_LDA)
sig_genes <- c()
p_vals <- c()
p_list <- p_list <- as.list(rep(NA, length(to_check)))
for (i in 1:length(to_check)) {
  my_gene <- to_check[i]
  long_fit <- survival::survdiff(survival::Surv(time = overall_survival, 
                                                event = deceased) ~ get(my_gene), data = res.cat)
  
  if (long_fit$pvalue > 0.05) {
    print(paste0(my_gene, " is not significant"))
  } else {
    fit <- survival::survfit(survival::Surv(time = overall_survival, 
                                            event = deceased) ~ get(my_gene), data = res.cat)
    label_1 = paste0("High (n =", as.character(fit$n[1]), ')')
    label_2 = paste0("Low (n =", as.character(fit$n[2]), ')')
    p<-ggsurvplot(fit, conf.int=F, pval=T, risk.table=F, 
                  legend.labs=c(label_1, label_2), legend.title="expression status",  
                  palette=c("dodgerblue2", "orchid2"), 
                  title=my_gene, 
                  risk.table.height=.15, legend='right')
    
    print(p)
    my_p <-long_fit$pvalue
    p_vals <- c(p_vals, my_p)
    sig_genes <- c(sig_genes, my_gene)
    p_list[[i]] <- p
  }
}
sig_genes

p_list <- p_list[!is.na(p_list)]
names(p_list) <- sig_genes

png("/mnt/Data/shameed/Olbretch/figures/PQBP1.png", width = 6.5, height = 4, units = 'in', res = 600)
p_list$PQBP1
dev.off()
final_sig <- data.frame(genes = sig_genes, p_values = p_vals )
write_csv(final_sig, 'KM_Significant_genes.csv')

#######boxplot to visualize genes in singlets and doublets#############################

library(tidyverse)
# Convert dgCMatrix to data frame
write.csv(as.data.frame(FCan_Data@assays$RNA@data), 'Olb_Data.csv')
data_df <-read_csv("/mnt/Data/shameed/Olbretch/Olb_Data.csv")
data_df <- column_to_rownames(data_df, "...1")
# Add stage_class or desired metadata
data_df<- t(data_df) 
data_df <- as.data.frame(data_df)
Olb_data$survClass <- paste0(Olb_data$survival, '_', Olb_data$neigb_type)
data_df$survival <- Olb_data$survClass
data_df$`Cell class` <-data_df$survival
#data_df$survival <- Olb_data$survival

#data_df <- data_df[Olb_data$neigb_class=='Doublet',]
library(ggplot2)
library(ggpubr)
# Plot the boxplot
ggplot(data_df, aes(survival, UBE2T, fill = survival)) + geom_boxplot()

my_comp <- list(c("CR_Cancer cells_Fibroblast/Stromal cells", "CR_Cancer cells"),
                c('PR_Cancer cells_Fibroblast/Stromal cells', 'PR_Cancer cells'), 
                c('PD_Cancer cells_Fibroblast/Stromal cells', 'PD_Cancer cells'),
                c("CR_Cancer cells_Fibroblast/Stromal cells", "CR_Fibroblast/Stromal cells"),
                c('PR_Cancer cells_Fibroblast/Stromal cells', 'PR_Fibroblast/Stromal cells'), 
                c('PD_Cancer cells_Fibroblast/Stromal cells', 'PD_Fibroblast/Stromal cells'),
                c("CR_Cancer cells_Fibroblast/Stromal cells", 'PR_Cancer cells_Fibroblast/Stromal cells'),
                c('CR_Cancer cells_Fibroblast/Stromal cells', 'PD_Cancer cells_Fibroblast/Stromal cells'), 
                c('PR_Cancer cells_Fibroblast/Stromal cells', 'PD_Cancer cells_Fibroblast/Stromal cells')) 

###plotting

p1<-ggplot(data_df, aes(`Cell class`, NEAT1, fill = `Cell class`)) + 
  geom_boxplot() +
  geom_signif(comparisons = my_comp, test = "wilcox.test", textsize = 4, 
              step_increase = 0.1, map_signif_level = T) +
  stat_compare_means(label.y = 10.5) + theme_bw() + theme(axis.text.x = element_blank()) +
  labs(y= 'Expression level') + ggtitle('NEAT1') + theme(title = element_text(face = 'bold', hjust =0.5 ),
                                                         plot.title = element_text(face = 'bold', hjust = 0.5))

p1
png("/mnt/Data/shameed/Olbretch/figures/NEAT1_box.png", width = 9, height = 5.5, units = 'in', res = 600)
p1
dev.off()

###repeat plotting code for HGMB3, PQBP1, UXT, RAN, SMARCA4, TXNRD1 & BTG3