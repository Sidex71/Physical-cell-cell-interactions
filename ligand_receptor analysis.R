#######LRP manual scoring approach###########
###subsetting doublets 
Olb_doub <- subset(FCan_Data, subset= neigb_class=='Doublet')

#####calculate average expression
average_expression <- as.list(rep(NA, length(unique(Olb_doub$TF_class))))
for(i in 1: length(unique(Olb_doub$TF_class))){
  clustername <- unique(Olb_doub$TF_class)[i]
  genecounts<- Olb_doub@assays$RNA@data[,Olb_doub$TF_class==clustername]
  gene_expression <- rowMeans(genecounts)
  average_expression[[i]] <- gene_expression
  names(average_expression)[i] <- unique(Olb_doub$TF_class)[i]
}
#####calculate enrichment Z-scores
overall_mean <- mean(unlist(Olb_doub@assays$RNA@data))
overall_sd <- sd(unlist(Olb_doub@assays$RNA@data))

overall_Zscores <- as.list(rep(NA, length(average_expression)))
for(i in 1:length(average_expression)){
  express <- average_expression[[i]]
  z_sco <- (express - overall_mean)/overall_sd
  overall_Zscores[[i]] <- z_sco
  names(overall_Zscores)[i]<- names(average_expression)[i]
}

#######ligand, receptor, and interaction scores 
interaction_list <- as.list(rep(NA, length(average_expression)))

for(i in 1:length(overall_Zscores)){
  my_Zscores<- overall_Zscores[[i]]
  my_Zscores<- as.data.frame(my_Zscores) %>% rownames_to_column('genes')
  avg_expression <- average_expression[[i]] %>% as.data.frame() %>% rownames_to_column('genes')
  colnames(avg_expression)[2] <- 'average_expression'
  
  ###ligand scores for doublets
  LR_Rowman <- read_excel("/mnt/Data/shameed/LR_Rowman.xlsx",sheet = "All.Pairs")
  ligands <-unique(LR_Rowman$Ligand.ApprovedSymbol)
  ligands<-  as.data.frame(ligands)
  ligand_scores <- ligands %>% inner_join(my_Zscores, by= c('ligands'='genes')) %>%
    rename('ligand_Zscore'=my_Zscores) %>% inner_join(avg_expression, by=c('ligands'='genes')) %>%
    rename('ligand_avg_expression'=average_expression)
  
  ###receptor scores for doublet
  receptors <-unique(LR_Rowman$Receptor.ApprovedSymbol)
  receptors <- as.data.frame(receptors)
  receptor_scores <- receptors %>% inner_join(my_Zscores, by= c('receptors'='genes'))%>%
    rename('receptor_Zscore'=my_Zscores) %>% inner_join(avg_expression, by=c('receptors'='genes')) %>%
    rename('receptor_avg_expression' =average_expression)
  
  ###LRP analysis for doublet
  LRP<- LR_Rowman[, c(1,2,4)]
  LR_merged <- left_join(receptor_scores, LRP, 
                         by= c('receptors'='Receptor.ApprovedSymbol'), 
                         multiple= 'all')
  LR_merged_B <- inner_join(ligand_scores, LR_merged, 
                            by= c('ligands'='Ligand.ApprovedSymbol'), 
                            multiple= 'all')
 
   ##Calculate interaction scores
  LR_merged_B <- LR_merged_B %>% mutate(interaction_score=
                                          sqrt(ligand_Zscore^2 + receptor_Zscore^2))
  #compute number of cells and percentages
  clustername <- unique(Olb_doub$TF_class)[i]
  total_cell <- sum(Olb_doub$TF_class==clustername)
  n_cell <- Olb_doub@assays$RNA@counts[, Olb_doub$TF_class==clustername] >0
  n_cell <- rowSums(n_cell)
  n_cell <- as.data.frame(n_cell)
  n_cell$percent_cell <- n_cell$n_cell*100/total_cell 
  n_cell <- rownames_to_column(n_cell, var = 'genes')
  
  LR_merged_B<- LR_merged_B %>% inner_join(n_cell, by =c('ligands'='genes')) %>%
    rename('percent_cell_ligand' =percent_cell,
           'n_cell_ligand' = n_cell) %>% inner_join(n_cell, by= c('receptors'= 'genes')) %>%
    rename('percent_cell_receptor' =percent_cell,
           'n_cell_receptor'= n_cell)
  
  interaction_list[[i]] <- LR_merged_B
  names(interaction_list)[i] <- unique(Olb_doub$TF_class)[i]
  
}
#saveRDS(interaction_list, 'interaction_list.rds')
interaction_list[[3]] %>% dim

#####afiltering the LRP list
interaction_filter_list <- as.list(rep(NA, length(interaction_list)))
for (i in 1:length(interaction_list)){
  LR_interaction_filt <- interaction_list[[i]] %>% 
    filter(receptor_avg_expression > 0 & ligand_avg_expression >0 &
             interaction_score > 1.5 & percent_cell_ligand >4.5 & percent_cell_receptor >4.5)
  interaction_filter_list[[i]] <- LR_interaction_filt
  names(interaction_filter_list)[i] <- names(interaction_list)[i] 
}

#saveRDS(interaction_filter_list, 'interaction_filter_list.rds')

interactions_CR <- interaction_filter_list[[1]] #%>% dim()
interactions_PD <- interaction_filter_list[[2]] #%>% dim()
interactions_PR <- interaction_filter_list[[3]] #%>% dim()

##########BulkSignalR approach############################################

Olb_doub <- subset(FCan_Data, subset= neigb_class=='Doublet')
doub_df <-as.data.frame(Olb_doub@assays$RNA@counts)
doub_CR<- doub_df[, Olb_doub@meta.data$survival=='CR']
doub_PR<- doub_df[, Olb_doub@meta.data$survival=='PR']
doub_PD<- doub_df[, Olb_doub@meta.data$survival=='PD']

library(BulkSignalR)
library(igraph)
library(dplyr)
library(doParallel)
n.proc <- 2
cl <- makeCluster(n.proc)
registerDoParallel(cl)
# To add at the end of your script
# stopCluster(cl)

bsrdm_doub <- prepareDataset(counts = doub_PD, normalize = F, method = 'UQ')

# print object
str(bsrdm_doub)
## adding null distribution
bsrdm_doub<- learnParameters(bsrdm_doub, filename = "sdc")
str(bsrdm_doub)
## Building a BSRInference object
bsrinf_doub <- initialInference(bsrdm_doub)
LRinter.dataframe_doub <- LRinter(bsrinf_doub)
LRinter.dataframe_doub$pairs <- paste0(LRinter.dataframe_doub$L, '_', LRinter.dataframe_doub$R)
LRinter.dataframe_filt <- LRinter.dataframe_doub %>% filter(qval <=0.01)
###check overlap between manual scoring LRPs and BulkSignalR LRPs
PD_pairs <- LRinter.dataframe_filt[LRinter.dataframe_filt$pairs %in% 
                                     interactions_PD$Pair.Name,] ###15 unique

saveRDS(LRinter.dataframe_doub, 'PD_LRinter.dataframe_doub.rds')
saveRDS(LRinter.dataframe_filt, 'PD_LRinter.dataframe_filt.rds')
saveRDS(PD_pairs, 'PD_pairs.rds')

bsrinf_doub@LRinter <- PD_pairs

####reductions
bsrinf.redP <- reduceToPathway(bsrinf_doub) #each pathway appear once (unique pathways)
bsrinf.redP_df <- LRinter(bsrinf.redP)
bsrinf.redBP    <- reduceToBestPathway(bsrinf) ##each ligand receptor appears once (unique L-R pairs)
bsrinf.redBP_df <- LRinter(bsrinf.redBP)
bsrinf.L    <- reduceToLigand(bsrinf) ### unique ligands
bsrinf.L_df <- LRinter(bsrinf.L)
bsrinf.R    <- reduceToReceptor(bsrinf) ##unique receptors
bsrinf.R_df <- LRinter(bsrinf.R)

####Scoring by ligand-receptor gene signatures

bsrsig.redBP <- getLRGeneSignatures(bsrinf.redP, qval.thres=max(LRinter.dataframe_doub$qval))

scoresLR <- scoreLRGeneSignatures(bsrdm_doub, bsrsig.redBP,
                                  name.by.pathway=T) ###set False to obtain L-R scores

colnames(scoresLR)<-  gsub("[^0-9]", "", colnames(scoresLR))
rownames(scoresLR) <- gsub("[{}]", "", rownames(scoresLR))

######################
pheatmap(scoresLR,
         scale = "row", 
         show_rownames = T, show_colnames = FALSE,
         border_color = NA,
         fontsize_row = 9,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cutree_cols = 2, cutree_rows = 2,
         #treeheight_col = 0,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         main = 'Downstream pathways PD')

############repeat codes for CR and PR groups!!
