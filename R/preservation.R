

source("~/devel/Diversity_Transcriptome/R/utils.R")

library(WGCNA)
library(tidyverse)
library(OCMSutility)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)  
library(pheatmap)
library(circlize)
library(corrr)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(viridis)
library(stringr)



ileum_normed_dat_corrected <- readRDS("export.dir/ileum_normed_corrected.Rds")
liver_normed_dat <- readRDS("export.dir/liver_normed.Rds")
results <- readRDS("export.dir/results.Rds")

# metadata
ileum_metadata <- read.table('/well/powrie/users/cqh100/work/mouse_rna_seq/splicing/ileum_metadata_with_PC_correction.tsv', sep = '\t', header = TRUE)

ileum_metadata <- subset(ileum_metadata, age == 'Young')
ileum_metadata <- subset(ileum_metadata, group == 'Solo')
ileum_metadata$status <- as.factor(ileum_metadata$status)
rownames(ileum_metadata) <- ileum_metadata$sample_id
colnames(ileum_metadata)[1] <- 'sampleID'
rownames(ileum_metadata) <- ileum_metadata$sampleID

## Liver
# metadata
metadata <- as.data.frame(readxl::read_xlsx(path = '/well/powrie/users/cqh100/work/mouse_rna_seq/metadata_fordeseq.xlsx'))

liver_metadata <- subset(metadata, tissue == 'liver')
liver_metadata <- subset(liver_metadata, age == 'Young')
liver_metadata <- subset(liver_metadata, group == 'Solo')
liver_metadata$status <- as.factor(liver_metadata$status)
rownames(liver_metadata) <- liver_metadata$sample_id
colnames(liver_metadata)[1] <- 'sampleID'
rownames(liver_metadata) <- liver_metadata$sampleID






# Ileum
ileum_counts_spf <- as.data.frame(ileum_normed_dat_corrected[, grep("^Ileum(1[0-5]|[1-9])$", colnames(ileum_normed_dat_corrected))])

## need to change the rownames to gene symbols 

ileum_counts_spf <- process_deseq_results(ileum_counts_spf, offline = TRUE)
ileum_counts_spf <- ileum_counts_spf$res[1:(length(ileum_counts_spf$res)-1)]
ileum_counts_spf_subset <- ileum_counts_spf[rownames(ileum_counts_spf) %in% rownames(results$Ileum_SPF_vs_GF$deseq), ]
ileum_counts_spf_subset <- ileum_counts_spf_subset[1:(length(ileum_counts_spf_subset)-1)]
ileum_counts_t_spf <- t(ileum_counts_spf_subset)


# Liver
liver_counts_spf <- as.data.frame(liver_normed_dat[, grep("^Liver(1[0-5]|[1-9])$", colnames(liver_normed_dat))])

## need to change the rownames to gene symbols 
liver_counts_spf <- process_deseq_results(liver_counts_spf, offline = TRUE)
liver_counts_spf <- liver_counts_spf$res[1:(length(liver_counts_spf$res)-1)]
liver_counts_spf_subset <- liver_counts_spf[rownames(liver_counts_spf) %in% rownames(results$Liver_SPF_vs_GF$deseq), ]
liver_counts_spf_subset <- liver_counts_spf_subset[1:(length(liver_counts_spf_subset)-1)]
liver_counts_t_spf <- t(liver_counts_spf_subset)



gsg_spf <- goodSamplesGenes(ileum_counts_t_spf)

if (!gsg_spf$allOK) {
  if (sum(!gsg_spf$goodGenes) > 0) 
    printFlush(paste("Removing genes:", paste(names(ileum_counts_t_spf)[!gsg_spf$goodGenes], collapse = ", ")))
  if (sum(!gsg_spf$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(ileum_counts_t_spf)[!gsg_spf$goodSamples], collapse = ", ")))
  
  ileum_counts_t_spf <- ileum_counts_t_spf[gsg_spf$goodSamples, gsg_spf$goodGenes]
}

spt_spf <- pickSoftThreshold(ileum_counts_t_spf)
softPower_spf <- 7
adjacency_spf <- adjacency(ileum_counts_t_spf, power = softPower_spf, type = 'signed')

TOM_spf <- TOMsimilarity(adjacency_spf)
TOM_spf.dissimilarity <- 1 - TOM_spf

geneTree_spf <- hclust(as.dist(TOM_spf.dissimilarity), method = "average")
Modules_spf <- cutreeDynamic(dendro = geneTree_spf, distM = TOM_spf.dissimilarity, deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = 100)
ModuleColors_spf <- labels2colors(Modules_spf)

MElist_spf <- moduleEigengenes(ileum_counts_t_spf, colors = ModuleColors_spf)
MEs_spf <- MElist_spf$eigengenes

ME_spf.dissimilarity <- 1 - cor(MEs_spf, use = "complete")
METree_spf <- hclust(as.dist(ME_spf.dissimilarity), method = "average")

merge_spf <- mergeCloseModules(ileum_counts_t_spf, ModuleColors_spf, cutHeight = 0.2)
mergedColors_spf <- merge_spf$colors
mergedMEs_spf <- merge_spf$newMEs

# create a simple mapping of module and gene
ileum_modules_df <- data.frame(gene = colnames(ileum_counts_t_spf),
                         module = mergedColors_spf)

# getkMEs for module membership - correlation with MEs
ileum_kMEs <- as.data.frame(WGCNA::cor(ileum_counts_t_spf, mergedMEs_spf, use = "p"))

# write out the module eigengenes
write.table(data.frame(sample = rownames(mergedMEs_spf), mergedMEs_spf),
            file = "export.dir/ileum_mes.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)



gsg_spf <- goodSamplesGenes(liver_counts_t_spf)

if (!gsg_spf$allOK) {
  if (sum(!gsg_spf$goodGenes) > 0) 
    printFlush(paste("Removing genes:", paste(names(liver_counts_t_spf)[!gsg_spf$goodGenes], collapse = ", ")))
  if (sum(!gsg_spf$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(liver_counts_t_spf)[!gsg_spf$goodSamples], collapse = ", ")))
  
  liver_counts_t_spf <- liver_counts_t_spf[gsg_spf$goodSamples, gsg_spf$goodGenes]
}

spt_spf_liver <- pickSoftThreshold(liver_counts_t_spf)
softPower_spf_liver <- 5
adjacency_spf_liver <- adjacency(liver_counts_t_spf, power = softPower_spf_liver, type = 'signed')

TOM_spf_liver <- TOMsimilarity(adjacency_spf_liver)
TOM_spf_liver.dissimilarity <- 1 - TOM_spf_liver

geneTree_spf_liver <- hclust(as.dist(TOM_spf_liver.dissimilarity), method = "average")
Modules_spf_liver <- cutreeDynamic(dendro = geneTree_spf_liver, distM = TOM_spf_liver.dissimilarity, deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = 200)
ModuleColors_spf_liver <- labels2colors(Modules_spf_liver)

MElist_spf_liver <- moduleEigengenes(liver_counts_t_spf, colors = ModuleColors_spf_liver)
MEs_spf_liver <- MElist_spf_liver$eigengenes

ME_spf_liver.dissimilarity <- 1 - cor(MEs_spf_liver, use = "complete")
METree_spf_liver <- hclust(as.dist(ME_spf_liver.dissimilarity), method = "average")

merge_spf_liver <- mergeCloseModules(liver_counts_t_spf, ModuleColors_spf_liver, cutHeight = 0.2)
mergedColors_spf_liver <- merge_spf_liver$colors
mergedMEs_spf_liver <- merge_spf_liver$newMEs

liver_modules_df <- data.frame(gene = colnames(liver_counts_t_spf),
                            module = mergedColors_spf_liver)

write.table(data.frame(sample = rownames(mergedMEs_spf_liver), mergedMEs_spf_liver),
            file = "export.dir/liver_mes.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


layout(matrix(1:6, nrow = 3, byrow = TRUE))
par(mar = c(5, 4, 4, 2))

# A
plot(spt_spf$fitIndices$Power, spt_spf$fitIndices$SFT.R.sq,
     xlab = "Soft Threshold", ylab = "R²",
     main = "Ileum: Scale-Free Topology Fit", type = "n",
     cex.lab = 2, cex.main = 2, cex.axis = 2)
text(spt_spf$fitIndices$Power, spt_spf$fitIndices$SFT.R.sq, labels = spt_spf$fitIndices$Power, col = "red", cex = 2)
abline(h = 0.70, col = "red")
#mtext("A", side = 3, adj = 0, line = 1.5, cex = 2.5, font = 2)

# B
plot(spt_spf_liver$fitIndices$Power, spt_spf_liver$fitIndices$SFT.R.sq,
     xlab = "Soft Threshold", ylab = "R²",
     main = "Liver: Scale-Free Topology Fit", type = "n",
     cex.lab = 2, cex.main = 2, cex.axis = 2)
text(spt_spf_liver$fitIndices$Power, spt_spf_liver$fitIndices$SFT.R.sq, labels = spt_spf_liver$fitIndices$Power, col = "red", cex = 2)
abline(h = 0.80, col = "red")
#mtext("B", side = 3, adj = 0, line = 1.5, cex = 2.5, font = 2)

# C
plot(spt_spf$fitIndices$Power, spt_spf$fitIndices$mean.k.,
     xlab = "Soft Threshold", ylab = "Connectivity",
     main = "Ileum: Mean Connectivity", type = "n",
     cex.lab = 2, cex.main = 2, cex.axis = 2)
text(spt_spf$fitIndices$Power, spt_spf$fitIndices$mean.k., labels = spt_spf$fitIndices$Power, col = "red", cex = 2)
#mtext("C", side = 3, adj = 0, line = 1.5, cex = 2.5, font = 2)

# D
plot(spt_spf_liver$fitIndices$Power, spt_spf_liver$fitIndices$mean.k.,
     xlab = "Soft Threshold", ylab = "Connectivity",
     main = "Liver: Mean Connectivity", type = "n",
     cex.lab = 2, cex.main = 2, cex.axis = 2)
text(spt_spf_liver$fitIndices$Power, spt_spf_liver$fitIndices$mean.k., labels = spt_spf_liver$fitIndices$Power, col = "red", cex = 2)
#mtext("D", side = 3, adj = 0, line = 1.5, cex = 2.5, font = 2)


### faceted version
combined_data <- bind_rows(
  spt_spf$fitIndices %>% mutate(tissue = "Ileum"),
  spt_spf_liver$fitIndices %>% mutate(tissue = "Liver")
) %>%
  select(Power, SFT.R.sq, mean.k., tissue) %>%
  pivot_longer(cols = c(SFT.R.sq, mean.k.), 
               names_to = "metric", 
               values_to = "value") %>%
  mutate(
    metric_label = case_when(
      metric == "SFT.R.sq" ~ "Scale-Free Topology Fit",
      metric == "mean.k." ~ "Mean Connectivity"
    ),
    y_label = case_when(
      metric == "SFT.R.sq" ~ "R²",
      metric == "mean.k." ~ "Connectivity"
    )
  )

ggplot(combined_data, aes(x = Power, y = value)) +
  geom_text(aes(label = Power), color = "red", size = 6) +
  geom_hline(data = filter(combined_data, metric == "SFT.R.sq", tissue == "Ileum"), 
             aes(yintercept = 0.70), color = "red") +
  geom_hline(data = filter(combined_data, metric == "SFT.R.sq", tissue == "Liver"), 
             aes(yintercept = 0.80), color = "red") +
  facet_grid(metric_label ~ tissue, scales = "free_y", 
             labeller = labeller(metric_label = label_value, tissue = label_value)) +
  labs(x = "Soft Threshold", y = "") +
  base_theme +
  theme(
    strip.text = element_text(size = 20, face = "bold")
  )



ileum_stats <- extract_wgcna_stats(ileum_counts_t_spf, mergedColors_spf, "ILE")
liver_stats <- extract_wgcna_stats(liver_counts_t_spf, mergedColors_spf_liver, "LIV")
knitr::kable(ileum_stats)
knitr::kable(liver_stats)


# For ileum data
ileum_results <- analyzeModuleData(mergedColors_spf, ileum_counts_t_spf, "Ileum")

# For liver data
liver_results <- analyzeModuleData(mergedColors_spf_liver, liver_counts_t_spf, "Liver")


### General module overview - number of gee

ileum_stats$Module <- factor(ileum_stats$Module, levels = ileum_stats$Module)
ggplot(ileum_stats, aes(x = Module, y = Number_of_Genes)) +
  geom_bar(stat="identity", fill = "#5E4FA2", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("export.dir/ileum_ngenes_modules.pdf", height = 2.5, width = 5)

liver_stats$Module <- factor(liver_stats$Module, levels = liver_stats$Module)
ggplot(liver_stats, aes(x = Module, y = Number_of_Genes)) +
  geom_bar(stat="identity", fill = "#9E0142", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("export.dir/liver_ngenes_modules.pdf", height = 2.5, width = 5)



### Plotting associated pathways

ileum_go <- ileum_results$go_terms
ileum_go <- ileum_go[order(ileum_go$Module),]
ileum_go$GOterm <- factor(ileum_go$GOterm,
                          levels = unique(ileum_go$GOterm))

module_colors <- as.character(unique(ileum_go$Module))
names(module_colors) <- unique(ileum_go$Module)

ggplot(ileum_go, aes(x = GOterm,
                     y = -log10(Pvalue),
                     fill = Module)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  scale_fill_manual(values = module_colors) +
  coord_flip() +
  facet_wrap(~Module, scales = "free", ncol = 1) +
  theme(legend.position = "none")

ggsave("export.dir/ileum_module_pathways.pdf", height = 20, width = 10)


# liver
liver_go <- liver_results$go_terms
liver_go <- liver_go[order(liver_go$Module),]
liver_go$GOterm <- factor(liver_go$GOterm,
                          levels = unique(liver_go$GOterm))

module_colors <- as.character(unique(liver_go$Module))
names(module_colors) <- unique(liver_go$Module)

ggplot(liver_go, aes(x = GOterm,
                     y = -log10(Pvalue),
                     fill = Module)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  scale_fill_manual(values = module_colors) +
  coord_flip() +
  facet_wrap(~Module, scales = "free", ncol = 1) +
  theme(legend.position = "none")

ggsave("export.dir/liver_module_pathways.pdf", height = 20, width = 10)



## Highest correlations in modules

ileum_cors <- list()
modules <- unique(ileum_modules_df$module)
for (i in 1:length(modules)){
  genes <- ileum_modules_df$gene[ileum_modules_df$module == modules[i]]
  all_cors <- cor(ileum_counts_t_spf[,genes])
  all_cors <- OCMSutility::sym_mat2df(all_cors)
  all_cors$module <- modules[i]
  all_cors <- all_cors[order(all_cors$value, decreasing = TRUE),]
  ileum_cors[[i]] <- all_cors
}
ileum_cors <- bind_rows(ileum_cors)

# liver
liver_cors <- list()
liver_modules <- unique(liver_modules_df$module)
for (i in 1:length(liver_modules)){
  genes <- liver_modules_df$gene[liver_modules_df$module == liver_modules[i]]
  all_cors <- cor(liver_counts_t_spf[,genes])
  all_cors <- OCMSutility::sym_mat2df(all_cors)
  all_cors$module <- liver_modules[i]
  all_cors <- all_cors[order(all_cors$value, decreasing = TRUE),]
  liver_cors[[i]] <- all_cors
}
liver_cors <- bind_rows(liver_cors)




# ileum
ileum_median_cors <- list()
for (i in modules){
  median_cor <- median(ileum_cors$value[ileum_cors$X1 != ileum_cors$X2 & ileum_cors$module == i])
  ileum_median_cors[[i]] <- median_cor
}

liver_median_cors <- list()
for (i in liver_modules){
  median_cor <- median(liver_cors$value[liver_cors$X1 != liver_cors$X2 & liver_cors$module == i])
  liver_median_cors[[i]] <- median_cor
}

## PC1 instead ie. variance explained for each module

# Ileum
ileum_pc1 <- list()
for (i in modules){
  genes <- ileum_modules_df$gene[ileum_modules_df$module == i]
  pc <- prcomp(ileum_counts_t_spf[,genes])
  var_explained <- pc$sdev^2
  pve_pc1 <- var_explained[1] / sum(var_explained)
  ileum_pc1[[i]] <- pve_pc1
}

# Liver
liver_pc1 <- list()
for (i in liver_modules){
  genes <- liver_modules_df$gene[liver_modules_df$module == i]
  pc <- prcomp(liver_counts_t_spf[,genes])
  var_explained <- pc$sdev^2
  pve_pc1 <- var_explained[1] / sum(var_explained)
  liver_pc1[[i]] <- pve_pc1
}

# add to gene set size
ileum_stats$Module <- as.character(ileum_stats$Module)
ileum_stats$module <- gsub(".*- ", "", ileum_stats$Module)
ileum_stats$PC1_VE <- unlist(ileum_pc1[ileum_stats$module], use.names = FALSE)

liver_stats$module <- gsub(".*- ", "", liver_stats$Module)
liver_stats$PC1_VE <- unlist(liver_pc1[liver_stats$module], use.names = FALSE)

ileum_stats$Module <- factor(ileum_stats$Module, levels = ileum_stats$Module)
ggplot(ileum_stats, aes(x = Module, y = Number_of_Genes, size = PC1_VE)) +
  geom_point(color = "#5E4FA2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("export.dir/ileum_ngenes_connectivity_modules.pdf", height = 2.5, width = 5)


liver_stats$Module <- factor(liver_stats$Module, levels = liver_stats$Module)
ggplot(liver_stats, aes(x = Module, y = Number_of_Genes, size = PC1_VE)) +
  geom_point(color = "#9E0142") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("export.dir/liver_ngenes_connectivity_modules.pdf", height = 2.5, width = 5)




## Preservation

allowWGCNAThreads(nThreads = 14)

if (!file.exists("export.dir/ileum_preservation.Rds")){

# ==============================================================================
# DATA PREPROCESSING
# ==============================================================================

ileum_counts_all <- as.data.frame(ileum_normed_dat_corrected)
ileum_counts_all <- process_deseq_results(ileum_counts_all, offline = TRUE)
ileum_counts_all <- ileum_counts_all$res[1:(length(ileum_counts_all$res)-1)]
ileum_counts_all <- ileum_counts_all[1:(length(ileum_counts_all)-1)]
ileum_counts_all <- ileum_counts_all[rowSums(ileum_counts_all >= 10) >= 15, ]

# ==============================================================================
# SAMPLE EXTRACTION
# ==============================================================================

spf_samples_ileum <- ileum_metadata$sampleID[ileum_metadata$age == 'Young' & 
                                       ileum_metadata$group == 'Solo' &
                                       ileum_metadata$status == "SPF"]

gf_samples_ileum <- ileum_metadata$sampleID[ileum_metadata$age == 'Young' & 
                                      ileum_metadata$group == 'Solo' &
                                      ileum_metadata$status == "GF"]

mm12_samples_ileum <- ileum_metadata$sampleID[ileum_metadata$age == 'Young' & 
                                        ileum_metadata$group == 'Solo' &
                                        ileum_metadata$status == "MM12"]

extract_group_data <- function(expr_data, sample_ids) {
  return(expr_data[, sample_ids])
}

spf_data_ileum  <- extract_group_data(ileum_counts_all, spf_samples_ileum)
gf_data_ileum  <- extract_group_data(ileum_counts_all, gf_samples_ileum)
mm12_data_ileum  <- extract_group_data(ileum_counts_all, mm12_samples_ileum)


# ==============================================================================
# MODULE PRESERVATION ANALYSIS
# ==============================================================================

# Test preservation in GF
multiExpr_gf_ileum <- list(
  SPF = list(data = ileum_counts_t_spf),
  GF = list(data = t(gf_data_ileum))
)
multiColor_spf <- list(SPF = mergedColors_spf)

mp_spf_to_gf_ileum <- modulePreservation(
  multiExpr_gf_ileum,
  multiColor_spf,
  referenceNetworks = 1,
  networkType = "signed",
  nPermutations = 500,
  randomSeed = 1,
  parallelCalculation = TRUE,
  verbose = 3
)

# Test preservation in MM12
multiExpr_mm12_ileum <- list(
  SPF = list(data = ileum_counts_t_spf),
  MM12 = list(data = t(mm12_data_ileum))
)

mp_spf_to_mm12_ileum <- modulePreservation(
  multiExpr_mm12_ileum,
  multiColor_spf,
  referenceNetworks = 1,
  networkType = "signed",
  nPermutations = 500,
  parallelCalculation = TRUE,
  randomSeed = 1,
  verbose = 3
)
# ==============================================================================
# EXTRACT PRESERVATION STATISTICS
# ==============================================================================

stats_gf_ileum <- mp_spf_to_gf_ileum$preservation$Z$ref.SPF$inColumnsAlsoPresentIn.GF
stats_mm12_ileum <- mp_spf_to_mm12_ileum$preservation$Z$ref.SPF$inColumnsAlsoPresentIn.MM12

# Create comprehensive preservation comparison dataframe
preservation_comparison_ileum <- data.frame(
  Module = rownames(stats_gf_ileum),
  ModuleSize = stats_gf_ileum[, "moduleSize"],
  Z_GF = stats_gf_ileum[, "Zsummary.pres"],
  Z_MM12 = stats_mm12_ileum[, "Zsummary.pres"],
  GF_Coexpression = stats_gf_ileum[, "Zdensity.pres"],
  GF_HubStability = stats_gf_ileum[, "Zconnectivity.pres"],
  MM12_Coexpression = stats_mm12_ileum[, "Zdensity.pres"],
  MM12_HubStability = stats_mm12_ileum[, "Zconnectivity.pres"]
)

# Add GO terms
module_descriptions_ileum <- setNames(ileum_results$module_terms$GOterm, ileum_results$module_terms$Module)
preservation_comparison_ileum$GOterm <- module_descriptions_ileum[preservation_comparison_ileum$Module]
preservation_comparison_ileum <- preservation_comparison_ileum[!is.na(preservation_comparison_ileum$GOterm), ]

saveRDS(preservation_comparison_ileum, file = "export.dir/ileum_preservation.Rds")

# ==============================================================================
# DATA PREPROCESSING
# ==============================================================================

liver_counts_all <- as.data.frame(liver_normed_dat)
liver_counts_all <- process_deseq_results(liver_counts_all, offline = TRUE)
liver_counts_all <- liver_counts_all$res[1:(length(liver_counts_all$res)-1)]
liver_counts_all <- liver_counts_all[1:(length(liver_counts_all)-1)]
liver_counts_all <- liver_counts_all[rowSums(liver_counts_all >= 10) >= 15, ]

# ==============================================================================
# SAMPLE EXTRACTION
# ==============================================================================

spf_samples_liver <- liver_metadata$sampleID[liver_metadata$age == 'Young' & 
                                       liver_metadata$group == 'Solo' &
                                       liver_metadata$status == "SPF"]

gf_samples_liver <- liver_metadata$sampleID[liver_metadata$age == 'Young' & 
                                      liver_metadata$group == 'Solo' &
                                      liver_metadata$status == "GF"]

mm12_samples_liver <- liver_metadata$sampleID[liver_metadata$age == 'Young' & 
                                        liver_metadata$group == 'Solo' &
                                        liver_metadata$status == "MM12"]

spf_data_liver <- extract_group_data(liver_counts_all, spf_samples_liver)
gf_data_liver <- extract_group_data(liver_counts_all, gf_samples_liver)
mm12_data_liver <- extract_group_data(liver_counts_all, mm12_samples_liver)


# ==============================================================================
# MODULE PRESERVATION ANALYSIS
# ==============================================================================

# Test preservation in GF
multiExpr_gf_liver <- list(
  SPF = list(data = liver_counts_t_spf),
  GF = list(data = t(gf_data_liver))
)
multiColor_spf <- list(SPF = mergedColors_spf_liver)

mp_spf_to_gf_liver <- modulePreservation(
  multiExpr_gf_liver,
  multiColor_spf,
  referenceNetworks = 1,
  networkType = "signed",
  nPermutations = 500,
  randomSeed = 1,
  parallelCalculation = TRUE,
  verbose = 3
)

# Test preservation in MM12
multiExpr_mm12_liver <- list(
  SPF = list(data = liver_counts_t_spf),
  MM12 = list(data = t(mm12_data_liver))
)

mp_spf_to_mm12_liver <- modulePreservation(
  multiExpr_mm12_liver,
  multiColor_spf,
  referenceNetworks = 1,
  networkType = "signed",
  nPermutations = 500,
  randomSeed = 1,
  parallelCalculation = TRUE,
  verbose = 3
)

# ==============================================================================
# EXTRACT PRESERVATION STATISTICS
# ==============================================================================

stats_gf_liver <- mp_spf_to_gf_liver$preservation$Z$ref.SPF$inColumnsAlsoPresentIn.GF
stats_mm12_liver <- mp_spf_to_mm12_liver$preservation$Z$ref.SPF$inColumnsAlsoPresentIn.MM12

# Create comprehensive preservation comparison dataframe
preservation_comparison_liver <- data.frame(
  Module = rownames(stats_gf_liver),
  ModuleSize = stats_gf_liver[, "moduleSize"],
  Z_GF = stats_gf_liver[, "Zsummary.pres"],
  Z_MM12 = stats_mm12_liver[, "Zsummary.pres"],
  GF_Coexpression = stats_gf_liver[, "Zdensity.pres"],
  GF_HubStability = stats_gf_liver[, "Zconnectivity.pres"],
  MM12_Coexpression = stats_mm12_liver[, "Zdensity.pres"],
  MM12_HubStability = stats_mm12_liver[, "Zconnectivity.pres"]
)

# Add GO terms
module_descriptions_liver <- setNames(liver_results$module_terms$GOterm, liver_results$module_terms$Module)
preservation_comparison_liver$GOterm <- module_descriptions_liver[preservation_comparison_liver$Module]
preservation_comparison_liver <- preservation_comparison_liver[!is.na(preservation_comparison_liver$GOterm), ]
saveRDS(preservation_comparison_liver, file = "export.dir/liver_preservation.Rds")

}else{
  preservation_comparison_ileum <- readRDS("export.dir/ileum_preservation.Rds")
  preservation_comparison_liver <- readRDS("export.dir/ileum_preservation.Rds")}
