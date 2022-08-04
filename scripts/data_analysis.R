################################################################################
# Donagh Egan
# POIAZ 
# Date: 03 August 2022
# 

# Description: Initial analysis of DIA mass spec data for IP proteomics. See email:
# Dia-nn output (PDCD1), 27 Jul 2022
################################################################################

## Library #### 
################################################################################

library(readxl)
library(readr)
library(tidyverse)
library(grid)
library(ggsci)
library(PhosR)

## Read in data ####
################################################################################

exp_design <- read_xlsx("/home/degan/dia_mass_spec/inputs/2022 07 14 exp85 exp details.xlsx", sheet = 2,skip = 3)
report_mat <- read_tsv("/home/degan/dia_mass_spec/inputs/report.gg_matrix.tsv")
report_mat <- report_mat %>% column_to_rownames("Genes")
colnames(report_mat) <-  substr(colnames(report_mat), 133, 147)

## formulating exp design ####
################################################################################

annotation_df <- data.frame(row.names = colnames(report_mat),
                            sample_now = rep(exp_design$no, each = 2),
                            biological_rep = rep(exp_design$...4, each = 2),
                            condition = rep(exp_design$description, each = 2),
                            type = rep(exp_design$...3, each = 2))


## Heatmap of missing values ####
################################################################################

missval <- report_mat[apply(report_mat, 1, function(x) any(is.na(x))), ] # select only proteins with NA values
missval <- ifelse(is.na(missval), 0, 1)

col_fun = pal_npg("nrc")(10)[1:6]
names(col_fun) <- na.omit(unique(annotation_df$condition))

heat_annotation <- ComplexHeatmap::HeatmapAnnotation(condition = annotation_df$condition,
                                                     show_annotation_name = FALSE,
                                                     col = list(condition = col_fun))

pdf("/home/degan/dia_mass_spec/figures/missingval_pd1.pdf", height = 5, width = 10)
ComplexHeatmap::Heatmap(missval, col = c("white", "black"), column_names_side = "top", 
                        show_row_names = FALSE,  top_annotation = heat_annotation, show_column_names = T, name = "Missing values pattern", 
                        column_names_gp = gpar(fontsize = 7), show_column_dend = T, cluster_columns = T, cluster_rows = T, 
                        heatmap_legend_param = list(at = c(0,1), labels = c("Missing value", "Valid value")))
dev.off()

## visualize number of proteins/sample ####
################################################################################

protein_data_number <- ifelse(is.na(report_mat), 0, 1)
protein_data_number <- as.data.frame(colSums(protein_data_number))
protein_data_number$experiment <- rownames(protein_data_number)
keep <- protein_data_number[apply(protein_data_number, MARGIN = 1, function(x) !all(is.na(x) == TRUE)), ]
protein_data_number <- cbind(protein_data_number, annotation_df)
protein_data_number <- protein_data_number %>% arrange(desc(abs(`colSums(protein_data_number)`)))

pdf("/home/degan/dia_mass_spec/figures/proteins_persample.pdf", width = 6, height = 5)
ggplot(protein_data_number, aes(x = experiment, y = `colSums(protein_data_number)`, fill = condition)) + geom_bar(stat = "identity", width = 0.7) + 
  geom_hline(yintercept = nrow(keep)) + labs(title = "Proteins per sample", 
                                             x = "", y = "Number of proteins")  + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=4)) + scale_fill_npg()
dev.off()

## Filter proteins ####
################################################################################

protein_data_filtered <- selectGrps(as.matrix(report_mat), annotation_df$condition, 5/6, n=1)

## normalize filtered protein data using global vsn transformation ####
################################################################################

protein_data_norm <- log2(protein_data_filtered) 
protein_data_norm <- tImpute(protein_data_norm)

## PCA ####
################################################################################

## Running PCA ####
pca <- prcomp(t(protein_data_norm), center = T, scale. = T)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

PCAvalues <- data.frame(pca$x)
PCAloadings <- data.frame(Variables = rownames(pca$rotation),  pca$rotation)

PCAvalues <- merge(PCAvalues, annotation_df, by=0)
PCAvalues <- PCAvalues %>% column_to_rownames("Row.names")

pdf("/home/degan/dia_mass_spec/figures/PCA_condition.pdf", height = 3, width = 4)
ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = condition, shape = type)) +
  geom_point(size = 1) + scale_color_npg() +
  theme_bw() +
  ggplot2::labs(
    x = paste0("PC1: ", round(percentVar[1] *  100), "% variance"),
    y = paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
dev.off()

## CV per antibody ####
annotation_df$type_cond <- paste(annotation_df$type, annotation_df$condition, sep = "_")
  
proteins_cv <- merge(t(protein_data_norm), annotation_df[, "type_cond", drop = FALSE], by=0)
proteins_cv <- proteins_cv %>% column_to_rownames("Row.names")

proteins_cv <- aggregate(.~ type_cond, data = proteins_cv, 
                         function(x) sd(x, na.rm = T) / mean(x, na.rm = T) * 100)

proteins_cv_melt <- reshape2::melt(proteins_cv)

pdf("/home/degan/dia_mass_spec/figures/coef_variation_condition.pdf", height = 5, width = 4)
ggplot(data = proteins_cv_melt, aes(x=reorder(type_cond, value), y=value, color = type_cond)) + geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 8)) + scale_color_npg() + labs(color = "Condition") + xlab("")
dev.off()

## within condition correlation ####
################################################################################

cor_protein_data <- cor(protein_data_norm)
cor_protein_data <- merge(cor_protein_data, annotation_df, by=0)
cor_protein_data <- column_to_rownames(cor_protein_data, "Row.names")

## subset correlations according to conditions ####
################################################################################
list_conditions <- list()
for (i in unique(cor_protein_data$type_cond)) {
  df_temp <- cor_protein_data[cor_protein_data$type_cond == i, , drop=F]
  df_temp <- df_temp[, which(colnames(df_temp) %in% rownames(df_temp))]
  diag(df_temp) <- NA
  list_conditions[[i]] <- df_temp
  
}

## find mean corrlation/conditon and add to list 
################################################################################
cor_values <- list()
for (i in names(list_conditions)) {
  values <- list_conditions[[i]][,1]
  values <- mean(values, na.rm=T)
  sd <- sd(values)
  cor_values[[i]] <- values
}

## Format and plot results ####
################################################################################

cor_df <- data.frame(avg_cor = matrix(unlist(cor_values), nrow=length(cor_values), byrow=TRUE),
                     row.names = names(cor_values))
cor_df$condition <- rownames(cor_df)

pdf("/home/degan/dia_mass_spec/figures/condition_correlation.pdf",height = 4, width = 3)
ggplot(data=cor_df, aes(x=reorder(condition, desc(avg_cor)), y=avg_cor, fill = condition)) +
  geom_bar(stat="identity", width = 0.5, fill=pal_npg("nrc")(10)[1]) +  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  geom_hline(aes(yintercept = mean(avg_cor)), linetype = "dashed") + ylim(0,1) +
  ylab("Average PCC")+ xlab("Condition") + theme(legend.position = "none" )
dev.off()

