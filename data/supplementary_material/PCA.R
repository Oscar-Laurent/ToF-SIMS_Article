"""
File: data_processing.py
Author: Oscar Laurent
Date: 2024-09-11
Description: A R script to perform PCA  of the processed csv data that came from the "data/ToF-SIMS_data/processed_sepctra" folder
"""

library(tidyverse)
library(readr)
library(FactoMineR)
library(factoextra)
library(stringr)
library(ggrepel)
library(ggpubr)
library(rstatix)
library(upstartr)
library(latex2exp)


theme_set(theme_bw() + 
            theme(axis.text=element_text(size=12, colour="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=13,face="bold"),
        strip.text.x = element_text(size=13,face="bold")))



# =======================================================
#  - DATA RETRIEVAL
# =======================================================

automatic_positive_df <- read.csv("data/ToF-SIMS_data/processed_spectra/AllSample_PositiveMode_ProcessedData.csv")
automatic_positive_df$replicat <- as.factor(automatic_positive_df$replicat)
automatic_positive_df$sample <- factor(automatic_positive_df$sample_name,
                                       levels = c("Ru", "Ru-Pt", "Ru-Pd", "Ru-Ir",
                                                  "Ru-Rh", "Ru-Pt-Pd", "Ru-Pt-Pd-Ir",
                                                  "Ru-Pt-Pd-Ir-Rh", "HEA19", "HEA_1.2"))
metadata_positive <- automatic_positive_df %>% select(!where(is.numeric))

# Scaling the spectra by their total counts
automatic_positive_df <- automatic_positive_df %>% 
  mutate("total_count" = rowSums(across(where(is.numeric))))

scaled_automatic_positive_df <- automatic_positive_df %>% mutate(across(where(is.numeric) &!total_count,  ~ .x/total_count))

automatic_negative_df <- read.csv("data/ToF-SIMS_data/processed_spectra/AllSample_NegativeMode_ProcessedData.csv")
automatic_negative_df$replicat <- as.factor(automatic_negative_df$replicat)
automatic_negative_df$sample <- factor(automatic_negative_df$sample_name,
                                       levels = c("Ru", "Ru-Pt", "Ru-Pd", "Ru-Ir",
                                                  "Ru-Rh", "Ru-Pt-Pd", "Ru-Pt-Pd-Ir",
                                                  "Ru-Pt-Pd-Ir-Rh", "HEA19", "HEA_1.2"))
metadata_negative <- automatic_negative_df %>% select(!where(is.numeric))

# Scaling the spectra by their total counts
automatic_negative_df <- automatic_negative_df %>% 
  mutate("total_count" = rowSums(across(where(is.numeric))))

scaled_automatic_negative_df <- automatic_negative_df %>% mutate(across(where(is.numeric) &!total_count,  ~ .x/total_count))



# =======================================================
#  - PCA ANALYSIS 
# =======================================================

PCA_df <- scaled_automatic_positive_df %>% 
  subset(sample == "Ru-Pt-Pd-Ir") %>% select(where(~ !all(. == 0)))

# Function to categorize variables
categorize_variable <- function(var_name) {
  num_elements <- str_count(var_name, pattern = "(Ru|Pd|Pt|Ir|Rh)")  # Count the number of elements
  if (num_elements == 4) {
    return("4 Metal Cluster")
  } else if (num_elements == 3) {
    return("3 Metal Cluster")
  } else if (num_elements == 2) {
    return("2 Metal Cluster")
  } else if (num_elements == 1) {
    return("1 Metal Cluster")
  } else {
    return("Other")
  }
}

# Apply categorization to the variable names
variable_categories <- PCA_df %>% select(where(is.numeric) & !total_count) %>% 
  colnames() %>% 
  map_chr(categorize_variable)

PCA_result <- PCA_df %>% 
  select(where(is.numeric) & !total_count) %>% 
  PCA(scale.unit = T, ncp = 4, graph = F)

fviz_pca_ind(PCA_result,
             axes = c(1, 2),
             geom.ind = "point", # Show points only (no text)
             col.ind = PCA_df$primary_ion, # Color by sample
             palette = "jco", 
             addEllipses = TRUE, 
             legend.title = "Sample")

# Extract variable loadings and contributions
var_loadings <- as.data.frame(PCA_result$var$coord)
var_loadings$variable <- rownames(var_loadings)
var_loadings$contrib <- PCA_result$var$contrib[, 1] + PCA_result$var$contrib[, 2]  # Contribution of PC1 + PC2
var_loadings$category <- variable_categories
var_loadings$cluster <- rownames(var_loadings)
var_loadings$emph <- var_loadings$Dim.1 < -0.55
var_loadings$emph <- var_loadings$emph | var_loadings$cluster %in% c("Ru2PtPd3","Ru3PtPd3", "Pt2", "Pt3", "Pt2Pd", "Pt3Pd", "RuPtPd3Ir3")

# Extract individual scores
ind_scores <- as.data.frame(PCA_result$ind$coord)
ind_scores$primary_ion <- PCA_df$primary_ion 

# Plot with ggplot2
ggplot(var_loadings, aes(x = Dim.1, y = Dim.2, color = category, label =cluster )) +
  geom_point(size = 3) +
  geom_text_repel(data = subset(var_loadings, emph),
                  size = 7,
                  max.overlaps = 10,
                  box.padding = 0.3,
                  force = 50,
                  seed = 1,
                  direction = "both",
                  aes(color = category)) +
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 3 , bycol = TRUE)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  scale_color_manual(values = c("4 Metal Cluster" = "#005f73",
                                "3 Metal Cluster" = "#00AFBB",
                                "2 Metal Cluster" = "#E7B800",
                                "1 Metal Cluster" = "#FC4E07",
                                "Other" = "grey")) +
  labs(x = paste0("PC 1 (", round(PCA_result$eig[1, 2], 1), "%)"),
       y = paste0("PC 2 (", round(PCA_result$eig[2, 2], 1), "%)"),
       color = "Category",) +
  theme_bw() +
  theme(legend.position = "top", 
        panel.grid = element_blank(),
        legend.text = element_text(size=30),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, 'cm'),
        axis.title = element_text(size=25), 
        axis.text = element_text(size=20),
        panel.border = element_rect(color = "#666666")) 

# Plot individual scores with ggplot2
ggplot(ind_scores, aes(x = Dim.1, y = Dim.2, color = primary_ion)) +
  geom_point(size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  #scale_color_manual(values = c("sample1_color", "sample2_color", "sample3_color")) +
  labs(x = paste0("PC 1 (", round(PCA_result$eig[1, 2], 1), "%)"),
       y = paste0("PC 2 (", round(PCA_result$eig[2], 2), "%)"),
       color = "Primary Ions") +
  theme_bw()+
  theme(legend.position = "top", 
        panel.grid = element_blank(),
        legend.text = element_text(size=30),
        legend.title = element_blank(),
        axis.title = element_text(size=25), 
        axis.text = element_text(size=20),
        panel.border = element_rect(color = "#666666"))
