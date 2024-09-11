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
