"""
File: data_processing.py
Author: Oscar Laurent
Date: 2024-09-11
Description: A R script to perform analysis of the processed csv data that came from the "data/ToF-SIMS_data/processed_sepctra" folder
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
#  - CLUSTER RATIO WITH ALL ELEMENTS
# =======================================================

str_all_primary_ion <- "^(Pt|Pd|Ru|Ir|Rh){1}$"

scaled_automatic_positive_df %>% select((matches("^(Ru\\d?|Pd\\d?|Pt\\d?){2,}$")))  %>% colnames()

sample_names <- c("Ru" = "Ru", "Ru-Pt" = "Ru-Pt", "Ru-Pd" = "Ru-Pd",
                  "Ru-Ir" = "Ru-Ir","Ru-Rh" = "Ru-Rh",
                  "Ru-Pt-Pd" = "Ru-Pt-Pd", 
                  "Ru-Pt-Pd-Ir" = "Ru-Pt-Pd-Ir", 
                  "Ru-Pt-Pd-Ir-Rh" = "HEA Mix", 
                  "HEA19" = "HEA Seg")

anova_df <- scaled_automatic_positive_df %>% 
  subset(sample %in% c("Ru", "Ru-Pt", "Ru-Pd", "Ru-Ir", "Ru-Rh",
                      "Ru-Pt-Pd", "Ru-Pt-Pd-Ir", "Ru-Pt-Pd-Ir-Rh", "HEA19")) %>% 
  mutate(HACS = rowSums(across(matches("^(Ru\\d?|Pd\\d?|Pt\\d?|Rh\\d?|Ir\\d?){2,}$"))),
         IACS = rowSums(across(matches("^(Pd\\d|Pt\\d|Ru\\d|Rh\\d|Ir\\d){1}$")))) %>% 
  select(sample, primary_ion, replicat, HACS, IACS, matches(str_all_primary_ion)) %>% 
  mutate(HACR = HACS/rowSums(across(matches(str_all_primary_ion))), 
         IACR = IACS/rowSums(across(matches(str_all_primary_ion)))) %>%
  mutate(CR = HACR/IACR) 

# Extracting the CR for each sample in a Df
CR_df <- anova_df %>% select(sample, primary_ion, replicat, CR)
  
# Selection of Bi5 and computing of the SE intra sample
  anova_df <- anova_df %>% 
  select(!replicat) %>% 
  subset(primary_ion == "Bi5") %>% 
  group_by(sample) %>% 
  mutate(SE_CR = sd(CR)/sqrt(n()), mean_CR = mean(CR)) %>% 
  select(sample, primary_ion, CR, SE_CR, mean_CR) 
  #pivot_longer(!c(sample, primary_ion)) 

p <- anova_df %>% 
  ggplot(aes(x = 1,  y = CR, color = sample)) + 
  geom_point(size = 6) + 
  scale_color_manual(values = c("#e55e57", "#e5a273", "#ffd479", "#3bc093", "#07b5bc" ,
                                 "#509bdd", "#8d85e5", "#b373d0", "#D452AA")) +
  stat_summary(geom = "point",fun = "mean", size = 7, fill = "red", color = "red",
              shape = 23) +
  facet_grid(~sample, switch = "x", scales = "free", labeller = labeller(sample = sample_names)) +
  scale_x_continuous(expand = c(0, 0.55))+
  geom_errorbar(aes(group = sample, ymin = mean_CR- 1.96*SE_CR, ymax = mean_CR + 1.96*SE_CR),
                color = "black", width = 0.3) +
  ylab("Cluster Ratio") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=30, color = "black"),
        axis.title.y=element_text(size=30),
        axis.title.x= element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size=25,face="bold"),
        panel.border = element_rect(color = "#666666"),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.background = element_rect(fill='transparent'))

# Define the color palette for the samples
sample_colors <- c("Ru" = "#e55e57" ,"Ru-Pt" = "#e5a273", "Ru-Pd" = "#ffd479",
                   "Ru-Ir" = "#3bc093", "Ru-Rh" = "#07b5bc", 
                   "Ru-Pt-Pd" = "#509bdd", "Ru-Pt-Pd-Ir" = "#8d85e5", 
                   "Ru-Pt-Pd-Ir-Rh" = "#b373d0", "HEA19" = "#D452AA")

# Function to get the strip background colors
strip_background_fill <- function(p) {
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-b', g$layout$name))
  fills <- sample_colors
  
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  return(g)
}

g <- strip_background_fill(p)

# Draw the plot
grid::grid.draw(g)

# Perform pairwise t-tests
pairwise.t.test(anova_df$CR, anova_df$sample, p.adjust.method = "bonferroni")



# =======================================================
#  - CLUSTER RATIO WITH ONLY Pt AND Pd
# =======================================================

str_primary_ion <- "^(Pt|Pd){1}$"

scaled_automatic_positive_df %>% select((matches("^(Pd\\d?|Pt\\d?){2,}$")))  %>% colnames()

sample_names <- c("Ru-Pt-Pd" = "Ru-Pt-Pd", 
                  "Ru-Pt-Pd-Ir" = "Ru-Pt-Pd-Ir", 
                  "Ru-Pt-Pd-Ir-Rh" = "HEA Mix", 
                  "HEA19" = "HEA Seg")

anova_df <- scaled_automatic_positive_df %>% 
  subset(sample %in% c("Ru-Pt-Pd", "Ru-Pt-Pd-Ir", "Ru-Pt-Pd-Ir-Rh", "HEA19")) %>% 
  mutate(HACS = rowSums(across(matches("^(Pd\\d?|Pt\\d?){2,}$"))),
         IACS = rowSums(across(matches("^(Pd\\d|Pt\\d){1}$")))) %>% 
  select(sample, primary_ion, replicat, HACS, IACS, matches(str_primary_ion)) %>% 
  mutate(HACR = HACS/rowSums(across(matches(str_primary_ion))), 
         IACR = IACS/rowSums(across(matches(str_primary_ion)))) %>%
  mutate(CR = HACR/IACR) 

# Extracting the CR for each sample in a Df
CR_df <- anova_df %>% select(sample, primary_ion, replicat, CR)
  
# Selection of Bi5 and computing of the SE intra sample
  anova_df <- anova_df %>% 
  select(!replicat) %>% 
  subset(primary_ion == "Bi5") %>% 
  group_by(sample) %>% 
  mutate(SE_CR = sd(CR)/sqrt(n()), mean_CR = mean(CR)) %>% 
  select(sample, primary_ion, CR, SE_CR, mean_CR) 

p <- anova_df %>% 
  ggplot(aes(x = 1,  y = CR, color = sample)) + 
  geom_point(size = 6) + 
  scale_color_manual(values = c("#509bdd", "#8d85e5", "#b373d0", "#D452AA")) +
  stat_summary(geom = "point",fun = "mean", size = 7, fill = "red", color = "red",
              shape = 23) +
  facet_grid(~sample, switch = "x", scales = "free", labeller = labeller(sample = sample_names)) +
  scale_x_continuous(expand = c(0, 0.55))+
  geom_errorbar(aes(group = sample, ymin = mean_CR- 1.96*SE_CR, ymax = mean_CR + 1.96*SE_CR),
                color = "black", width = 0.3) +
  ylab("Cluster Ratio") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=30, color = "black"),
        axis.title.y=element_text(size=30),
        axis.title.x= element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size=25,face="bold"),
        panel.border = element_rect(color = "#666666"),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.background = element_rect(fill='transparent'))

# Define the color palette for the samples
sample_colors <- c("Ru-Pt-Pd" = "#509bdd", "Ru-Pt-Pd-Ir" = "#8d85e5", 
                   "Ru-Pt-Pd-Ir-Rh" = "#b373d0", "HEA19" = "#D452AA")

# Function to get the strip background colors
strip_background_fill <- function(p) {
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-b', g$layout$name))
  fills <- sample_colors
  
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  return(g)
}

g <- strip_background_fill(p)

# Draw the plot
grid::grid.draw(g)

# Perfomr pairwise t-tests
pairwise.t.test(anova_df$CR, anova_df$sample, p.adjust.method = "bonferroni")
