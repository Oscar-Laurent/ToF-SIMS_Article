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
  mutate(CR = HACS/IACS) 

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

sample_names <- c("Ru-Pt-Pd" = "Ru-Pt-Pd", 
                  "Ru-Pt-Pd-Ir" = "Ru-Pt-Pd-Ir", 
                  "Ru-Pt-Pd-Ir-Rh" = "HEA Mix", 
                  "HEA19" = "HEA Seg")

anova_df <- scaled_automatic_positive_df %>% 
  subset(sample %in% c("Ru-Pt-Pd", "Ru-Pt-Pd-Ir", "Ru-Pt-Pd-Ir-Rh", "HEA19")) %>% 
  mutate(HACS = rowSums(across(matches("^(Pd\\d?|Pt\\d?){2,}$"))),
         IACS = rowSums(across(matches("^(Pd\\d|Pt\\d){1}$")))) %>% 
  select(sample, primary_ion, replicat, HACS, IACS) %>% 
  mutate(CR = HACS/IACS) 

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



# =======================================================
#  - CLUSTER RATIO WITH ONLY Pt Pd AND Ir
# =======================================================

sample_names <- c("Ru-Pt-Pd-Ir" = "Ru-Pt-Pd-Ir", 
                  "Ru-Pt-Pd-Ir-Rh" = "HEA Mix", 
                  "HEA19" = "HEA Seg")

anova_df <- scaled_automatic_positive_df %>% 
  subset(sample %in% c("Ru-Pt-Pd-Ir", "Ru-Pt-Pd-Ir-Rh", "HEA19")) %>% 
  mutate(HACS = rowSums(across(matches("^(Pd\\d?|Pt\\d?|Ir\\d?){2,}$"))),
         IACS = rowSums(across(matches("^(Pd\\d|Pt\\d|Ir\\d){1}$")))) %>% 
  select(sample, primary_ion, replicat, HACS, IACS) %>%
  mutate(CR = HACS/IACS) 

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
  scale_color_manual(values = c("#8d85e5", "#b373d0", "#D452AA")) +
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
sample_colors <- c("Ru-Pt-Pd-Ir" = "#8d85e5", 
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
#  - CLUSTER RATIO WITH ONLY Pt Pd Ir AND Rh
# =======================================================

sample_names <- c("Ru-Pt-Pd-Ir" = "Ru-Pt-Pd-Ir", 
                  "Ru-Pt-Pd-Ir-Rh" = "HEA Mix", 
                  "HEA19" = "HEA Seg")

anova_df <- scaled_automatic_positive_df %>% 
  subset(sample %in% c("Ru-Pt-Pd-Ir-Rh", "HEA19")) %>% 
  mutate(HACS = rowSums(across(matches("^(Pd\\d?|Pt\\d?|Ir\\d?|Rh\\d?){2,}$"))),
         IACS = rowSums(across(matches("^(Pd\\d|Pt\\d|Ir\\d|Rh\\d){1}$")))) %>% 
  select(sample, primary_ion, replicat, HACS, IACS) %>% 
  mutate(CR = HACS/IACS) 

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
  scale_color_manual(values = c("#b373d0", "#D452AA")) +
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
sample_colors <- c("Ru-Pt-Pd-Ir-Rh" = "#b373d0", "HEA19" = "#D452AA")

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
#  - OXIDATION RATIO WITH ALL THE ELEMENTS
# =======================================================

str_all_primary_ion <- "^(Pt|Pd|Ru|Ir|Rh){1}$"
oxide_df <- scaled_automatic_negative_df %>% mutate(
                        RuO_sum = rowSums(across(matches("RuO\\d?"))),
                        IrO_sum = rowSums(across(matches("IrO\\d?"))),
                        PdO_sum = rowSums(across(matches("PdO\\d?"))),
                        PtO_sum = rowSums(across(matches("PtO\\d?"))),
                        RhO_sum = rowSums(across(matches("RhO\\d?"))),) %>% 
  mutate(across(c(contains("sum")), ~ .x/rowSums(across(matches(str_all_primary_ion))),  
                .names = "{gsub('_sum$', '_ratio', .col)}")) %>% 
  select(sample , primary_ion, replicat, contains("ratio"), contains("sum")) %>% 
  pivot_longer(!c(sample, primary_ion, replicat) , names_to = c("oxides", "metric"),
                         values_to = "value", names_sep = "_") %>% 
  subset(metric == "ratio" 
         & sample %in% c("Ru", "Ru-Pt", "Ru-Pd", "Ru-Ir", "Ru-Rh", "Ru-Pt-Pd",
                         "Ru-Pt-Pd-Ir", "Ru-Pt-Pd-Ir-Rh", "HEA19")) 

# Calculate summary statistics
OR_df <- oxide_df %>%
  subset(value > 0) %>%
  group_by(sample, primary_ion, replicat) %>%
  reframe(OR = sum(value)) 

summary_df <- OR_df %>% 
  subset(primary_ion == "Bi5") %>% 
  group_by(sample) 

# Computation of the mean and sd of the OR for each sample
plot_df <- oxide_df %>%
  subset(primary_ion == "Bi5") %>% 
  subset(value > 0) %>%
  group_by(sample, oxides) %>%
  reframe(mean_ratio = mean(value)) %>%
  inner_join(summary_df %>% 
    reframe(OR_mean = mean(OR), OR_SE = sd(OR) / sqrt(n())),
  by = "sample") %>%
  mutate(sample = fct_reorder(sample, OR_mean, .desc = TRUE))
  
p <- plot_df %>% 
  ggplot(aes(x = sample, y = mean_ratio)) + 
  geom_bar(stat = "identity" , aes(fill = oxides)) + 
  scale_fill_manual(values=c("#3bc093", "#ffd479", "#e5a273", "#07b5bc", "#e55e57"),
                    labels = c(expression(IrO[x]), expression(PdO[x]), 
                               expression(PtO[x]),expression(RhO[x]),expression(RuO[x]))) +
  scale_x_discrete(labels=c("Ru-Ir" = "Ru-Ir", "Ru-Pd" = "Ru-Pd", "RuPtPd" = "Ru-Pt-Pd",
                               "RuPtPdIr" = "Ru-Pt-Pd-Ir", "Ru-Pt-Pd-Ir-Rh" = "HEA Mix",
                            "HEA19" = "HEA Seg")) +
  geom_errorbar(aes(ymin = OR_mean- 1.96*OR_SE, ymax = OR_mean + 1.96*OR_SE), 
                width = 0.2) +
  labs(x = NULL, y ="OR", fill = "")  +
  theme_bw() +
  theme(legend.position = "top",
        legend.text = element_text(size=35),
        axis.text = element_text(size=20),
        axis.title.y=element_text(size=30),
        axis.text.x = element_text(size=30, color = "black", angle = 90, hjust = 1, vjust = 0.5 ),
        axis.text.y = element_text(size=30, color = "black"),
        panel.border = element_rect(color = "grey"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)
    )
p

#Perform pairwise t-tests
pairwise.t.test(summary_df$OR, summary_df$sample, p.adjust.method = "bonferroni", pool.sd = T)



# =======================================================
#  - OXIDATION RATIO WITH ONE ELEMENT AT THE TIME
# =======================================================

anova_df <- scaled_automatic_negative_df |> mutate(
                        IrO_sum = rowSums(across(matches("^IrO"))),
                        PdO_sum = rowSums(across(matches("^PdO"))),
                        PtO_sum = rowSums(across(matches("^PtO"))),
                        RhO_sum = rowSums(across(matches("^RhO"))),)

str_all_primary_ion <- "^(Ru|Pt|Pd|Ir|Rh){1}$"

# Define colors for each oxide
oxide_colors <- c("PtO" = "#e5a273", "RhO" = "#07b5bc",
                  "IrO" = "#3bc093", "PdO" = "#ffd479")

# Define the color palette for the samples
sample_colors <- c("Ru" = "#e55e57" ,"Ru-Pt" = "#e5a273", "Ru-Pd" = "#ffd479",
                   "Ru-Ir" = "#3bc093", "Ru-Rh" = "#07b5bc", 
                   "Ru-Pt-Pd" = "#509bdd", "Ru-Pt-Pd-Ir" = "#8d85e5", 
                   "HEA Mix" = "#b373d0", "HEA Seg" = "#D452AA")

sample_names <- c("Ru" = "Ru", "Ru-Pt" = "Ru-Pt", "Ru-Pd" = "Ru-Pd",
                  "Ru-Ir" = "Ru-Ir","Ru-Rh" = "Ru-Rh",
                  "Ru-Pt-Pd" = "Ru-Pt-Pd", 
                  "Ru-Pt-Pd-Ir" = "Ru-Pt-Pd-Ir", 
                  "Ru-Pt-Pd-Ir-Rh" = "HEA Mix", 
                  "HEA19" = "HEA Seg")

# function to create the TEX Ratio
appender <- function(string){ 
  M <- substr(string, 1,2)
    TeX(paste0("$\\bf{OR_{" , M, "}", "}$"))
}

anova_df <- anova_df |> 
  mutate(
    IrO_ratio = IrO_sum/rowSums(across(matches(str_all_primary_ion))),
    PdO_ratio = PdO_sum/rowSums(across(matches(str_all_primary_ion))),
    PtO_ratio = PtO_sum/rowSums(across(matches(str_all_primary_ion))),
    RhO_ratio = RhO_sum/rowSums(across(matches(str_all_primary_ion)))) |>
  select(sample , primary_ion, contains("ratio"), contains("sum")) |> 
  pivot_longer(!c(sample, primary_ion), names_to = c("oxides", "fct"),
                         values_to = "value", names_sep = "_") |> 
  subset(value != 0 
         & fct == "ratio"
         & primary_ion == "Bi5" 
         & sample %in% c("Ru-Ir","Ru-Pt", "Ru-Pd", "Ru-Pt-Pd", "Ru-Rh",
                         "Ru-Pt-Pd-Ir", "Ru-Pt-Pd-Ir-Rh", "HEA19")) 

# ---------------- HISTOGRAM

p <- anova_df |> 
  group_by(sample, oxides) |> 
  summarise(mean_sum = mean(value), SE = sd(value)/sqrt(n()), .groups = "drop") %>% 
  mutate(sample = factor(sample, levels = names(sample_names),
                               labels = sample_names)) %>% 
  ggplot(aes(x = reorder_within(sample, -mean_sum, oxides, .fun = "sum") , 
             fill = sample,
             y = mean_sum,)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = mean_sum- 1.96*SE, ymax = mean_sum + 1.96*SE), 
                width = 0.2)+
  labs( y = NULL, fill = "") +
  scale_fill_manual(values = sample_colors) +
  scale_x_reordered() +
  facet_wrap(~oxides, scales = "free",
             labeller = as_labeller(appender, 
                            default = label_parsed)) + 
  theme(legend.position = "none",
        axis.text.x = element_text(size=25, color = "black", angle = 45,
                                   hjust = 1, vjust = 1 , face = "plain"),
        axis.text.y = element_text(size = 25),
        strip.text.x = element_text(size=25, face="bold"),
        panel.border = element_rect(color = "#666666")) 
  
  
# Function to get the strip background colors
strip_background_fill <- function(p) {
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-t', g$layout$name))
  fills <- oxide_colors
  
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k + 1
  }
  return(g)
}

axis_text_color <- function(plot, col = "fill") {
  c <- ggplot_build(plot)$data[[1]]
  plot +
    theme(axis.text.x = element_text(colour = c[[col]]))
}

axis_text_color(p)

g <- strip_background_fill(p)

# Draw the plot
grid::grid.draw(g)

#--------------- Pairwise Test 
# Function to perform pairwise t-test for each oxide
pairwise_test <- function(df) {
  test <- pairwise.t.test(df$value, df$sample, p.adjust.method = "none", pool.sd = F) %>% 
    tidy(test) %>% mutate(significant = ifelse(p.value < 0.05, "*", "ns"))
}

# Perform pairwise t-tests for each oxide group
results <- anova_df %>%
  group_by(oxides) %>%
  nest() %>%
  mutate(pairwise_results = map(data, pairwise_test))

# Extract and print the results
results <- results %>%
  select(oxides, pairwise_results) %>% 
  unnest(pairwise_results) %>%
  print(n = Inf)


# =========================================================
#  - Reviewer Comments: Plot of differents oxidation states
# =========================================================

scaled_automatic_negative_df %>% 
  subset(primary_ion == "Bi5" &
  sample %in% c("Ru-Pd", "Ru-Pt-Pd", "Ru-Pt-Pd-Ir", "Ru-Pt-Pd-Ir-Rh", "HEA19")) %>%
  select(sample, Pt, matches("^PdO\\d?$")) %>% # Change here the Oxide Spicies that you want to analyse
  group_by(sample) %>% 
  summarise(across(where(is.numeric), mean)) %>%
  pivot_longer(!c(sample), names_to = "Pd_oxide") %>% 
  ggplot(aes(x = sample, y = value, color = Pd_oxide, group = Pd_oxide)) +
  geom_point() + 
  geom_line() + 
  ylab("Scaled Intensities") + 
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=14, face = "bold")) 
