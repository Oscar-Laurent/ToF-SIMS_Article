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
