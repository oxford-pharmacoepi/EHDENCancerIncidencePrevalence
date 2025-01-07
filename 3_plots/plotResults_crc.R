library(dplyr)
library(tidyr)
library(stringr)
library(here)
library(ggplot2)
library(scales)
library("cowplot")
library("gridExtra")
library("ggpubr")
library(ggh4x)
library(patchwork)

pathResults <- "C:/Users/dnewby/Desktop/"

datapath <- "C:/Users/dnewby/Documents/GitHub/EHDENCancerIncidencePrevalence/CancerIncidencePrevalanceShiny/shiny/data"

# read in files
incidence_estimates <- readRDS(paste0(datapath ,"/incidence_estimates.rds")) %>% 
  mutate(database_name = ifelse(database_name == "CPRDGoldUpdate2", "CPRD GOLD", database_name))

