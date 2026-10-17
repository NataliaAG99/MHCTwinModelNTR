#########################
#
# MHC heritability study
# 
# AE comparison
#
# 24/11/2025; R version 4.5.1
#
#########################
# Preparing the environment
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls(all=TRUE)) #Delete every previous object in the environment


#### Install and load packages ####
## First specify the packages of interest
packages = c("OpenMx", "psych", "umx", "haven", "readxl", "dplyr", "ggplot2",
             "lme4", "grid", "reshape2", "here", "moments", "foreign", "gee", "tidyr",
             "ISLR", "polycor", "tidyverse", "broom", "mosaic", "ggcorrplot")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


# Set the directory
#
#
# setwd("C:/Users/vgs055/OneDrive - Vrije Universiteit Amsterdam/2_Projecten/2025_Empathy heritability/2_Data analysis/2.2_Input")
#
# setwd("C:/Users/naa840/OneDrive - Vrije Universiteit Amsterdam/Nueva carpeta/MHC/03. Data/")
setwd("/data/natalia/MHC/")

fit5v <- readRDS("AE_5v_model.rds")
fit7v <- readRDS("AE_7v_model.rds")
fitcP <- readRDS("AE_cP_model.rds")

mxCompare(fit5v, fit7v)
mxCompare(fit7v, fitcP)
mxCompare(fit5v, fitcP)
  
summary(fit5v)

fitiP <- readRDS("AE_iP_model.rds")
mxCompare(fit7v, fitiP)

fitSAT <- readRDS("Sat_model.rds")
mxCompare(fitSAT, fit7v)
