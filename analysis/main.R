# Main file to run all versions of PartSA and PSA
###############################################################################
# Author G. Cupples, E.Krebs. 10/10/2023

## Install/Load Required Packages
install.packages(c('pacman','remotes','survival','readxl','MCMCpack','prevalence',
'xlsx', 'openxlsx','tidyverse','dampack','truncnorm'))

library('pacman')
pacman::p_load('remotes','survival','readxl','MCMCpack','prevalence',
               'xlsx', 'openxlsx','tidyverse','dampack','truncnorm')


# library(remotes)
# library(survival)
# library(readxl)
# library(dplyr)
# library(truncnorm)
# library(MCMCpack)
# library(prevalence)
# library(xlsx)
# library(openxlsx)
# library(tidyverse)
# library(ggplot2)
# library(dampack)


## Load workspace and parameters
source('data-raw/01_parameter_load.R')

## Generate survival curves
source('R/02_curve_fit_functions.R')
source('data-raw/02_curve_fit.R')

## Load plotting functions
source('R/03_plotting_functions.R')

## Run deterministic analysis
source('R/04a_PartSA_functions.R')
source('analysis/04a_deterministic_analysis.R')

## Run probabilistic analysis
source('R/04a_PartSA_functions.R')
source('R/04b_PSA_functions.R')
source('analysis/04b_probabilistic_analysis.R')

## Manuscript summary outcomes
source('analysis/05a_summary_outcomes.R')
source('analysis/05b_CEplane.R')
