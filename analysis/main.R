# Main file to run all versions of PartSA and PSA
###############################################################################
# Author G. Cupples, E.Krebs. 10/10/2023

## Install/Load Required Packages
library('dampack')
library('pacman')
# package 'xlsx' requires Java installation
# MCMCpack requires JAGS installation
pacman::p_load('remotes','survival','prevalence', 'xlsx', 'readxl','MCMCpack', 'openxlsx',
               'tidyverse','truncnorm', 'ggplot2', 'ggrepel', 'scales', 'Hmisc')

## Load workspace and parameters
source('data-raw/01_parameter_load.R')

## Generate survival curve
source('R/02_curve_fit_functions.R')
source('data-raw/02_curve_fit.R')

## Load plotting functions for deterministic analysis 
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
# Cost-effectiveness plane
source('analysis/05b_CEplane.R')
# Cost-effectiveness frontier
source('R/05c_Frontier_Combinations.R')
source('analysis/05c_Frontier.R')
# Expected value of perfect information
source("analysis/05d_EVPI.R")


