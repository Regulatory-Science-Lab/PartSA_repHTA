# Main file to run all versions of PartSA and PSA
###############################################################################
# Author G. Cupples, E.Krebs. 10/10/2023

## Install/Load Required Packages
if (!require('remotes'))    install.packages('remotes'  , repos = "http://cran.us.r-project.org"); library(remotes)
if (!require('survival'))   install.packages('survival' , repos = "http://cran.us.r-project.org"); library(survival)
if (!require('readxl'))     install.packages('readxl'   , repos = "http://cran.us.r-project.org"); library(readxl)
if (!require('dplyr'))      install.packages('dplyr'    , repos = "http://cran.us.r-project.org"); library(dplyr)
if (!require('truncnorm'))  install.packages('truncnorm', repos = "http://cran.us.r-project.org"); library(truncnorm)
if (!require('MCMCpack'))   install.packages('MCMCpack'  , repos = "http://cran.us.r-project.org"); library(MCMCpack)
if (!require('prevalence')) install.packages('prevalence', repos = "http://cran.us.r-project.org"); library(prevalence)
if (!require('dampack'))    install_github('DARTH-git/dampack'); library(dampack)

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
source('analysis/05_summary_outcomes.R')
