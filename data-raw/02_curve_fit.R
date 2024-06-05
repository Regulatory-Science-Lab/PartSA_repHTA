# Script to construct parametric curves from median survival data
###############################################################################
# Author G. Cupples, E.Krebs. 10/10/2023

## t_max for treatment
df_data_treat <- l_parameters$`3.2_MedianOSPFS_treat`[-c(16)] # without notes column
df_data_treat <- df_data_treat[order(df_data_treat$tumour),]
df_data_treat$t_max_pfs <- as.numeric(df_data_treat$t_max_pfs) /12
df_data_treat$t_max_os  <- as.numeric(df_data_treat$t_max_os ) /12


## Survival data

# convert to years
df_data_treat$PFS <- df_data_treat$PFS / 12
df_data_treat$OS  <- df_data_treat$OS  / 12
df_data_treat     <- df_data_treat[order(df_data_treat$tumour),]

# convert standard error to standard deviation for truncated normal, in cases where
# the efficacy is given by the trial and not the comparator
# pfs
df_data_treat[(df_data_treat$dist_pfs == "truncated normal" & !is.na(df_data_treat$t_max_pfs)),'par2_pfs'] <- 
  sqrt(df_weights[(df_data_treat$dist_pfs == "truncated normal" & !is.na(df_data_treat$t_max_pfs)),]$Trial) *
  df_data_treat[(df_data_treat$dist_pfs == "truncated normal" & !is.na(df_data_treat$t_max_pfs)),'par2_pfs']
# os
df_data_treat[(df_data_treat$dist_os == "truncated normal" & !is.na(df_data_treat$t_max_os)),'par2_os'] <- 
  sqrt(df_weights[(df_data_treat$dist_os == "truncated normal" & !is.na(df_data_treat$t_max_os)),]$Trial) *
  df_data_treat[(df_data_treat$dist_os == "truncated normal" & !is.na(df_data_treat$t_max_os)),'par2_os']

## Standard of Care
df_data_SoC     <- l_parameters$`3.1_MedianOSPFS_SoC`[-c(12)] # without notes column
df_data_SoC    <- df_data_SoC[order(df_data_SoC$tumour),]
df_data_SoC$OS  <- df_data_SoC$OS /12  # scale to years
df_data_SoC$PFS <- df_data_SoC$PFS/12  # scale to years

## Construct survival curves
l_curves    <- list()
l_model_fit <- list()
for (tumour in v_tumour){
  t_max_pfs <- df_data_treat[df_data_treat$tumour == tumour,]$t_max_pfs
  t_max_os  <- df_data_treat[df_data_treat$tumour == tumour,]$t_max_os
  
  # median data
  df_surv_treat <- df_data_treat[df_data_treat$tumour == tumour,]
  df_surv_SoC   <- df_data_SoC[df_data_SoC$tumour == tumour,]
  
  
  temp <- Generate_Survival_Curves(df_surv_SoC, df_surv_treat, v_times, 
                                   t_max_os, t_max_pfs, switch_observed)
  l_curves[[tumour]] <- list(os_SoC   = temp$os_SoC[[2]]  , pfs_SoC   = temp$pfs_SoC[[2]],
                             os_treat = temp$os_treat[[2]], pfs_treat = temp$pfs_treat[[2]])
  # store fit functions for PSA
  l_model_fit[[tumour]] <- list(os_SoC   = temp$os_SoC[[1]]  , pfs_SoC   = temp$pfs_SoC[[1]],
                                os_treat = temp$os_treat[[1]], pfs_treat = temp$pfs_treat[[1]])
}



