# Function file containing all functions needed to run a construct curves from
# median survival data
################################################################################
# Author G.Cupples, E.Krebs. 10/10/2023

Generate_Survival_Curves <- function(df_data_SoC, df_data_treat, v_times, t_max_os = NULL, 
                                     t_max_pfs = NULL, switch_observed = 0){
  ## Generate Exponential Curves From Median Survival Data
  # INPUTS
  # df_data_SoC: tumour specific survival data for SoC
  # df_data_treat: tumour specific survival data for entrectinib
  # v_times: vector of time points to calculate curves on
  # t_max_os: maximum observation time for OS to switch from treatment to SoC
  # t_max_pfs: maximum observation time for PFS to switch from treatment to SoC
  # switch_observed: switch to either extrapolate to the end of time period or
  #                  switch to SoC after an observation period
  # OUTPUTS
  # results: list containing curve points and model fit in certain cases
  

  # outputs curve and exponential fit function
  #SoC
  l_curve_os_SoC  <- Exponential_Curve(df_data_SoC$OS , v_times)
  l_curve_pfs_SoC <- Exponential_Curve(df_data_SoC$PFS, v_times)
  
  #entrec
  if (switch_observed == 0){
    # full extrapolation
    l_curve_os_treat  <- Exponential_Curve(df_data_treat$OS , v_times, df_data_treat$at_risk_os)
    l_curve_pfs_treat <- Exponential_Curve(df_data_treat$PFS, v_times, df_data_treat$at_risk_pfs)
  } else if (switch_observed == 1){
    # restricted extrapolation - switch to SoC exponential rates after observed period
    l_curve_os_treat  <- Exponential_Curve_Restricted(df_data_treat$OS, df_data_treat$at_risk_os, 
                                                      t_max_os, v_times, l_curve_os_SoC[[1]])
    l_curve_pfs_treat <- Exponential_Curve_Restricted(df_data_treat$PFS, df_data_treat$at_risk_pfs, 
                                                      t_max_pfs, v_times, l_curve_pfs_SoC[[1]])
  }
  # outputs
  results <- list(os_SoC = l_curve_os_SoC, pfs_SoC = l_curve_pfs_SoC, 
                  os_treat = l_curve_os_treat, pfs_treat = l_curve_pfs_treat)
  
  return(results)
  
  
}



Exponential_Curve <- function(med_surv, times, at_risk_prop = 0.5){
  # construct curve from median survival data using an exponential extrapolation
  # INPUTS
  # med_surv: median survival value
  # times: vector of time points for each cycle in the model
  # OUTPUTS
  # l_curve_fit: list containing curve and fit function
  
  # set data points for extrapolation (median survival at 50%)
  x <- c(0, med_surv)
  y <- c(1, at_risk_prop)
  
  # create dataframe of points
  df <- data.frame(x, y)
  
  # fit exponential model
  exp_fit <- lm(log(y)~x , data = df)
  
  # calculate curves by tumour
  v_curve_fit  <- exp(predict(exp_fit , data.frame(x  = times)))
  
  # output the curve points and fit function if needed for extrapolation period
  l_curve_fit <- list(exp_fit, v_curve_fit)
  
  return(l_curve_fit)
  
}

Exponential_Curve_Restricted <- function(med_surv, at_risk_prop, t_max, times, 
                                         exp_comp){
  # construct exponential curve using median survival up until the end of the 
  # 'observed period' then use comparator data to extrapolate
  # INPUTS
  # med_surv: median survival value
  # at_risk_prop: the proportion of patients who remain at med_surv time
  # t_max: maximum follow-up time in the observed period
  # times: vector of time points for each cycle in the model
  # exp_comp: model fit for extrapolated calculation
  # OUTPUTS
  # l_curve_fit: list containing curve and fit function
  
  # set data points for extrapolation (median survival at 50%)
  x <- c(0, med_surv)
  y <- c(1, at_risk_prop)
  
  # create dataframe of points
  df  <- data.frame(x, y)
  
  # fit exponential model
  exp_fit  <- lm(log(y) ~ x , data = df)
  
  # Construct exponential curves
  if (is.na(t_max)){
    # if there is no t_max specified then set to full time period
    l_curve_fit <- list(exp_fit, exp(predict(exp_fit , data.frame(x = times))))
    
  } else {
    # if there is a t_max, split times to observed and extrapolated vectors
    
    t_obs <- seq(from = 0, to = t_max, by = cycle_length)              # observation time
    t_ext <-  seq(max(t_obs) + cycle_length, max(times), cycle_length) # extrapolation time
    
    
    # observation curve
    v_curve_fit <- exp(predict(exp_fit, data.frame(x = t_obs)))
    
    v_exp_ext  <- min(v_curve_fit) * exp(exp_comp$coefficients[[1]] +  
                                           exp_comp$coefficients[[2]] *(t_ext - min(t_ext) + cycle_length))
    
    # combine curves
    l_curve_fit  <- list(exp_fit, c(v_curve_fit, v_exp_ext))
    
  }
  
  
  return(l_curve_fit)
  
  
}