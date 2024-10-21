# Function file containing all functions needed to run a probabilistic
# tumour-agnostic partitioned survival analysis (PartSA) with a variety of 
# input survival data 
################################################################################
# Author G.Cupples, E.Krebs. 10/10/2023

run_PSA <- function(df_data_treat, df_data_SoC, n_psa, seed, 
                    switch_observed){
  ## function to run PSA of PartSA framework
  # INPUTS
  # df_data_treat: treatment survival data
  # df_data_SoC: SoC survival data
  # n_psa: number of PSA runs to do
  # seed: seed for randomization
  # switch_observed: if the efficacy of the treatment is restricted beyond the 
  #                  observed period
  # OUTPUTS
  # l_PSA_results: List containing stratified and pooled outcomes for each PSA run
  
  ## generate parameter space
  l_parameter_space <- Generate_Parameter_Space(n_psa, seed)
  
  
  ## run PartSA for each parameter set
  l_PSA_stratified <- list()
  df_PSA_weighted  <- data.frame()
  for (i in 1:n_psa){
    temp <- Parameter_Setup(df_data_treat, l_parameter_space, i)
    l_c_treat_psa  <- temp$l_c_treat
    l_c_SoC_psa    <- temp$l_c_SoC
    l_u_treat_psa  <- temp$l_u_treat
    l_u_SoC_psa    <- temp$l_u_SoC
    l_curves_psa   <- temp$l_curves_psa
    df_weights_psa <- temp$weights
    
    ## Run Analysis
    # summary outputs only and no plots 
    l_outcomes <- run_PartSA(df_data_treat, l_curves_psa, l_c_treat_psa, l_c_SoC_psa, 
                           l_u_treat_psa, l_u_SoC_psa, v_times, df_weights_psa, 
                           switch_observed, switch_summary = 1, 
                           switch_plot = 0)
    
    
    ## store individual and weighted outcomes
    
    # individual
    for (tumour in v_tumour){
      l_PSA_stratified[[tumour]] <- rbind(l_PSA_stratified[[tumour]], 
                                          l_outcomes$stratified_outcomes[tumour,])
    }

    # weighted
    df_PSA_weighted <- rbind(df_PSA_weighted, l_outcomes$weighted_outcomes)
  }
  
  # create PSA object per tumour indication
  l_PSA_obj <- list()
  for (tumour in v_tumour){
    l_PSA_obj[[tumour]] <- make_psa_obj(cost          = as.data.frame(l_PSA_stratified[[tumour]][,c(1,4)]),
                                        effectiveness = as.data.frame(l_PSA_stratified[[tumour]][,c(3,6)]),
                                        parameters    = l_outcomes$parameter_space[[tumour]],
                                        strategies    = v_treat_names)
  }
  
  
  l_PSA_results <- list(stratified_outcomes = l_PSA_stratified, stratified_obj = l_PSA_obj, weighted_outcomes = df_PSA_weighted)
  return(l_PSA_results)
  
}

Parameter_Setup <- function(df_data_treat, l_parameter_space, ii){
  ## function to convert PSA parameters into PartSA input format
  # INPUTS
    # df_data_treat: survival data for treatment
    # l_parameter_space: set of PSA parameters
    # ii: number corresponding to the PSA run
  # OUTPUTS
    # Outputs: list containing costs, utilities, and survival data for treatment and SoC
  
  l_c_treat       <- list()
  l_c_SoC         <- list()
  l_u_treat       <- list()
  l_u_SoC         <- list()
  l_curves_psa    <- list()
  df_weights_psa  <- data.frame(tumour = df_data_treat$tumour, 
                                Clinical = rep(NA, length(df_data_treat$tumour)))
  for (tumour in df_data_treat$tumour){
    df_param_space <- l_parameter_space[[tumour]][ii,]
    df_weights_psa[df_weights_psa$tumour == tumour,]$Clinical = df_param_space$weights
    
    
    ## Costs
    
    # SoC
    l_c_SoC[[tumour]] <- list(
      c_pfs  = df_param_space$c_pfs_SoC,
      c_P    = df_param_space$c_P_SoC,
      c_D    = df_param_space$c_D_SoC,
      c_AE   = df_param_space$c_AE_SoC,
      c_test = df_param_space$c_test_SoC
    )
    
    # treatment
    # progression-free
    v_c_pfs_treat <- rep(df_param_space$c_pfs_treat, length(v_times))
    
    # treatment costs up until time to off treatment
    v_c_pack_treat_P <- rep(0, length(v_times))
    v_c_pack_treat_P[1:which.min(abs(v_times - ttot_fixed))] <- c_pack_treat_fixed
    # treatment administration costs up until time to off treatment
    v_c_admin_treat_P <- rep(0, length(v_times))
    v_c_admin_treat_P[1:which.min(abs(v_times - ttot_fixed))] <- c_admin_treat_fixed
    # progressed
    v_c_P_treat <- df_param_space$c_P_treat + v_c_pack_treat_P + v_c_admin_treat_P
    
    # death
    v_c_D_treat <- rep(df_param_space$c_D_treat, length(v_times))
    
    l_c_treat[[tumour]] <- list(
      c_pfs  = v_c_pfs_treat,
      c_P    = v_c_P_treat,
      c_D    = v_c_D_treat,
      c_AE   = df_param_space$c_AE_treat,
      c_test = df_param_space$c_test_treat
    )
    
    ## Uilities
    l_u_treat[[tumour]] <- list(
      u_pfs = df_param_space$u_pfs_treat,
      u_P   = df_param_space$u_P_treat,
      u_D   = df_param_space$u_D
    )
    l_u_SoC[[tumour]] <- list(
      u_pfs = df_param_space$u_pfs_SoC,
      u_P   = df_param_space$u_P_SoC,
      u_D   = df_param_space$u_D
    )
    
    ## Surv
    SoC_vec   <- data.frame(tumour      = tumour,
                            PFS         = df_param_space$PFS_SoC,
                            OS          = df_param_space$OS_SoC)
    treat_vec <- data.frame(tumour      = tumour,
                            at_risk_pfs = df_param_space$at_risk_pfs,
                            PFS         = df_param_space$PFS_treat,
                            t_max_pfs   = df_data_treat[df_data_treat$tumour == tumour,]$t_max_pfs,
                            at_risk_os  = df_param_space$at_risk_os,
                            OS          = df_param_space$OS_treat,
                            t_max_os    = df_data_treat[df_data_treat$tumour == tumour,]$t_max_os)
    
    #SoC
    l_curve_os_SoC  <- Exponential_Curve(SoC_vec$OS , v_times)
    l_curve_pfs_SoC <- Exponential_Curve(SoC_vec$PFS, v_times)
    
    # treatment
    if (switch_observed == 0){
      # full extrapolation
      l_curve_os_treat   <- Exponential_Curve(treat_vec$OS , v_times)
      l_curve_pfs_treat  <- Exponential_Curve(treat_vec$PFS, v_times)
    } else if (switch_observed == 1){
      # restricted extrapolation - switch to SoC exponential rates after observed period
      l_curve_os_treat  <- Exponential_Curve_Restricted(treat_vec$OS, treat_vec$at_risk_os, 
                                                        treat_vec$t_max_os, v_times, l_curve_os_SoC[[1]])
      l_curve_pfs_treat <- Exponential_Curve_Restricted(treat_vec$PFS, treat_vec$at_risk_pfs, 
                                                        treat_vec$t_max_pfs, v_times, l_curve_pfs_SoC[[1]])
    }
    
    l_curves_psa[[tumour]] <- list(os_SoC   = l_curve_os_SoC[[2]]  , pfs_SoC   = l_curve_pfs_SoC[[2]], 
                                   os_treat = l_curve_os_treat[[2]], pfs_treat = l_curve_pfs_treat[[2]])
    
  }
  
  
  Outputs <- list(l_c_SoC = l_c_SoC, l_c_treat = l_c_treat,
                  l_u_SoC = l_u_SoC, l_u_treat = l_u_treat,
                  l_curves_psa = l_curves_psa, weights = df_weights_psa)
  
  return(Outputs)
  
}

Generate_Parameter_Space <- function(n_psa, seed = 50){
  # function to generate parameter space for PSA
  # INPUTS
  # n_psa: number of PSA runs
  # seed: seed for randomization
  # OUTPUTS
  # l_param_space: list per tumour containing all PSA parameters
  
  
  set.seed(seed)
  
  ## Weights
  # vary overall prevalence with Dirichlet distribution
  df_weights_psa <- MCMCpack::rdirichlet(n_psa - 1, df_weights$Clinical)
  df_weights_psa <- rbind((df_weights$Clinical / sum(df_weights$Clinical)), df_weights_psa)
  df_weights_psa <- as.data.frame(df_weights_psa)
  colnames(df_weights_psa) <- df_weights$tumour
  
  ## TTOT
  ttot_psa <- param_dist(l_parameters$`3.6_TTOT`, n_psa)
  
  ## Utilities - tumour agnostic
  for (var in l_parameters$`2.1_Utilities`$parameter){
    name <- paste("v_", var, "_psa", sep = "")
    dist <- param_dist(df_utilities[df_utilities$parameter == var,], n_psa)
    assign(name, dist)
  }
  
  ## Costs - tumour agnostic
  for (var in l_parameters$`1.1_TumourAgnosticCosts`$parameter){
    name <- paste("v_", var, "_psa", sep = "")
    dist <- param_dist(df_costs_agnostic[df_costs_agnostic$parameter == var,], n_psa)
    assign(name, dist)
  }

  ## Tumour specific parameters
  l_param_space <- list()
  for (tumour in v_tumour){
    
    ## Costs
    v_c_admin_SoC_psa <- param_dist(df_c_admin_SoC[df_c_admin_SoC$tumour == tumour,], n_psa)
    v_c_pack_SoC_psa  <- param_dist(df_c_pack_SoC[df_c_pack_SoC$tumour == tumour,], n_psa)
    v_c_noncancer_psa <- param_dist(df_c_noncancer[df_c_noncancer$tumour == tumour,], n_psa)
    
    ## Utility Check
    if (is.na(df_data_treat[df_data_treat$tumour == tumour,]$t_max_pfs)){
      # if SoC efficacy is being assumed then we have no incremental benefits
      temp_pfs <- v_u_pfs_SoC_psa
      temp_P   <- v_u_P_SoC_psa
    } else{
      temp_pfs <- v_u_pfs_treat_psa
      temp_P   <- v_u_P_treat_psa
    }
    
    
    ## Dataframe
    df_param_space <- data.frame(
      weights      = df_weights_psa[[tumour]], 
      c_pack_treat = v_c_pack_treat_psa,  
      c_test_SoC   = v_c_test_SoC_psa,     
      c_D_treat    = v_c_D_psa,            
      c_pfs_treat  = v_c_admin_treat_psa + v_c_care_pfs_treat_psa + 
                     v_c_noncancer_psa + v_c_pack_treat_psa,
      c_P_treat    = v_c_care_P_psa + v_c_noncancer_psa,
      c_D_SoC      = v_c_D_psa,
      c_pfs_SoC    = v_c_admin_SoC_psa + v_c_care_pfs_SoC_psa + 
                     v_c_pack_SoC_psa + v_c_noncancer_psa,
      c_P_SoC      = v_c_care_P_psa + v_c_noncancer_psa,
      c_AE_SoC     = v_c_AE_SoC_psa,
      c_AE_treat   = v_c_AE_treat_psa,
      
      u_D          = v_u_D_psa,            
      u_pfs_SoC    = v_u_pfs_SoC_psa,
      u_P_SoC      = v_u_P_SoC_psa,
      u_pfs_treat  = temp_pfs,
      u_P_treat    = temp_P
      
    )
    
    l_param_space[[tumour]] <- df_param_space
    
  }
  
  # testing costs
  if (switch_test == 1){
    df_c_test_treat_psa <- data.frame(matrix(NA, nrow = n_psa, ncol = length(df_data_treat$tumour)))
    colnames(df_c_test_treat_psa) <- v_tumour
    for (tumour in v_tumour[v_tumour != "Other"]){
      df_c_test_treat_psa[[tumour]] <- param_dist(df_test_temp[df_test_temp$tumour == tumour,], n_psa)
    }
    # other costs 
    other_min <- 0.8 * df_c_test_treat$Other
    other_max <- 1.2 * df_c_test_treat$Other
    df_c_test_treat_psa$Other <- c(df_c_test_treat$Other,
                                   runif(n_psa - 1, other_min, other_max))
  } else {
    # no testing costs
    df_c_test_treat_psa <- data.frame(matrix(0, nrow = n_psa, ncol = length(df_data_treat$tumour)))
    colnames(df_c_test_treat_psa) <- df_data_treat$tumour
  }
  
  for (tumour in v_tumour){
    l_param_space[[tumour]]$c_test_treat <- df_c_test_treat_psa[[tumour]]
  }
  
  ## Survival
  set.seed(seed)
  for (tumour in v_tumour){
    # survival parameters based off median values
    Surv <- Survival_PSA_Parameters(df_data_SoC[df_data_SoC$tumour == tumour,], 
                                    df_data_treat[df_data_treat$tumour == tumour,], 
                                    n_psa)
    
    l_param_space[[tumour]]$OS_SoC      <- Surv$v_OS_SoC
    l_param_space[[tumour]]$PFS_SoC     <- Surv$v_PFS_SoC
    
    l_param_space[[tumour]]$OS_treat    <- Surv$v_OS_treat
    l_param_space[[tumour]]$at_risk_os  <- Surv$v_at_risk_os
    l_param_space[[tumour]]$PFS_treat   <- Surv$v_PFS_treat
    l_param_space[[tumour]]$at_risk_pfs <- Surv$v_at_risk_pfs
    
  }
  
  return(l_param_space)
}

param_dist <- function(df_params, n_psa, upper = NULL, lower = NULL){
  ## construct distribution based on distribution type and parameters
  # INPUTS
  # df_params: data.frame with details of uncertainty for parameter of interest
  # n_psa: number of parameters to produce
  # upper: upper limit for truncated normal, where needed
  # lower: lower limit for truncated normal, where needed
  # OUTPUTS
  # v_distribution: distribution for the chosen parameter
  
  v_distribution <- rep(NA, n_psa)
  if (df_params$dist == "normal"){
    v_distribution[2:n_psa] <- rnorm(n_psa - 1, as.numeric(df_params$par1),
                                   as.numeric(df_params$par2))
  }
  else if (df_params$dist == "truncated normal"){
    v_distribution[2:n_psa] <- rtruncnorm(n_psa - 1, a = 0, mean = as.numeric(df_params$par1),
                                        sd = as.numeric(df_params$par2))
  } else if (df_params$dist == "uniform"){
    v_distribution[2:n_psa] <- runif(n_psa - 1, as.numeric(df_params$par1),
                                   as.numeric(df_params$par2))
  } else if(df_params$dist == "fixed"){
    v_distribution[2:n_psa] <- rep(as.numeric(df_params$value), n_psa - 1)
  } else if (df_params$dist == "beta"){
    v_distribution[2:n_psa] <- rbeta(n_psa - 1, as.numeric(df_params$par1),
                                   as.numeric(df_params$par2))
  } else if (df_params$dist == "beta estimate") {
    # estimate beta distribution parameters from 95% CI
    temp   <- prevalence::betaExpert(best  = df_params$value, lower = as.numeric(df_params$par1), 
                                     upper = as.numeric(df_params$par2), p = 0.95, method = "mean")
    shape1 <- temp$alpha
    shape2 <- temp$beta
    # construct beta distribution
    v_distribution[2:n_psa] <- rbeta(n_psa - 1, shape1, shape2)
  }
  v_distribution[1] <- as.numeric(df_params$value) # set first point to base case
  
  return(v_distribution)
}


Survival_PSA_Parameters <- function(df_data_SoC, df_data_treat, n_psa){
  # function to calculate PSA parameters for median/critical point survival data
  # INPUTS
    # df_data_SoC: survival data for SoC, contains distribution information
    # df_data_treat: survival data for treatment, contains distribution information
    # n_psa: number of PSA runs
  # OUTPUTS
    # df_results: data frame containing OS, PFS, and at risk proportion for SoC and treatment
  
  ## SoC
  os_params_SoC  <- c(df_data_SoC$par1_os,  df_data_SoC$par2_os )
  pfs_params_SoC <- c(df_data_SoC$par1_pfs, df_data_SoC$par2_pfs)
  
  v_OS_SoC <- surv_dist(os_params_SoC, df_data_SoC$dist_os, n_psa - 1)
  v_OS_SoC <- c(df_data_SoC$par1_os, v_OS_SoC) # set deterministic as first parameter set
  
  # PFS
  v_PFS_SoC <- rep(0, n_psa)
  v_PFS_SoC[1] <- df_data_SoC$par1_pfs # set deterministic as first parameter set
  for (i in 2:n_psa){
    # use OS as max for truncated normal
    v_PFS_SoC[i] <- surv_dist(pfs_params_SoC, df_data_SoC$dist_pfs, 1, max = v_OS_SoC[i])
  }

  ## Treatment
  if (is.na(df_data_treat$t_max_os)){
    # if t_max is NA, then SoC efficacy is assumed throughout
    v_OS_treat      <- v_OS_SoC
    v_at_risk_os  <- df_data_treat$at_risk_os
    v_PFS_treat     <- v_PFS_SoC
    v_at_risk_pfs <- df_data_treat$at_risk_pfs
  } else{
    os_params_treat  <- c(df_data_treat$par1_os,  df_data_treat$par2_os )
    pfs_params_treat <- c(df_data_treat$par1_pfs, df_data_treat$par2_pfs)
    if (df_data_treat$dist_os == "truncated normal"){
      # OS
      v_OS_treat   <- surv_dist(os_params_treat, df_data_treat$dist_os, n_psa - 1)
      v_OS_treat   <- c(df_data_treat$par1_os, v_OS_treat) # set deterministic as first parameter set
      v_at_risk_os <- rep(df_data_treat$at_risk_os, n_psa) # fixed  
      
      # PFS
      if (df_data_treat$dist_pfs == "truncated normal"){
        v_PFS_treat    <- rep(0, n_psa)
        v_PFS_treat[1] <- df_data_treat$par1_pfs # set deterministic as first parameter set
        for (i in 2:n_psa){
          v_PFS_treat[i] <- surv_dist(pfs_params_treat, df_data_treat$dist_pfs, 1, max = v_OS_treat[i])
        }
        v_at_risk_pfs <- rep(df_data_treat$at_risk_pfs, n_psa) # fixed
      } else if (df_data_treat$dist_pfs == "beta"){
        v_PFS_treat   <- rep(df_data_treat$PFS * 12, n_psa) # fixed
        v_at_risk_pfs <- surv_dist(pfs_params_treat, df_data_treat$dist_pfs, n_psa - 1)
        v_at_risk_pfs <- c(df_data_treat$at_risk_pfs, v_at_risk_pfs) # set deterministic as first parameter set
        
      }
      
    } else if (df_data_treat$dist_os == "beta"){
      # OS
      v_OS_treat   <- df_data_treat$OS * 12 # fixed
      v_at_risk_os <- surv_dist(os_params_treat, df_data_treat$dist_os, n_psa - 1)
      v_at_risk_os <- c(df_data_treat$at_risk_os, v_at_risk_os) # set deterministic as first parameter set
      
      # PFS
      if (df_data_treat$dist_pfs == "truncated normal"){
        # since OS is at a 'critical' point instead of the median, it cannot be used
        # as an upper limit for pfs
        v_PFS_treat   <- surv_dist(pfs_params_treat, df_data_treat$dist_pfs, n_psa - 1)
        v_PFS_treat   <- c(df_data_treat$par1_pfs, v_PFS_treat) # set deterministic as first parameter set
        v_at_risk_pfs <- df_data_treat$at_risk_pfs # fixed
        
      } else if (df_data_treat$dist_pfs == "beta"){
        v_PFS_treat   <- df_data_treat$PFS * 12 # fixed
        v_at_risk_pfs <- surv_dist(pfs_params_treat, df_data_treat$dist_pfs, n_psa - 1)
        v_at_risk_pfs <- c(df_data_treat$at_risk_pfs, v_at_risk_pfs) # set deterministic as first parameter set
      }
    }
    
  }
  # convert to years
  v_OS_SoC    <- v_OS_SoC   /12
  v_PFS_SoC   <- v_PFS_SoC  /12
  v_OS_treat  <- v_OS_treat /12
  v_PFS_treat <- v_PFS_treat/12
  
  # dataframe containing 6 columns with OS/PFS and at risk for treatment
  df_results <- data.frame(v_OS_SoC, v_PFS_SoC, v_OS_treat, v_at_risk_os, v_PFS_treat, v_at_risk_pfs)
  
  return(df_results)
}

surv_dist <- function(v_params, dist, n_psa, max = ""){
  ## construct distribution based on distribution type and parameters
  # INPUTS
    # v_params: vector containing the parameters needed to fit distribution
    # dist: distribution type (e.g. normal, beta)
    # n_psa: number of parameters to produce
    # max: upper limit for truncated normal for PFS, where needed
  # OUTPUTS
    # v_dist: distribution for the chosen parameter
  
  if (dist == "normal") {
    v_dist <- rnorm(n_psa, mean = v_params[1], sd = v_params[2])
  } else if (dist == "truncated normal") {
    if (max == ""){
      v_dist <- rtruncnorm(n_psa, a = 0.5, mean = v_params[1], sd = v_params[2])
    } else if (max != "") {
      # if there is an upper bound for the truncated normal
      v_dist <- rtruncnorm(n_psa, a = 0.5, b = max, mean = v_params[1], sd = v_params[2])
    }
  } else if (dist == "beta") {
    v_dist <- rbeta(n_psa, shape1 = v_params[1], shape2 = v_params[2])
  } else if (dist == "gamma") {
    v_dist <- rgamma(n_psa, shape = v_params[1], rate = v_params[2])
  } else if (dist == "uniform") {
    v_dist <- runif(n_psa, min = v_params[1], max = v_params[2])
  } else if (dist == "fixed") {
    v_dist <- v_params[1]
  } else if (dist == "exponential"){
    v_dist <- rexp(n_psa, rate = 1/v_params[1])
  }
  return(v_dist)
}

PSA_frontier <- function(v_tumour_groups, outcomes_PSA, save_name){
  # construct and save frontier for PSA outputs
  # INPUTS
  # v_tumour_groups: vector of tumour indications included in the analysis
  # outcomes_PSA: list of outcomes for the PSA
  # save_name: filename for csv output
  # OUTPUTS
  # df_summary_totals: data frame of incremental outputs for each tumour combination
  
  # all possible combinations of tumour type
  combs <- lapply(1:length(v_tumour_groups), function(y){combn(v_tumour_groups, y)})
  
  save_ind <- 0
  l_totals <- list()
  for (ind in seq(1,length(combs))){
    for (comb_ind in seq(1, dim(combs[[ind]])[2])){
      # select specific tumour from combination
      comb_sub   <- combs[[ind]][,comb_ind]
      save_ind <- save_ind + 1
      # sum the costs/qalys/icer for each run for all tumour combinations
      v_cost <- 0
      v_qaly <- 0
      v_icer <- 0
      for (tumour in comb_sub){
        v_cost <- v_cost + mean(outcomes_PSA$stratified_outcomes[[tumour]]$Inc_Cost)
        v_qaly <- v_qaly + mean(outcomes_PSA$stratified_outcomes[[tumour]]$Inc_QALYs)
      }
      v_icer <- v_cost / v_qaly

      l_totals[[save_ind]] <- data.frame(Combination = paste(comb_sub, collapse = ", "),
                                         Costs = v_cost, 
                                         QALYs = v_qaly,
                                         ICERs = v_icer)
    }
    
  }
  
  df_summary_totals <- do.call("rbind", l_totals)
  
  write.csv(df_summary_totals, save_name)
  
  return(df_summary_totals)
}

