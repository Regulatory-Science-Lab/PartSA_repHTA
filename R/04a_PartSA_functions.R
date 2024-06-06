# Function file containing all functions needed to run a deterministic
# tumour-agnostic partitioned survival analysis (PartSA) with a variety of 
# input survival data
################################################################################
# Author G.Cupples, E.Krebs. 10/10/2023

run_PartSA <- function(df_data_treat, l_curves, l_c_treat, l_c_SoC, l_u_treat, 
                       l_u_SoC, v_times, df_weights, switch_observed = 0,
                       switch_summary = 1, switch_plot = 0){
  # Main function to run PartSA framework
  # INPUTS
  # df_data_treat: dataframe with information about the treatment arm
  # l_curves: list of parametric curve datapoints per tumour
  # l_c_: list of costs per tumour
  # l_u_: list of utilities per tumour
  # v_times: vector of time points at each cycle in the model
  # df_weights: dataframe of weights for each tumour
  # switch_observed: determines if efficacy is restricted for treatment group
  # switch_summary: determines vector or summary outputs
  # switch_plot: determines whether to create and save plots
  # OUTPUTS
  # l_results: a list containing treatment specific and incremental outcomes per tumour
  
  ## Calculate trace
  l_trace_SoC   <- lapply(l_curves, function(x) Part_Surv(x$pfs_SoC  , x$os_SoC))
  l_trace_treat <- lapply(l_curves, function(x) Part_Surv(x$pfs_treat, x$os_treat))
  
  ## Calculate costs and QALYs
  l_out_SoC <- mapply(Calculate_Cost_QALY, l_trace_SoC, l_c_SoC, l_u_SoC,
                      MoreArgs = list(v_discount_cost = v_discount_c, 
                                      v_discount_qaly = v_discount_q), SIMPLIFY = FALSE)
  l_out_treat <- mapply(Calculate_Cost_QALY, l_trace_treat, l_c_treat, l_u_treat,
                        MoreArgs = list(v_discount_cost = v_discount_c, 
                                        v_discount_qaly = v_discount_q,
                                        switch_vec_cost = 1), SIMPLIFY = FALSE)
  
  ## Calculate incremental
  # create list of t_max
  l_t_max <- as.list(df_data_treat[match(v_tumour, df_data_treat$tumour),]$t_max_pfs)
  names(l_t_max) <- v_tumour
  
  l_incremental <- mapply(Calculate_Incremental, l_out_treat, l_out_SoC, l_t_max,
                          MoreArgs = list(v_times = v_times), SIMPLIFY = FALSE)
  
  ## Plotting
  if (switch_plot == 1){
    for (tumour in df_data_treat$tumour){
      Create_Plots(l_curves[[tumour]], l_trace_treat[[tumour]], l_trace_SoC[[tumour]], v_times,
                   tumour, t_max = df_data_treat[df_data_treat$tumour == tumour,]$t_max_pfs)
    }
  }

  
  
  ## Results
  
  # summary outcomes per tumour indication
  df_sum_treat <- as.data.frame(do.call('rbind', lapply(l_out_treat, colSums)))
  colnames(df_sum_treat) <- c("Cost_treat", "LY_treat","QALY_treat")
  df_sum_SoC   <- as.data.frame(do.call('rbind', lapply(l_out_SoC, colSums)))
  colnames(df_sum_SoC)   <- c("Cost_SoC", "LY_SoC","QALY_SoC")
  df_sum_inc   <- as.data.frame(do.call('rbind', lapply(l_incremental, function(x) colSums(x[,1:3]))))
  # collate NMB, icer and extrapolated effect data
  df_inc_additional <- as.data.frame(do.call('rbind', lapply(l_incremental, function(x) x[1,4:8])))
  
  # weighted outcomes
  df_weighted_outcomes <- Weighted_Outcomes(df_sum_SoC, df_sum_treat, df_inc_additional, df_weights)
  
  
  if (switch_summary ==  1){
    # summary only
    df_stratified_outcomes <- cbind(df_sum_SoC, df_sum_treat, df_sum_inc, df_inc_additional)
    
    l_results <- list(stratified_outcomes = df_stratified_outcomes,
                      weighted_outcomes   = df_weighted_outcomes)
    
  } else if (switch_summary == 0){
    # detailed results - includes summary too
    df_stratified_outcomes <- cbind(df_sum_SoC, df_sum_treat, df_sum_inc, df_inc_additional)
    
    l_results <- list(stratified_outcomes = df_stratified_outcomes,
                      weighted_outcomes   = df_weighted_outcomes,
                      out_SoC   = list(outcomes = l_out_SoC, trace = l_trace_SoC),
                      out_treat = list(outcomes = l_out_treat, trace = l_trace_treat),
                      incremental = l_incremental)
  }
  
  return(l_results)
}

Part_Surv <- function(v_fit_pfs, v_fit_os){
  ## Calcualte state occupancy over time
  # INPUTS
  # v_fit_pfs: fitted curve for pfs
  # v_fit_os: fitted curve for os
  # OUTPUTS:
  # df_trace: a data frame of proportions associated with health state occupation.
  
  prog                 <- v_fit_os - v_fit_pfs # estimate the probability of remaining in the progressed state
  prog[prog < 0]       <- 0                    # in cases where the probability is negative replace with zero
  stable               <- v_fit_pfs            # probability of remaining stable
  dead                 <- 1 - v_fit_os         # probability of being dead
  
  # create output trace
  df_trace <- data.frame(stable = stable, prog = prog, dead = dead)
  return(df_trace)
}


Calculate_Cost_QALY <- function(df_trace, l_c, l_u, v_discount_cost, 
                                v_discount_qaly, switch_vec_cost = 0) {
  ## Calculate the costs and qalys from being in a specific state
  # INPUTS
  # df_trace: data frame with 3 columns determining state occupancy over time
  # v_discount_cost: vector of discount rates for costs
  # v_discount_qaly: vector of discount rates for utilities
  # l_c: list of relevant costs 
  # l_u: list of per-state utilities
  # switch_vec_cost: Switch to determine if costs are vector or scalar inputs
  # OUTPUTS
  # df_results: a data frame containing total discounted costs, life years and QALYs
  #          for a single treatment strategy
  
  
  # create trace matrix
  M_tr <- as.matrix(df_trace)
  
  ## Costs
  if (switch_vec_cost == 0){
    # single cost per state
    v_c  <- M_tr %*% c(l_c$c_pfs, l_c$c_P, l_c$c_D)
  } else if (switch_vec_cost == 1){
    # time dependent costs per state - element wise multiplication
    v_c  <- M_tr * matrix(c(l_c$c_pfs, l_c$c_P, l_c$c_D), ncol = 3)
    # sum over rows to get total cost per time point
    v_c  <- rowSums(v_c) 
  }
  
  # additional testing and AE costs accrued at first time point
  v_c[1] <- v_c[1] + l_c$c_test + l_c$c_AE
  
  ## Life years
  v_ly    <- M_tr %*% c(1,1,0)
  
  ## QALYs
  v_q     <- M_tr %*% c(l_u$u_pfs, l_u$u_P, l_u$u_D)
  
  ## Discounting
  v_c_discount  <- v_c  * v_discount_c                  # Discount costs 
  # LYs and QALYs determined per state so scaled by cycle length
  v_ly_discount <- v_ly * v_discount_q * cycle_length   # Discount LYs
  v_q_discount  <- v_q  * v_discount_q * cycle_length   # Discount QALYS
  
  ## data frame outputs
  df_results <- data.frame( "Discounted_Cost"  = v_c_discount,
                            "Discounted_LYs"   = v_ly_discount,
                            "Discounted_QALYs" = v_q_discount, 
                            check.names = F)
  
  return(df_results)
  
  
}

Calculate_Incremental <- function(df_out_treat, df_out_SoC, t_max, v_times){
  # Calcualte incremental costs and effects, along with the ICER and proportion of
  # QALYs accrued in the extrapolated period
  # INPUTS
  # df_out_treat: dataframe of costs and effects over time for entrectinib
  # df_out_SoC: dataframe of costs and effects over time for SoC
  # v_times: vector of time points where costs and effects have been calculated
  # t_max: time point indicating the end of the observation period for entrectinib
  # OUTPUTS
  # df_incremental: dataframe containing incremental costs and effects, ICER, and proportion 
  #                 of QALYs accrued in the extrapolation period
  
  # Calculate incremental
  
  df_incremental <- data.frame(Inc_Cost  = df_out_treat$Discounted_Cost  - df_out_SoC$Discounted_Cost,
                               Inc_LYs   = df_out_treat$Discounted_LYs   - df_out_SoC$Discounted_LYs,
                               Inc_QALYs = df_out_treat$Discounted_QALYs - df_out_SoC$Discounted_QALYs
                               )
  
  # NMB
  df_incremental$NMB_50000 <- (sum(df_incremental$Inc_QALYs) * 50000) - sum(df_incremental$Inc_Cost)
  df_incremental$NMB_100000 <- (sum(df_incremental$Inc_QALYs) * 100000) - sum(df_incremental$Inc_Cost)
  
  # ICER
  df_incremental$ICER <- sum(df_incremental$Inc_Cost) / sum(df_incremental$Inc_QALYs)
  
  
  ## Effects accrued in extrapolated period
  
  # if no t_max then set to end of run
  if (is.na(t_max)){t_max <- time_horizon - cycle_length}
  # index for end of observation period
  temp_ind <- which.min(abs(v_times - t_max))
  
  # effects accrued in extrapolated period
  df_incremental$Extrap_Inc_QALYs <- sum(df_out_treat$Discounted_QALYs[(temp_ind + 1):length(v_times)]) - 
                                      sum(df_out_SoC$Discounted_QALYs[(temp_ind + 1):length(v_times)]) 
  
  # calculate proportion
  if (sum(df_incremental$Inc_QALYs) < 0 || df_incremental$Extrap_Inc_QALYs[1] < 0){
    # if incremental qalys are negative then SoC dominates
    df_incremental$Proportion_QALYs <- NA
  } else {
    df_incremental$Proportion_QALYs <- df_incremental$Extrap_Inc_QALYs / sum(df_incremental$Inc_QALYs)
  }
  
  return(df_incremental)
  
}


Weighted_Outcomes <- function(df_sum_SoC, df_sum_treat, df_inc_extrap, df_weights){
  ## Calculate pooled results
  # INPUTS
  # df_sum_SoC: data frame of outcome per tumour for SoC
  # df_sum_treat: data frame of outcomes per tumour for the treatment arm
  # df_inc_extrap: data frame of extrapolated outcomes (incremental qalys and proportion)
  # df_weights: data frame of weights per tumour
  # OUTPUTS
  # df_weighted: data frame for weighted costs, qalys and incremental outcomes

  # order by tumour indication
  df_sum_treat <- df_sum_treat[order(row.names(df_sum_treat)),]
  df_sum_SoC   <- df_sum_SoC[order(row.names(df_sum_SoC)),]
  df_weights <- df_weights[order(df_weights$tumour),]

  # set up data frame for weighted outcomes
  df_weighted <- data.frame(Cost_SoC = 0, LY_SoC = 0, QALY_SoC = 0, Cost_treat = 0, 
                            LY_treat = 0, QALY_treat = 0, Inc_Cost = 0, Inc_LYs = 0, 
                            Inc_QALYs = 0, Extrap_Inc_QALYs = 0, 
                            Proportion_QALYs = 0, ICER = 0)

  # weighted function
  weighted_av <- function(outcome, weights){
    sum(outcome * weights, na.rm = TRUE)/sum(weights)
  }


  # individual outcomes
  df_weighted$Cost_treat <- weighted_av(df_sum_treat$Cost_treat, df_weights$Clinical)
  df_weighted$LY_treat   <- weighted_av(df_sum_treat$LY_treat  , df_weights$Clinical)
  df_weighted$QALY_treat <- weighted_av(df_sum_treat$QALY_treat, df_weights$Clinical)
  df_weighted$Cost_SoC   <- weighted_av(df_sum_SoC$Cost_SoC    , df_weights$Clinical)
  df_weighted$LY_SoC     <- weighted_av(df_sum_SoC$LY_SoC      , df_weights$Clinical)
  df_weighted$QALY_SoC   <- weighted_av(df_sum_SoC$QALY_SoC    , df_weights$Clinical)

  # incremental outcomes
  df_weighted$Inc_Cost  <- df_weighted$Cost_treat - df_weighted$Cost_SoC
  df_weighted$Inc_LYs   <- df_weighted$LY_treat   - df_weighted$LY_SoC
  df_weighted$Inc_QALYs <- df_weighted$QALY_treat - df_weighted$QALY_SoC


  ## extrapolated and proportions

  # weighted extrapolated QALYs
  df_inc_extrap <- df_inc_extrap[order(row.names(df_inc_extrap)),]
  df_weighted$Extrap_Inc_QALYs <- weighted_av(df_inc_extrap$Extrap_Inc_QALYs, df_weights$Clinical)

  # proportion of QALYs accrued in extrapolated period
  if (df_weighted$Inc_QALYs <0 ||  df_weighted$Extrap_Inc_QALYs < 0 || is.na(df_weighted$Extrap_Inc_QALYs)){
    df_weighted$Proportion_QALYs <- NA
  } else{
    df_weighted$Proportion_QALYs <- df_weighted$Extrap_Inc_QALYs / df_weighted$Inc_QALY
  }
  df_weighted$NMB_50000 <- (df_weighted$Inc_QALY * 50000) - df_weighted$Inc_Cost
  df_weighted$NMB_100000 <- (df_weighted$Inc_QALY * 100000) - df_weighted$Inc_Cost
  
  df_weighted$ICER   <- df_weighted$Inc_Cost   / df_weighted$Inc_QALY


  return(df_weighted)

}


Deterministic_Frontier <- function(v_tumour_groups, df_outcomes, save_name){
  # Calculate and save frontier
  # INPUTS
  # v_tumour_groups: vector of tumour indications included in the analysis
  # df_outcomes: data frame containing outcomes per tumour
  # save_name: string for filename
  # OUTPUTS
  # df_summary_totals: dataframe with incremental costs and QALYs for each 
  #                    tumour indication combination
  
  # all possible combinations of tumour type
  combs    <- lapply(1:length(v_tumour_groups), function(y){combn(v_tumour_groups,y)})
  
  save_ind <- 0
  totals   <- list()
  for (ind in seq(1,length(combs))){
    for (comb_ind in seq(1, dim(combs[[ind]])[2])){
      # tumour indication in combination list
      comb_sub   <- combs[[ind]][,comb_ind]
      
      save_ind <- save_ind + 1
      # sum incremental
      sum_costs <- sum(df_outcomes[comb_sub,]$Inc_Cost)
      sum_QALYs <- sum(df_outcomes[comb_sub,]$Inc_QALY)
      
      totals[[save_ind]] <- data.frame(Combination = paste(comb_sub, collapse = ", "),
                                       Costs = sum_costs, 
                                       QALYs = sum_QALYs)
    }
    
  }
  
  df_summary_totals <- do.call("rbind",totals)
  # output data frame
  write.csv(df_summary_totals, save_name)
  
  return(df_summary_totals)
  
}