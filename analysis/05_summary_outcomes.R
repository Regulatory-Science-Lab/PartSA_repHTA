# Probabilistic summary outcomes for publication
################################################################################
# Author G.Cupples, E.Krebs. 12/03/2024

## Functions ##

Summary_Outcomes <- function(v_psa_outcome){
  # create summary outcomes for PSA
  # INPUTS
  # v_psa_outcome: vector of outcome for summary, obtained from PSA
  # OUTPUTS
  # df_summary: dataframe containing summary information
  
  # mean and 95%CI
  mean_out  <- mean(v_psa_outcome)
  sd_out    <- sd(v_psa_outcome)
  
  ci_l <- quantile(v_psa_outcome, c(0.025, 0.975))[[1]]
  ci_u <- quantile(v_psa_outcome, c(0.025, 0.975))[[2]]
  
  # median and IQR
  median_out <- median(v_psa_outcome)
  iqr_l      <- quantile(v_psa_outcome, c(0.25))[[1]]
  iqr_u      <- quantile(v_psa_outcome, c(0.75))[[1]]
  
  # create dataframe of outcomes
  df_summary <- data.frame(mean = mean_out, sd = sd_out, low_95 = ci_l, up_95 = ci_u, 
                           median = median_out, low_iqr = iqr_l, up_iqr = iqr_u)
  return(df_summary)
  
}

Prop_Dom <- function(psa_out){
  # get proportion of negative runs from ICER
  prop_dominated <- length(psa_out$ICER[psa_out$ICER < 0 & psa_out$Inc_Cost > 0])/length(psa_out$ICER) * 100
  return(prop_dominated)
}

Prop_CE <- function(ICER, WTP){
  # get proportion of runs that are Cost-effective at certain willingness-to-pay
  # threshold
  prop_CE <- length(ICER[ICER < WTP & ICER > 0])/ length(ICER) * 100
  return(prop_CE)
}



## Pooled ##

# columns to find summary for
sum_cols <- c("Inc_Cost","Inc_QALYs","ICER","NMB_100000", "NMB_200000")
# get summary statistics for all included columns
df_pooled_summary <- sapply(l_outcomes_PSA$weighted_outcomes[,sum_cols], Summary_Outcomes)

# calculation proportion of dominated runs
prop_dominated <- Prop_Dom(l_outcomes_PSA$weighted_outcomes)
df_pooled_summary <- cbind(df_pooled_summary, prop_dominated)

# proportion of Cost-effective runs at WTP $100,000 and $200,000
prop_ce_1 <- Prop_CE(l_outcomes_PSA$weighted_outcomes$ICER, WTP = 100000)
prop_ce_2 <- Prop_CE(l_outcomes_PSA$weighted_outcomes$ICER, WTP = 200000)
df_pooled_summary <- cbind(df_pooled_summary, prop_ce_1)
df_pooled_summary <- cbind(df_pooled_summary, prop_ce_2)

## Stratified ##

# mean, 95% CI, median, IQR - Inc Cost, Inc QALY, ICER
l_strat_summary <- list()
for (tumour in v_tumour){
  # columns to find summary for
  sum_cols <- c("Inc_Cost","Inc_QALYs","ICER","NMB_100000", "NMB_200000")
  # get summary statistics for all included columns
  df_summary <- sapply(l_outcomes_PSA$stratified_outcomes[[tumour]][,sum_cols], Summary_Outcomes)
  
  # calculation proportion of dominated runs
  prop_dom_ts <- Prop_Dom(l_outcomes_PSA$stratified_outcomes[[tumour]])
  df_summary <- cbind(df_summary, prop_dom_ts)
  
  # proportion of Cost-effective runs at WTP $100,000 and $200,000
  prop_ce_1 <- Prop_CE(l_outcomes_PSA$stratified_outcomes[[tumour]]$ICER, WTP = 100000)
  prop_ce_2 <- Prop_CE(l_outcomes_PSA$stratified_outcomes[[tumour]]$ICER, WTP = 200000)
  df_summary <- cbind(df_summary, prop_ce_1)
  df_summary <- cbind(df_summary, prop_ce_2)

  l_strat_summary[[tumour]] <- df_summary
}
