# Script to run a probabilistic tumour-agnostic partitioned survival analysis 
# (PartSA) with a variety of input survival data
################################################################################
# Author G.Cupples, E.Krebs. 10/10/2023

# run analysis
l_outcomes_PSA <- run_PSA(df_data_treat, df_data_SoC, n_psa, seed, switch_observed)

# Create output file directory if it does not yet exist
dir_name <- file.path("outputs")
if (!dir.exists(dir_name)) dir.create(dir_name)

# save outputs
psa_save <- paste("prob_observed_", switch_observed, 
                  "_test_", switch_test,".rds", sep = "")
save(l_outcomes_PSA, file = file.path(dir_name, psa_save))

# comb PSA obj
PSA_obj_weighted <- make_psa_obj(cost          = as.data.frame(l_outcomes_PSA$weighted_outcomes[,c(1,4)]),
                                 effectiveness = as.data.frame(l_outcomes_PSA$weighted_outcomes[,c(3,6)]),
                                 strategies    = v_treat_names)

# Frontier
save_name_frontier <- paste("frontier_PSA_observed_", switch_observed, 
                            "_test_", switch_test,".csv", sep = "")
df_frontier_PSA <- PSA_frontier(v_tumour, l_outcomes_PSA, 
                                file.path(dir_name, save_name_frontier))


