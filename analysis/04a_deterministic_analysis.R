# Script to run a deterministic # tumour-agnostic partitioned survival analysis 
# (PartSA) with a variety of input survival data
################################################################################
# Author G.Cupples, E.Krebs. 10/10/2023

# run analysis
l_outcomes <- run_PartSA(df_data_treat, l_curves, l_c_treat, l_c_SoC, l_u_treat, 
                         l_u_SoC, v_times, df_weights, switch_observed,
                         switch_summary, switch_plot)

# Create output file directory if it does not yet exist
dir_name <- file.path("outputs")
if (!dir.exists(dir_name)) dir.create(dir_name)

# save output
det_save <- paste("det_observed_", switch_observed, "_test_", 
                  switch_test,".rds", sep = "")
save(l_outcomes, file = file.path(dir_name, det_save))

# frontier
save_name_det <- paste("frontier_det_observed_", switch_observed, 
                       "_test_", switch_test, ".csv", sep = "")
df_frontier_det <- Deterministic_Frontier(v_tumour, l_outcomes$stratified_outcomes, 
                                          file.path(dir_name,save_name_det))
