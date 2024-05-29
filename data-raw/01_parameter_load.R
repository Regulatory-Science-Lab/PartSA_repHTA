# Script to load parameters related to costs, and effects for the treatment arm
# and the comparator
###############################################################################
# Author G. Cupples, E. Krebs. 10/10/2023



## load data
# Publicly available data set
filename <- "data-raw/01_public_parameters.xlsx"

# assign each sheet in parameter file to dataframe entry
sheets <- excel_sheets(filename)
l_parameters <- lapply(sheets, function(X) read_excel(filename, sheet = X, .name_repair = "minimal"))
# name the sheets with corresponding sheet name
names(l_parameters) <- sheets

## Workspace
Workspace <- l_parameters$Workspace
list2env(setNames(as.list(Workspace$value), Workspace$Parameter), .GlobalEnv)

Workspace_names <- l_parameters$Workspace_names
for (i in 1:length(Workspace_names$Parameter)){
  temp <- strsplit(Workspace_names[[i,2]], ", ")
  list2env(setNames(temp, Workspace_names[[i,1]]),.GlobalEnv)
}

n_states <- length(v_state_names) # number of states

v_times <- seq(from = 0, to = time_horizon, by = cycle_length) # vector of time points
# calculate discount weights vector
v_discount_c <- 1 / (1 + discount_cost   * cycle_length) ^ (v_times / cycle_length) # cost
v_discount_q <- 1 / (1 + discount_effect * cycle_length) ^ (v_times / cycle_length) # effect


## tumour groups
v_tumour <- l_parameters$`3.2_MedianOSPFS_treat`$tumour

## weights
# set data.frame of weights to match tumour groups included in analysis
df_weights <- data.frame(tumour   = v_tumour, 
                         Clinical = rep(0,length(v_tumour)), 
                         Trial    = rep(0,length(v_tumour)))

# calculate prevalence for all tumour indications
df_trial_weight <- l_parameters$`3.5_TrialWeights`

# prevalence taken from CADTH reanalysis
df_clinical_presentation <- l_parameters$`3.4_CADTHWeights`


df_clinical_presentation$counted <- "N"
for (tumour in v_tumour[v_tumour != "Other"]){
  ind_tumour <- grep(tumour, df_clinical_presentation$tumour, ignore.case = TRUE)
  # update counted to Y
  df_clinical_presentation[ind_tumour,]$counted <- "Y"
  # sum and store in stratified dataframe - ignoring NA's from summary
  df_weights[df_weights$tumour == tumour,]$Clinical <- sum(df_clinical_presentation[ind_tumour,]$value, na.rm = TRUE)
  # trial weights 
  df_weights[df_weights$tumour == tumour,]$Trial    <- df_trial_weight[df_trial_weight$tumour == tumour,]$value
  
}
# group 'other' tumour groups together
other_Clinical <- sum(df_clinical_presentation[df_clinical_presentation$counted =="N",]$value, na.rm = TRUE)
df_weights[df_weights$tumour == "Other",]$Clinical <- other_Clinical
df_weights[df_weights$tumour == "Other",]$Trial    <- df_trial_weight[df_trial_weight$tumour == "Other",]$value



### Other Parameters ###
weights_no_other <- df_weights[df_weights$tumour != "Other",]
weights_no_other <- weights_no_other[order(weights_no_other$tumour),]

## Comparator treatment
df_c_pack_SoC <- l_parameters$`1.3_TreatmentCost`
df_c_pack_SoC <- df_c_pack_SoC[order(df_c_pack_SoC$tumour),]
# calculate the average of all tumour groups with costs
c_treat_av <- sum(df_c_pack_SoC$value * weights_no_other$Clinical) / sum(weights_no_other$Clinical)
# assign average to "Other"
df_c_pack_SoC <- rbind(df_c_pack_SoC, 
                       c("Other", c_treat_av, "fixed", c_treat_av, "na", "na", NA))
df_c_pack_SoC$value <- as.numeric(df_c_pack_SoC$value)


## Comparator administration cost
df_c_admin_SoC <- l_parameters$`1.4_AdminCost`
df_c_admin_SoC <- df_c_admin_SoC[order(df_c_admin_SoC$tumour),]
# assign costs to 'Other' tumour groups
c_admin_av     <- sum(df_c_admin_SoC$value * weights_no_other$Clinical) / sum(weights_no_other$Clinical)
df_c_admin_SoC <- rbind(df_c_admin_SoC, 
                        c("Other", c_admin_av, "uniform", c_admin_av*0.8, c_admin_av*1.2, "na", NA))
df_c_admin_SoC$value <- as.numeric(df_c_admin_SoC$value)


## Non-cancer care
df_c_noncancer    <- l_parameters$`1.2_NonCancerCost`
df_c_noncancer    <- df_c_noncancer[order(df_c_noncancer$tumour),]
# calculate the average of all tumour groups with costs
c_noncancer_av <- sum(df_c_noncancer$value * weights_no_other$Clinical) / sum(weights_no_other$Clinical)
# assign average to "Other"
df_c_noncancer <- rbind(df_c_noncancer, 
                        c("Other", c_noncancer_av, "uniform", c_noncancer_av*0.8, c_noncancer_av*1.2, "na", NA))
df_c_noncancer$value <- as.numeric(df_c_noncancer$value)

# Tumour agnostic costs
df_costs_agnostic <- l_parameters$`1.1_TumourAgnosticCosts`
# extract as individual parameters
list2env(setNames(as.list(df_costs_agnostic$value), paste(df_costs_agnostic$parameter,'fixed', sep = "_")), .GlobalEnv)

# Time to off treatment
ttot_fixed <- l_parameters$`3.6_TTOT`$value


## Testing
if (switch_test == 0){
  # no entrectinib testing costs
  df_c_test_treat <- data.frame()
  df_c_test_treat <- rbind(df_c_test_treat, rep(0,length(v_tumour)))
  colnames(df_c_test_treat) <- v_tumour
} else if (switch_test == 1){
  # entrectinib testing costs
  df_test_temp <- l_parameters$`1.5_TestingCost`
  
  df_test_temp$counted <- "N"
  
  v_c_test_treat <- c()
  for (tumour in v_tumour[v_tumour != "Other"]){
    # assign costs
    v_c_test_treat <- cbind(v_c_test_treat, df_test_temp[df_test_temp$tumour == tumour,]$value)
    # check off counted
    df_test_temp$counted[df_test_temp$tumour == tumour] <- "Y"
  }
  # other costs are weighted average of remaining tumour indications
  other_cost <- sum(df_test_temp$value[df_test_temp$counted == "N"]  * 
                      df_clinical_presentation$value[df_test_temp$counted == "N"], na.rm = TRUE)/
    sum(df_clinical_presentation$value[df_test_temp$counted =="N"], na.rm = TRUE)
  
  
  v_c_test_treat <- cbind(v_c_test_treat, other_cost)
  df_c_test_treat <- as.data.frame(v_c_test_treat)
  colnames(df_c_test_treat) <- v_tumour
  
  # set the excluded tumour groups costs to 0
  if (paste(v_test_exclude, collapse = "|") != "") {
    df_c_test_treat[,grep(paste(v_test_exclude, collapse = "|"), x = colnames(df_c_test_treat))] <- 0
  }
  
  
}


## Utilities
df_utilities <- l_parameters$`2.1_Utilities`
# extract as individual parameters
list2env(setNames(as.list(df_utilities$value), paste(df_utilities$parameter, 'fixed', sep = "_")), .GlobalEnv)


## Cost per state
# costs for death state the same for each tumour group
c_D_SoC   <- c_D_fixed
c_D_treat <- c_D_fixed

l_c_SoC   <- list()
l_c_treat <- list()
l_u_SoC   <- list()
l_u_treat <- list()
for (tumour in v_tumour){
  ## SoC
  
  # progression-free
  c_pfs_SoC <- df_c_noncancer[df_c_noncancer$tumour == tumour,]$value +
               df_c_pack_SoC[df_c_pack_SoC$tumour == tumour,]$value + 
               df_c_admin_SoC[df_c_admin_SoC$tumour == tumour,]$value + c_care_pfs_SoC_fixed
  # progressed
  c_P_SoC   <- df_c_noncancer[df_c_noncancer$tumour == tumour,]$value +
                c_care_P_fixed
  
  # create list entry
  l_c_SoC[[tumour]] <- list(
    c_pfs  = c_pfs_SoC,
    c_P    = c_P_SoC,
    c_D    = c_D_SoC,
    c_AE   = c_AE_SoC_fixed,
    c_test = c_test_SoC_fixed # SoC patients were not tested - unknown NTRK status
  )
  
  ## entrectinib
  
  # vector costs
  c_pfs_treat  <- df_c_noncancer[df_c_noncancer$tumour == tumour,]$value + 
    c_pack_treat_fixed + c_admin_treat_fixed + c_care_pfs_treat_fixed
  
  # progression-free
  c_pfs_treat <- rep(c_pfs_treat, length(v_times))
  
  # treatment costs up until time to off treatment (11 months)
  c_pack_treat_P <- rep(0, length(v_times))
  c_pack_treat_P[1:which.min(abs(v_times - ttot_fixed))] <- c_pack_treat_fixed
  # treatment administration costs up until time to off treatment (11 months)
  c_admin_treat_P <- rep(0, length(v_times))
  c_admin_treat_P[1:which.min(abs(v_times - ttot_fixed))] <- c_admin_treat_fixed
  # progressed
  c_P_treat <- df_c_noncancer[df_c_noncancer$tumour == tumour,]$value +
    c_pack_treat_P + c_admin_treat_P + c_care_P_fixed
  
  # death
  c_D_treat <- rep(0,length(v_times))


  # create list entry
  l_c_treat[[tumour]] <- list(
    c_pfs  = c_pfs_treat,
    c_P    = c_P_treat,
    c_D    = c_D_treat,
    c_AE   = c_AE_treat_fixed,
    c_test = df_c_test_treat[[tumour]]
  )
  
  
  # Utilities - not currently varying by tumour but framework set up to enable this
  l_u_treat[[tumour]] <- list(
    u_pfs = u_pfs_treat_fixed,
    u_P   = u_P_treat_fixed,
    u_D   = u_D_fixed
  )
  l_u_SoC[[tumour]] <- list(
    u_pfs = u_pfs_SoC_fixed,
    u_P   = u_P_SoC_fixed,
    u_D   = u_D_fixed
  )

}

