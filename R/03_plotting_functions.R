# Function file containing plotting functions
###############################################################################
# Author G. Cupples, E.Krebs. 08/09/2023

Create_Plots <- function(df_curves, df_trace_treat, df_trace_SoC, v_times, tumour, 
                         t_max = NULL){
  # Create all deterministic plots
  # INPUTS
  # df_curves: data frame containing points defining survival curve
  # df_trace_treat: dataframe with 3 columns describing state occupancy for treatment arm
  # df_trace_SoC: dataframe with 3 columns describing state occupancy for SoC
  # v_times: vector of time points for each cycle in the model
  # tumour: tumour indication to create plots for
  # t_max: maximum follow-up time in observation period
  
  ## Save names
    save_exponential <- paste("figs/curve_fit/observed_", switch_observed, 
                              "_test_", switch_test, "_tumour_",tumour,".png", sep = "")
    save_trace       <- paste("figs/trace/observed_", switch_observed, "_test_",
                              switch_test, "_tumour_",tumour,".png", sep = "")
  
  
  if (switch_observed == 0){
    # plot without t_max as no switch at end of observation period
    Plot_Exponential(df_curves, tumour, v_times, save_exponential)
    Plot_Trace(df_trace_treat, df_trace_SoC, tumour, v_times, save_trace)
  } else if (switch_observed == 1){
    # plot including end of observation period
    Plot_Exponential(df_curves, tumour, v_times, save_exponential, 
                     t_max)
    Plot_Trace(df_trace_treat, df_trace_SoC, tumour, v_times, save_trace, 
               t_max)
  }
  

}


Plot_Exponential <- function(l_curve, tumour, v_times, savename, t_max = ""){
  ## plot exponential curves constructed from median survival data
  # INPUTS
  # l_curve: dataframe containing curve points for entrectinib and SoC
  # tumour: tumour indication to plot
  # v_times: vector of time points where the curves have been calculated
  # savename: save name for plot
  # t_max: end of observation period
  
  png(file=savename, width=800, height=800, res = 90)
  par(cex = 1.5)
  # entrectinib curves
  plot(v_times, l_curve$os_treat, type = 'l', lwd = 3, col = "green",
       ylab = paste("Survival Probability", tumour, sep = " "), xlab = "Time",
         main = paste("Exponential curves,", tumour, sep = " "))
  lines(v_times, l_curve$pfs_treat,          lwd = 3, col = "red"  )
  # SoC curves
  lines(v_times, l_curve$os_SoC,  lty = 2, lwd = 3, col = "green")
  lines(v_times, l_curve$pfs_SoC, lty = 2, lwd = 3, col = "red"  )
  # observed point
  if (switch_observed == 1){
    abline(v = t_max, lwd = 3, col = "grey")
  }
  # legend
  legend("right", c("OS", "PFS", "entrectinib", "SoC"), 
         col = c("green", "red", "black", "black"), lty=c(1,1,1,2), lwd = 3, bty='n')
  dev.off()
  
  
}


Plot_Trace <- function(df_trace_treat, df_trace_SoC, tumour, v_times, savename, t_max = ""){
  ## Plot trace determined by PartSA
  # INPUTS
    # df_trace_treat: dataframe with 3 columns desctribing state occupancy for treatment arm
    # df_trace_SoC: dataframe with 3 columns desctribing state occupancy for SoC
    # tumour: tumour group to plot
    # v_times: vector of time points that the trace has been calculated at
    # savename: plot name
    # t_max: end of observation point
  
  # trace as matrix
  M_treat <- as.matrix(df_trace_treat)
  M_SoC   <- as.matrix(df_trace_SoC)
  
  png(file=savename, width=800, height=800, res = 90)
  par(cex = 1.5)
  # entrectinib
  matplot(M_treat, type = 'l', lty=1, lwd = 3, col = c(3,2,1), xlab = "Time", 
          ylab = paste("Proportion",tumour, sep = " "))
  # SoC
  lines(M_SoC[,"stable"], lty=2, lwd = 3, col = 3)
  lines(M_SoC[,"prog"  ], lty=2, lwd = 3, col = 2)
  lines(M_SoC[,"dead"  ], lty=2, lwd = 3, col = 1)
  # observed point
  if (switch_observed == 1){
    ind_switch <- which.min(abs(v_times - t_max))
    abline(v = ind_switch, lwd = 3, col = "grey")
  }
  # legend
  legend("right", c(v_state_names, v_treat_names), col=c(3,2,1,1,1), lty=c(rep(1,n_states),2,1), lwd = 3, bty='n')
  dev.off()
  
  
}

