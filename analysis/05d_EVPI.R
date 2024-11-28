# Script to generate stratified and pooled expected value of perfect information (EVPI)

## Prep data

# get strat obj
PSA.outputs <- l_outcomes_PSA$stratified_obj

# Replace NAs
for (s in 1:length(PSA.outputs)) {
  PSA.outputs[[s]]$cost[is.na(PSA.outputs[[s]]$cost)] <- 0
  PSA.outputs[[s]]$effectiveness[is.na(PSA.outputs[[s]]$effectiveness)] <- 0
}

# set willingness to pay (WTP) limits 
wtp <- seq(0, 500000, by = 10000)

# check folder exists
dir_PSA <- file.path("figs", "EVPI")
if (!dir.exists(dir_PSA)) dir.create(dir_PSA, recursive = TRUE)



## EVPI - stratified
# some tumour indications only accrue costs (not QALYs) and so are not included 
# in the frontier
tumours <- c("Breast","CRC","MASC","NSCLC","Sarcoma","Thyroid")

# create plot
savename <- paste('figs/EVPI/stratified_observed_', switch_observed, "_test_", 
                              switch_test, ".png", sep = "")
  
png(filename = savename, width = 15, height = 10, units = "in", res = 300)
par(mfrow = c(2,3))
par(cex = 1.3)
x_ax_points <- wtp[seq(1, length(wtp), by = 10)]
x_ax_labels <- wtp[seq(1, length(wtp), by = 10)] / 1000
y_ax_points <- c(0, 50000, 100000, 150000, 200000)
y_ax_labels <- y_ax_points / 1000
temp <- data.frame(WTP = wtp)
i <- 0
for (tumour in tumours){
  i <- i+1
  evpi_obj <- calc_evpi(PSA.outputs[[tumour]], wtp = wtp, pop = 1)
  evpi_title <- tumour


  plot.new()
  grid(nx = NULL, ny = NULL,
       lty = 1,      # Grid line type
       col = "gray85", # Grid line color
       lwd = 2)      # Grid line width
  par(new = TRUE)
  plot(evpi_obj$WTP, evpi_obj$EVPI, type = "l", lwd = 2.5, ylim = c(0, 200000),
       xlim = c(0, 500000), xlab = "WTP (Thousand $/QALY)", xaxt='n',
       ylab = "EVPI (Thousand $)", yaxt = 'n', cex = 1.5,
       col = "#006bb6")
  title(tumour, line = 0.5)
  axis(1, at = x_ax_points, labels = x_ax_labels)
  axis(2, at = y_ax_points, labels = y_ax_labels, las = 1)

}
dev.off()

## EVPI - pooled
evpi_df <- calc_evpi(PSA.outputs$Breast, wtp = wtp, pop = 1)
colnames(evpi_df) <- c("WTP", "Breast")
evpi_df <- cbind(evpi_df, CRC = calc_evpi(PSA.outputs$CRC, wtp = wtp, pop = 1)$EVPI)
evpi_df <- cbind(evpi_df, MASC = calc_evpi(PSA.outputs$MASC, wtp = wtp, pop = 1)$EVPI)
evpi_df <- cbind(evpi_df, Neuroendocrine = calc_evpi(PSA.outputs$Neuroendocrine, wtp = wtp, pop = 1)$EVPI)
evpi_df <- cbind(evpi_df, NSCLC = calc_evpi(PSA.outputs$NSCLC, wtp = wtp, pop = 1)$EVPI)
evpi_df <- cbind(evpi_df, Other = calc_evpi(PSA.outputs$Other, wtp = wtp, pop = 1)$EVPI)
evpi_df <- cbind(evpi_df, Pancreatic = calc_evpi(PSA.outputs$Pancreatic, wtp = wtp, pop = 1)$EVPI)
evpi_df <- cbind(evpi_df, Sarcoma = calc_evpi(PSA.outputs$Sarcoma, wtp = wtp, pop = 1)$EVPI)
evpi_df <- cbind(evpi_df, Thyroid = calc_evpi(PSA.outputs$Thyroid, wtp = wtp, pop = 1)$EVPI)

weighted_evpi <- evpi_df
for (s in 2:ncol(evpi_df)) {
  weighted_evpi[, s] <- (evpi_df[, s] * df_weights$Clinical[s - 1]) / sum(df_weights$Clinical)
}
weighted_evpi$row_sums <- rowSums(weighted_evpi[-1])

weighted_evpi_sum <- evpi_obj
weighted_evpi_sum$EVPI <- weighted_evpi$row_sums



# create plot
y_ax_points <- c(0, 10000, 20000, 30000, 40000, 50000)
y_ax_labels <- c(0, '10', '20', '30', '40', '50')

savename <- paste('figs/EVPI/weighted_observed_', switch_observed, "_test_", 
                  switch_test, ".png", sep = "")

png(filename = savename, width = 8, height = 8, units = "in", 
    res = 300)
par(mfrow = c(1,1))
par(cex = 1.3)
plot.new()
grid(nx = NULL, ny = NULL,
     lty = 1,      # Grid line type
     col = "gray85", # Grid line colorx
     lwd = 2)      # Grid line width
par(new = TRUE)
plot(weighted_evpi_sum$WTP, weighted_evpi_sum$EVPI, ylim = c(0, 50000), type = "l", lwd = 2,
     xlim = c(0, 500000), xlab = "WTP (Thousand $/QALY)", xaxt='n', yaxt='n',
     ylab = "EVPI (Thousand $)", main = "EVPI: tumour-agnostic", cex = 1.5, las = 1)
axis(1, at = x_ax_points, labels = x_ax_labels)
axis(2, at = y_ax_points, labels = y_ax_labels, las = 1)
dev.off()
