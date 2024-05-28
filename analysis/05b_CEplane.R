# Script to plot CE outcomes and save plot for publication
################################################################################
# Author G.Cupples, E.Krebs. 28/05/2024

# select outcomes
PSA.outputs <- l_outcomes_PSA$weighted_outcomes


# Get % dominated
dom_percent <- round(sum(PSA.outputs$Inc_QALYs < 0) / length(PSA.outputs$Inc_QALYs) * 100)

PSA.outputs %>% mutate(dominated = ifelse(Inc_QALYs>0, "#0482c9", "#acadad")) %>%
  ggplot(aes(x = Inc_QALYs, y = Inc_Cost, color = dominated)) +
  geom_point() +
  annotate(geom = "text", label = paste(dom_percent,"% dominated PSA runs", sep = ""), 
           size = 6, x = -0.5, y = max(PSA.outputs$Inc_Cost), color = "grey40", hjust = 0) +
  scale_color_identity() +
  geom_vline(xintercept = 0, linewidth = 1) +
  geom_hline(yintercept = 0, linewidth = 1) +
  xlim(-0.5, 0.5) +
  geom_abline(intercept = 0, slope = 100000, color = "#662506", lwd = 1.2) +
  annotate(geom = "text", label = "ICER: $100,000/QALY", size = 5,
           x = .21, y = 120000, color = "#662506", hjust = 0) +
  geom_abline(intercept = 0, slope = 50000, color = "#EA6E13", lwd = 1.2) +
  annotate(geom = "text", label = "ICER: $50,000/QALY", size = 5,
           x = .21, y = 114000, color = "#EA6E13", hjust = 0) +
  xlab("Incremental quality-adjusted life years (QALYs)") +
  scale_y_continuous(labels = scales::dollar_format(), limits = c(0, max(PSA.outputs$Inc_Cost)*1.15)) +
  ylab("Incremental costs (2021 CAD)") +
  # ggtitle("Cost-effectiveness plane") +
  theme_minimal() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 15))

## save plot

# check folder exists
dir_PSA <- file.path("figs", "PSA")
if (!dir.exists(dir_PSA)) dir.create(dir_PSA, recursive = TRUE)

# save
savename <- paste("figs/PSA/CEplane_observed_", switch_observed, "_test_", switch_test,
                  ".png", sep = "")

ggsave(filename = savename,
       plot = last_plot(),
       width = 8, height = 6, units = "in", scale = 1,
       dpi = 300)
