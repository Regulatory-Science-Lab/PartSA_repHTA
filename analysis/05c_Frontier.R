options(scipen = 999)

totals <- df_frontier_PSA %>% 
  mutate(model = "PartSA") %>%
  add_column(number = 1:nrow(df_frontier_PSA), .before = 1)


# determine frontier
frontier.df <- Frontier_Combinations(totals)
frontier <- totals %>% 
  filter(number %in% frontier.df$combination) %>%
  arrange(Costs) %>%
  add_column(ICER_front = frontier.df$ICER_front, .after = 6) %>%
  mutate(labels = paste(Combination, ", ", "$",
                        format(round(as.numeric(ICER_front), 0), nsmall=0, big.mark=","),
                        "/QALY", sep = ""))


# create plot
pf <- ggplot(totals, aes(Costs, QALYs)) +
  geom_point(show.legend = FALSE, size = 0.9, color = "#00A08A", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey50", linewidth = 0.7) +
  geom_point(data = frontier, aes(Costs, QALYs), colour = "#00A08A", size = 3, alpha = 1) +
  geom_segment(data = frontier, aes(x = 0, y = 0, xend = Costs[1], yend = QALYs[1])) +
  geom_segment(data = frontier, aes(x = Costs[1], y = QALYs[1],
                                    xend = Costs[2], yend = QALYs[2])) +
  geom_segment(data = frontier, aes(x = Costs[2], y = QALYs[2],
                                    xend = Costs[3], yend = QALYs[3])) +
  geom_segment(data = frontier, aes(x = Costs[3], y = QALYs[3],
                                    xend = Costs[4], yend = QALYs[4])) +
  geom_segment(data = frontier, aes(x = Costs[4], y = QALYs[4],
                                    xend = Costs[5], yend = QALYs[5])) +
  # tumour indications
  geom_text(data = frontier[c(1), ], aes(Costs, QALYs, label = Combination),
            hjust = 0.5, vjust = -0.8, colour = "#006d5e") +
  geom_text(data = frontier[c(2), ], aes(Costs, QALYs, label = Combination),
            hjust = 0.75, vjust = -0.6, colour = "#006d5e") +
  geom_text(data = frontier[c(3), ], aes(Costs, QALYs, label = Combination),
            hjust = 1.05, vjust = -0.1, colour = "#006d5e") +
  geom_text(data = frontier[c(4), ], aes(Costs, QALYs, label = Combination),
            hjust = 1.02, vjust = -0.8, colour = "#006d5e") +
  geom_text(data = frontier[c(5), ], aes(Costs, QALYs, label = Combination),
            hjust = 1.02, vjust = -0.8, colour = "#006d5e") +
  # ICERs
  geom_text(data = frontier[1:3, ], 
            aes(Costs, QALYs, label = scales::label_dollar(suffix = "/QALY")(frontier[1:3, ]$ICER_front)),
            hjust = 0, vjust = 1.6, colour = "black", fontface = "bold") +
  geom_text(data = frontier[4:5, ],
            aes(Costs, QALYs, label = scales::label_dollar(suffix = "/QALY")(frontier[4:5, ]$ICER_front)),
            hjust = -0.05, vjust = 0.7, colour = "black", fontface = "bold") +
  scale_x_continuous(breaks = pretty_breaks(n=4),
                     labels = scales::dollar_format(),
                     limit = c(0, max(totals$Costs) * 1.15)) +
  scale_y_continuous(breaks = pretty_breaks(),
                     limit = c(0 , max(totals$QALYs) * 1.15),
                     expand = c(0, 0)) +
  xlab("Incremental Costs (2021 CAD)") +
  ylab("Incremental QALYs") +
  ggtitle("") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.3,1,0.3,0.3), "cm"))
pf


filename <- paste("figs/frontier/PSA_observed_", switch_observed, "_test_", switch_test,
                  ".png", sep = "")

ggsave(file = filename,
       pf,
       width = 8, height = 6, units = "cm", scale = 1.8,
       dpi = 300)

