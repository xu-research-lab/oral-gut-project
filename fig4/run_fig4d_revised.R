library(ggplot2)
library(dplyr)
library(tidyr)
load("fig4d.RData")
p_data <- Species_shared_trans_diff %>%
  pivot_longer(
    cols = shared_rate:tran_rate,
    names_to = "rate_type",
    values_to = "rate"
  )
p_data <- p_data %>% mutate(
  rate_type = factor(ifelse(rate_type == "shared_rate", "shared rate", "Trans rate"),
                     levels = c("shared rate", "Trans rate")),
  sig = ifelse(adjusted_pvalue_shared <= 0.25 | adjusted_pvalue_trans <= 0.25, "*", NA),
  sig = ifelse(is.na(adjusted_pvalue_trans) & rate_type == "Trans rate", NA, sig))
p <- ggplot(subset(p_data, group == "non-PMA")) +
  geom_point(aes(x = rate, y = id, color = Response_6), size = 2.4) +
  geom_line(aes(group = interaction(SGB, rate_type), x = rate, y = id), color = "gray") +
  geom_text(data = subset(p_data, rate_type == "shared rate"), aes(x = 0.45, y = id, label = sig, vjust = 0.7), size = 6, color = "black") +
  geom_text(data = subset(p_data, rate_type == "Trans rate"), aes(x = 0.45, y = id, label = sig, vjust = 0.7), size = 6, color = "black") +
  scale_color_manual(name = "Response", values = c("R" = "#377EB8", "NR" = "#E41A1C")) +
  scale_y_discrete(position = "right") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4)) +
  facet_grid(cols = vars(rate_type)) +
  theme_bw(base_size = 12) +
  xlab(NULL) + ylab(NULL) +
  expand_limits(x = 0) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 11),
        plot.margin = margin(2, 2, 2, 20),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.4),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        legend.title = element_text(size = 11), legend.text = element_text(size = 11),
        legend.key.size = unit(0.5, "cm"))
ggsave("fig4_adjusted/fig4d_revised.pdf", p, height = 2.5, width = 9)
ggsave("fig4_adjusted/fig4d_revised.png", p, height = 2.5, width = 9, dpi = 300)
cat("done fig4d_revised\n")
