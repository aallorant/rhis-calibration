## =====================================================================
## fig_transport_provinces.R -- the per-province transport test statistic
## (S^2_p - tau_p) / sd(S^2_p), standard normal under transport, one-sided: one point per
## province, labelled with the number of districts compared. Mozambique is the boundary-dropped
## run on design-based inputs. Drawn at 6.30in so \includegraphics[width=\linewidth] applies
## no scaling. Run after R/transport_diagnostic.R.
## Reads:  output/canonical/estimates/transport_diagnostic_by_province.csv
## Writes: figures/fig_transport_provinces.pdf, output/canonical/estimates/fig_transport_provinces.csv
## Usage:  Rscript R/fig_transport_provinces.R
## =====================================================================
suppressMessages({library(dplyr); library(ggplot2)})
fix_uesc <- function(x) {
  m <- gregexpr("<U\\+([0-9A-Fa-f]{4})>", x)
  regmatches(x, m) <- lapply(regmatches(x, m), function(s)
    if (length(s)) intToUtf8(strtoi(sub("<U\\+([0-9A-Fa-f]{4})>", "\\1", s), 16L), multiple = TRUE) else s)
  x
}
per <- read.csv("output/canonical/estimates/transport_diagnostic_by_province.csv")
per <- per[!(per$country == "Mozambique" & !per$boundary_dropped), ]
crit <- qnorm(0.95)
d <- per |>
  mutate(country = factor(country, levels = c("Mozambique", "Zimbabwe", "Malawi")),
         prov = fix_uesc(prov),
         lab = sprintf("%s (%d)", prov, n_dist),
         result = factor(ifelse(p_one_sided < 0.05, "Rejects transport", "Does not reject"),
                         levels = c("Does not reject", "Rejects transport"))) |>
  arrange(country, zstat) |>
  mutate(lab = factor(lab, levels = unique(lab)))
stopifnot(nrow(d) == 23)
g <- ggplot(d, aes(zstat, lab)) +
  geom_vline(xintercept = 0, colour = "grey75") +
  geom_vline(aes(xintercept = crit, linetype = "One-sided 5% critical value"), colour = "grey35") +
  geom_point(aes(shape = result), size = 2.2, colour = "grey15") +
  facet_grid(country ~ ., scales = "free_y", space = "free_y") +
  scale_shape_manual(values = c(`Does not reject` = 1, `Rejects transport` = 16), name = NULL, drop = FALSE) +
  scale_linetype_manual(values = c(`One-sided 5% critical value` = "22"), name = NULL) +
  scale_x_continuous(breaks = seq(-1, 2, 1), limits = c(-1.3, 2.1)) +
  labs(x = "Test statistic for each province, standard normal under transport", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = "grey93"),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
        axis.text.y = element_text(size = 8.5), axis.title.x = element_text(size = 10),
        legend.position = "none")    # keys are described in the figure caption
dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig_transport_provinces.pdf", g, width = 6.3, height = 5.2)
write.csv(d |> transmute(country, prov, n_dist, S2, tau, t, zstat, p_one_sided, result),
          "output/canonical/estimates/fig_transport_provinces.csv", row.names = FALSE)
cat(sprintf("provinces %d | rejecting at one-sided 0.05: %d (%s) | z range %.2f to %.2f\n",
            nrow(d), sum(d$p_one_sided < 0.05), paste(d$prov[d$p_one_sided < 0.05], collapse = ", "),
            min(d$zstat), max(d$zstat)))
