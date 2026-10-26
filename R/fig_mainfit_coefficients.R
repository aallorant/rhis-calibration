## =====================================================================
## fig_mainfit_coefficients.R -- four-country forest plot of the main fit (block-covariance
## EBLUP at area x age x sex): province offsets (age-sex adjusted), national age-sex offsets
## (EBLUP rows) and a sigma_u^2 panel with 95% profile intervals.
## Drawn at 6.30in wide, the text width of the paper.
## Reads:  output/canonical/estimates/mainfit_{province_offset,national_by_stratum,sigma_u_profile}.csv
## Writes: figures/fig_coef_fourcountry.pdf
## Usage:  Rscript R/fig_mainfit_coefficients.R
## =====================================================================
suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
suppressMessages({library(dplyr); library(ggplot2)})

EST <- "output/canonical/estimates"

## decode the literal "<U+00E9>" escapes of the canonical CSVs (encoding= on read.csv does not)
fix_uesc <- function(x){
  m <- gregexpr("<U\\+([0-9A-Fa-f]{4})>", x)
  regmatches(x, m) <- lapply(regmatches(x, m), function(s)
    if(length(s)) intToUtf8(strtoi(sub("<U\\+([0-9A-Fa-f]{4})>", "\\1", s), 16L), multiple=TRUE) else s)
  x
}

po_csv <- read.csv(file.path(EST, "mainfit_province_offset.csv"),     stringsAsFactors=FALSE)
st_csv <- read.csv(file.path(EST, "mainfit_national_by_stratum.csv"), stringsAsFactors=FALSE)
su_csv <- read.csv(file.path(EST, "mainfit_sigma_u_profile.csv"),     stringsAsFactors=FALSE)

cts4 <- c("Mozambique","Zimbabwe","Malawi","Lesotho")

canon <- function(cty){
  pv <- po_csv %>% filter(country==cty) %>%
          transmute(country=cty, term=fix_uesc(province), type="province", OR, lo, hi)
  st <- st_csv %>% filter(country==cty, estimator=="EBLUP") %>%
          transmute(country=cty, term=trimws(gsub(" \\| ", " ", stratum)), type="age-sex", OR, lo, hi)
  bind_rows(pv, st)
}

typ_cols <- c(province="#E69F00", `age-sex`="#009E73")            # Okabe-Ito
## legend labels: the two blocks are marginal summaries of the same province x stratum cells
typ_labs <- c(province = "Province (marginal)", `age-sex` = "Age and sex (marginal)")

## age-sex rows in a fixed order: Female block above Male, age decreasing
agesex_rank <- function(term){
  sex <- ifelse(grepl("Female", term), 1, 2)
  age <- ifelse(grepl("50\\+", term), 1,
         ifelse(grepl("25-49", term), 2,
         ifelse(grepl("15-24", term), 3, 0)))                      # 0 = sex-only (Malawi)
  sex*10 + age
}

plot_forest <- function(d, xlim=NULL){
  d$sig  <- factor(ifelse(d$lo>1 | d$hi<1,"significant","not significant"), levels=c("significant","not significant"))
  d$typ2 <- factor(d$type, levels=c("age-sex","province"))         # province block on top
  d$ord  <- ifelse(d$type=="age-sex", -agesex_rank(d$term), d$OR)  # provinces by estimate
  d <- d %>% arrange(typ2, ord); d$term <- factor(d$term, levels=d$term)
  ggplot(d, aes(OR, term, colour=type)) +
    geom_vline(xintercept=1, linetype=2, colour="grey65") +
    geom_pointrange(aes(xmin=lo, xmax=hi, shape=sig, linetype=sig), size=.42) +
    scale_x_log10(limits=xlim, breaks=c(0.5,1,2,3,5)) +
    scale_colour_manual(values=typ_cols, labels=typ_labs, name=NULL) +
    scale_shape_manual(values=c(significant=16,`not significant`=1), name=NULL, drop=FALSE) +
    scale_linetype_manual(values=c(significant="solid",`not significant`="22"), name=NULL, drop=FALSE) +
    labs(x="Selection OR (log scale)", y=NULL, title=d$country[1]) +
    theme_minimal(base_size=11) +
    theme(plot.background=element_rect(fill="white",colour=NA), panel.background=element_rect(fill="white",colour=NA),
          panel.grid.minor=element_blank(), plot.title=element_text(face="bold",size=12),
          legend.position="none", axis.text.y=element_text(size=8.5), axis.title.x=element_text(size=10))
}

cdat  <- lapply(cts4, canon)
xr    <- range(unlist(lapply(cdat, function(d) c(d$lo, d$hi))), na.rm=TRUE) * c(0.92, 1.08)
plots <- lapply(cdat, function(d) plot_forest(d, xlim=xr))


sig_all <- su_csv %>% filter(estimable) %>%
  mutate(country=factor(country, levels=rev(cts4[cts4 %in% country])))
inset <- ggplot(sig_all, aes(sigma_u2, country)) +
  geom_vline(xintercept=0, linetype=2, colour="grey65") +
  geom_pointrange(aes(xmin=sigma_u2_lo, xmax=sigma_u2_hi), colour="grey25", linewidth=.5, size=.4) +
  scale_x_continuous(limits=c(0, NA), expand=expansion(mult=c(0.02, 0.08))) +
  scale_y_discrete(expand=expansion(add=0.5)) +
  labs(x=expression(sigma[u]^2), y=NULL, title="Within-province variance (95% profile CI)") +
  theme_minimal(base_size=10) +
  theme(plot.background=element_rect(fill="white",colour=NA), panel.background=element_rect(fill="white",colour=NA),
        panel.grid.minor=element_blank(), plot.title=element_text(size=9.5, face="bold"),
        axis.title.x=element_text(size=9), axis.text=element_text(size=9))

main   <- cowplot::plot_grid(plotlist=plots, ncol=2, labels="AUTO", label_size=13)
insrow <- cowplot::plot_grid(inset, NULL, ncol=2, rel_widths=c(1.1, 1), labels=c("E", ""), label_size=13)
## keys are described in the figure caption
g      <- cowplot::plot_grid(main, insrow, ncol=1, rel_heights=c(6.05, 1.25))
ggsave("figures/fig_coef_fourcountry.pdf", g, width=6.3, height=7.36)

cat("wrote fig_coef_fourcountry.pdf - ALL markers from the MAIN FIT\n")
cat("  province offsets (age-sex adjusted):\n")
for (cn in cts4) { d <- po_csv %>% filter(country==cn)
  cat(sprintf("    %-11s %2d provinces, OR %.2f-%.2f, %d exclude 1\n",
      cn, nrow(d), min(d$OR), max(d$OR), sum(d$lo>1 | d$hi<1))) }
cat("  sigma_u^2 (95% profile CI):\n")
for (i in seq_len(nrow(su_csv))) { r <- su_csv[i,]
  cat(sprintf("    %-11s %s\n", r$country,
      if (!r$estimable) "not estimable (1 area per province x stratum)"
      else sprintf("%.3f (%.3f, %.3f)", r$sigma_u2, r$sigma_u2_lo, r$sigma_u2_hi))) }
