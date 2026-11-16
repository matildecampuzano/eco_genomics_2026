library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

get_top_freq <- function(cov, n = 5) {
  top <- bf %>%
    filter(covariate == cov) %>%
    slice_max(BF_dB, n = n, with_ties = FALSE) %>%
    mutate(snp = sprintf("%s:%s\nBF = %.1f dB", CHROM, POS, BF_dB))
  
  map_dfr(seq_len(nrow(top)), function(i) {
    counts <- as.numeric(geno[top$MRK[i], ])
    ref <- counts[seq(1, length(counts), by = 2)]
    alt <- counts[seq(2, length(counts), by = 2)]
    tibble(snp = top$snp[i], pop = pops, alt_freq = alt / (ref + alt))
  }) %>%
    mutate(snp = factor(snp, levels = top$snp)) %>%   # keeps rank order (strongest first)
    left_join(env, by = "pop") %>%
    mutate(pop = factor(pop, levels = unique(pop[order(-lat)])))
}

library(patchwork)

plot_top_snps <- function(cov, xlab, n = 5) {
  dat <- get_top_freq(cov, n)
  ggplot(dat, aes(x = .data[[cov]], y = alt_freq, color = pop)) +
    geom_point(size = 3) +
    scale_colour_viridis_d(end = 0.9) +
    scale_x_continuous(n.breaks = 3) +                 # avoids crowded x labels in narrow panels
    ylim(0, 1) +
    facet_wrap(~ snp, nrow = 1) +
    labs(x = xlab, y = "Allele frequency", color = "Population",
         title = paste("Top", n, "SNPs associated with", cov)) +
    theme_classic() +
    theme(
      panel.border  = element_rect(color = "grey40", fill = NA, linewidth = 0.8),  # box around each panel
      panel.spacing = unit(0.5, "lines"),                                          # gap between panels
      strip.background = element_rect(fill = "grey92", color = "grey40"),
      strip.text = element_text(face = "bold")
    )
}

p_sst <- plot_top_snps("thetao_max", "Max sea surface temperature (°C)")
p_ph  <- plot_top_snps("ph_min", "Minimum pH")

(p_sst / p_ph) + plot_layout(guides = "collect")
