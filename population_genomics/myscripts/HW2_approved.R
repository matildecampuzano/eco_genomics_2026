#Background info PCA and map####

library(tidyverse)

setwd("/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics/mydata")   # the folder with your plink output

#set stuff up to be read for pca
pca <- read_table("Lg_LD_0.2.eigenvec", col_names = FALSE)
eigenval <- scan("Lg_LD_0.2.eigenval")
pca <- pca[, -2]
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

head(pca)

# convert the eigenvalues to a percentage variance explained (pve)
pve <- data.frame(PC = seq_along(eigenval),
                  pve = eigenval / sum(eigenval) * 100)

#which pc to use?
pve_plot <- ggplot(pve, aes(PC, pve)) +
  geom_col() +
  labs(x = "Principal component", y = "Percentage variance explained") +
  theme_bw()
pve_plot

meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (4.2%)"
pc1_lab <- paste0("PC1 (", round(pve$pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)")

pca1_2_plot <- ggplot(pca, aes(PC1, PC2, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, colour = "Site (north to south)") +
  theme_bw()
pca1_2_plot

library(patchwork) ##combined plots
library(maps)

world <- map_data("world")
sites <- pca %>% distinct(Site, Latitude, Longitude)

map_plot <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey90", colour = "grey60") +
  geom_point(data = sites, aes(Longitude, Latitude, colour = Site), size = 3) +
  scale_colour_viridis_d(end = 0.9) +
  # zoom to our sites, with a 2-degree margin all round
  coord_quickmap(xlim = range(sites$Longitude) + c(-2, 2),
                 ylim = range(sites$Latitude) + c(-2, 2)) +
  labs(x = "Longitude", y = "Latitude", colour = "Site (north to south)") +
  theme_bw()
map_plot

# with patchwork, + puts plots side by side
# we hide the map's legend, as the PCA's legend already shows the sites
map_pca <- map_plot + theme(legend.position = "none") + pca1_2_plot
map_pca


#Now do BayPass w/ SSTmax and pHmin####
library(tidyverse)

setwd("~/Projects/eco_genomics_2026/population_genomics/mydata")

pops <- readLines("baypass/Lottia.baypass.pops.txt")
npop <- length(pops)

env <- read_csv("pop_level_biooracle.csv")
pop_order<-env$pop[order(env$lat)]

omega <- as.matrix(read.table("BayPassHomework/Lottia_bp_2_mat_omega.out"))
dimnames(omega) <- list(pops, pops)

omega_cor <- cov2cor(omega)
round(omega_cor, 2)

omega_tree <- hclust(as.dist(1 - omega_cor), method = "average")

### we are not using ggplot so it is saved in a different way:
plot(omega_tree, main = "Populations clustered by BayPass Omega", xlab = "", sub = "")
dev.off()

#Generate heatmapppp of dem
omega_df <- as.data.frame(as.table(omega_cor)) %>%
  rename(pop1 = Var1, pop2 = Var2, correlation = Freq) %>%
  mutate(pop1 = factor(pop1, levels = pop_order),
         pop2 = factor(pop2, levels = pop_order))

###if you get any errors for visuals, try dev.off()
omega_hw <- ggplot(omega_df, aes(x = pop1, y = pop2, fill = correlation)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "black", mid = "darkorchid4", high = "darkorange2",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL,
       title = "BayPass Omega as correlations (south to north)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
omega_hw
#look only at snps more differentiated than expected based on omega value
snps <- read_tsv("baypass/Lottia.baypass.snps.txt")

snps$MRK<-rownames(snps)
snps$MRK<-as.numeric(snps$MRK)

npop <- 19

xtx <- read.table("BayPassHomework/Lottia_bp_2_summary_pi_xtx.out", header = TRUE, check.names = FALSE) %>%
  left_join(snps, by = "MRK") %>%
  mutate(lp_high   = pchisq(XtXst, npop, lower.tail = FALSE, log.p = TRUE),   #
         lp_low    = pchisq(XtXst, npop, lower.tail = TRUE,  log.p = TRUE),   
         direction = if_else(lp_high < lp_low, "high", "low"),
         log10_p   = -(pmin(lp_high, lp_low) + log(2)) / log(10),             
         pval      = 10^(-log10_p),
         qval      = p.adjust(pval, method = "BH"))
#dont trust p values immedieately, first look at their distribution:
ggplot(xtx, aes(x = pval)) +
  geom_histogram(breaks = seq(0, 1, by = 0.025), fill = "violetred4", colour = "white") +
  labs(x = "p-value", y = "no. of SNPs") +
  theme_classic()



#correct genome position to make manhattan plot
##same helper tool as earlier with FST
add_genome_pos <- function(df) {
  scaffolds <- df %>% 
    group_by(CHROM) %>% ##group by CHROM
    summarise(length = max(POS)) %>% ## length = max position within a chromosome
    arrange(match(CHROM, unique(df$CHROM))) %>% ##make sure each is unique, arranged
    mutate(offset = cumsum(as.numeric(length)) - length, ###cumulatively number (helps with labels)
           shade  = factor(row_number() %% 2))    ##value for alternating colour 
  df %>%
    left_join(scaffolds, by = "CHROM") %>% ##join this data by given data with CHROM
    mutate(genome_pos = offset + POS) ###cumalative genome pos 
}


xtx <- xtx %>%
  add_genome_pos()     ##recall above if reloaded R

sum(xtx$outlier)                           

threshold <- quantile(xtx$qval, 0.001)   # top 0.1% (smallest q)

ggplot(xtx, aes(genome_pos / 1e6, -log10(qval))) +
  geom_point(aes(colour = shade), size = 0.4) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", colour = "red") +
  scale_colour_manual(values = c("grey20", "grey80"), guide = "none") +
  labs(x = "Genome pos (Mb)", y = "-log10(q)") +
  theme_classic()

xtx_outliers <- xtx %>%
  filter(qval < threshold) %>%
  arrange(desc(qval)) %>%
  select(MRK, CHROM, POS, M_XtX, XtXst, pval, qval)

head(xtx_outliers, 10)
write_tsv(xtx_outliers, "baypass/Lottia_bp_2_xtx_outliers.txt")

#genotype envrionment asscooation
covariates <- readLines("baypass/Lottia.baypass.covariates.txt")

bf <- read.table("BayPassHomework/Lottia_bp_2_summary_betai_reg.out", header = TRUE, check.names = FALSE) %>%
  rename(BF_dB = `BF(dB)`) %>%
  mutate(covariate = covariates[COVARIABLE]) %>%
  left_join(snps, by = "MRK") %>%
  add_genome_pos() %>%
  drop_na(
  )

bf %>%
  group_by(covariate) %>%
  summarise(n_decisive = sum(BF_dB > 20))

#maps stuff to environmental factors
ggplot(bf, aes(x = genome_pos / 1e6, y = BF_dB)) +
  geom_point(aes(colour = shade), size = 0.6) +
  geom_point(data = filter(bf, BF_dB > 20), colour = "red", size = 1) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  facet_wrap(~ covariate, ncol = 1) +
  scale_colour_manual(values = c("grey35", "grey70"), guide = "none") +
  labs(x = "Genome pos (Mb)",
       y = "Bayes factor (dB)") +
  theme_classic()

#double check stufff
geno <- as.matrix(read.table("baypass/Lottia.baypass.geno"))

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

library(dplyr)
library(ggplot2)
library(patchwork)

plot_top_snps <- function(cov, xlab, n = 5, method = "pearson") {
  dat <- get_top_freq(cov, n)
  
  # one correlation per SNP panel
  cors <- dat %>%
    group_by(snp) %>%
    summarise(
      r = cor(.data[[cov]], alt_freq, method = method, use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("r = %.2f", r))
  
  ggplot(dat, aes(x = .data[[cov]], y = alt_freq)) +
    geom_smooth(method = "lm", se = FALSE, color = "grey30",
                linewidth = 0.8, linetype = "dashed") +
    geom_point(aes(color = pop), size = 2) +
    geom_text(data = cors, aes(label = label),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
              size = 3.3, inherit.aes = FALSE) +
    scale_colour_viridis_d(end = 0.9) +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_wrap(~ snp, nrow = 1) +
    labs(x = xlab, y = "Allele frequency", color = "Population",
         title = paste("Top", n, "SNPs associated with", cov)) +
    theme_classic() +
    theme(
      panel.border  = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_rect(fill = "grey92", color = "grey40"),
      strip.text = element_text(face = "bold")
    )
}

p_sst <- plot_top_snps("thetao_max", "Max sea surface temperature (°C)")
p_ph  <- plot_top_snps("ph_min", "Minimum pH")

(p_sst / p_ph) + plot_layout(guides = "collect")