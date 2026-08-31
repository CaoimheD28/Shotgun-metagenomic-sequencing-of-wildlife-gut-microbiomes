setwd("C:/Users/cedoy/Thesis/Chapter4/bracken_reads_all")
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(tools)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(vegan)
library(zCompositions)
library(compositions)
library(stringr)



# Rarefaction curves using MICROBIAL-ONLY read counts (Chordata excluded)
# This redoes the original rarefaction analysis, but now using the read
# counts that actually matter for phylum-level microbiome diversity -
# i.e. with host/dietary Chordata reads removed first.


## Figure S4.1 

# --- 1. Load the full phylum count matrix (blank-row-safe read) ---
raw_matrix <- read.csv("phylum_count_matrix.csv", check.names = FALSE)
raw_matrix <- raw_matrix[raw_matrix[[1]] != "" & !is.na(raw_matrix[[1]]), ]
count_matrix_full <- as.matrix(raw_matrix[, -1])
rownames(count_matrix_full) <- raw_matrix[[1]]

cat("Full matrix (all 57 samples, including Chordata):",
    nrow(count_matrix_full), "samples x", ncol(count_matrix_full), "phyla\n")

# --- 2. Remove Chordata BEFORE assessing depth ---
# Starting fresh from all 57 samples here - the original 7 exclusions were
# based on total (Chordata-inclusive) depth, which we now know was
# sometimes a poor proxy for microbial depth. This re-evaluates everyone.
count_matrix_microbial <- count_matrix_full[, colnames(count_matrix_full) != "Chordata"]

# Remove any samples that now have zero microbial reads at all (can't be
# rarefied) - shouldn't happen based on what we've seen, but a safe check
zero_read_samples <- rownames(count_matrix_microbial)[rowSums(count_matrix_microbial) == 0]
if (length(zero_read_samples) > 0) {
  cat("\nWarning: these samples have ZERO microbial reads and will be dropped:\n")
  print(zero_read_samples)
  count_matrix_microbial <- count_matrix_microbial[rowSums(count_matrix_microbial) > 0, ]
}

cat("\nMicrobial-only matrix:", nrow(count_matrix_microbial), "samples x",
    ncol(count_matrix_microbial), "phyla\n")

# --- 3. Plot rarefaction curves for ALL samples, colour-coded by species ---
sample_names <- rownames(count_matrix_microbial)
host_species <- ifelse(grepl("fox", sample_names, ignore.case = TRUE), "fox",
                       ifelse(grepl("badger", sample_names, ignore.case = TRUE), "badger", "unknown"))
curve_colors <- ifelse(host_species == "fox", "#1D9E75",
                       ifelse(host_species == "badger", "#7F77DD", "grey"))

png("rarefaction_curves_microbial_only.png", width = 1200, height = 900, res = 120)
rarecurve(count_matrix_microbial, step = 500, col = curve_colors, label = FALSE,
          xlab = "Number of microbial (non-Chordata) reads",
          ylab = "Number of phyla detected",
          main = "Rarefaction curves: microbial-only read depth")
legend("bottomright", legend = c("fox", "badger"), col = c("#1D9E75", "#7F77DD"), lty = 1)
dev.off()
cat("\nSaved rarefaction_curves_microbial_only.png\n")

# --- 4. Zoomed plot: the 7 excluded samples specifically, labelled ---
excluded_samples <- c(
  "0121002_bracken_phylum_badger", "0227016_bracken_phylum_badger",
  "0155007_bracken_phylum_badger", "0227007_bracken_phylum_badger",
  "0178008_bracken_phylum_badger", "0149008_bracken_phylum_badger",
  "0155005_bracken_phylum_badger"
)

excluded_matrix <- count_matrix_microbial[excluded_samples, , drop = FALSE]
excluded_host <- ifelse(grepl("fox", excluded_samples, ignore.case = TRUE), "fox",
                        ifelse(grepl("badger", excluded_samples, ignore.case = TRUE), "badger", "unknown"))
excluded_colors <- ifelse(excluded_host == "fox", "#1D9E75",
                          ifelse(excluded_host == "badger", "#7F77DD", "grey"))

if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)

# Give each of the 8 samples its own distinct colour (rather than just
# fox/badger colour-coding), so a side legend can identify each curve
# unambiguously without needing in-plot labels at all.
sample_palette <- brewer.pal(7, "Dark2")
names(sample_palette) <- excluded_samples
short_labels <- sub("_bracken_phylum_badger|_bracken_phylum_fox", "", excluded_samples)

png("rarefaction_curves_excluded.png", width = 1600, height = 900, res = 130)
layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))  # plot area + legend area

par(mar = c(5, 5, 4, 1))
rarecurve(excluded_matrix, step = 20, col = sample_palette[rownames(excluded_matrix)],
          label = FALSE, lwd = 1.8,
          xlab = "Number of microbial (non-Chordata) reads",
          ylab = "Number of phyla detected",
          main = "Rarefaction curves: 7 excluded low-depth samples")

par(mar = c(5, 0, 4, 1))
plot.new()
legend("center", legend = short_labels, col = sample_palette[excluded_samples],
       lwd = 2.5, cex = 0.85, bty = "n", title = "Sample")

dev.off()
cat("Saved rarefaction_curves_excluded.png\n")

# --- 5. Print the 8 excluded samples' depth for quick reference ---
cat("\nExcluded samples and their microbial read depth:\n")
print(sort(rowSums(excluded_matrix)))

# --- 6. Save the microbial-only matrix for downstream re-use ---
write.csv(count_matrix_microbial, "phylum_count_matrix_microbial_only.csv")



# consistent depth assessment: Chao1 completeness for ALL 57 samples
# (microbial-only reads, Chordata already excluded)
#
# This replaces the piecemeal approach of the last several steps with one
# single, consistently-applied check across the whole dataset, so the
# final exclusion list is decided the same way for every sample.


# --- 1. Load microbial-only matrix ---
count_matrix_microbial <- as.matrix(read.csv("phylum_count_matrix_microbial_only.csv",
                                             row.names = 1, check.names = FALSE))

cat("Assessing", nrow(count_matrix_microbial), "samples\n\n")

# --- 2. Calculate Chao1 completeness for every sample ---
all_results <- data.frame(
  sample = rownames(count_matrix_microbial),
  total_microbial_reads = rowSums(count_matrix_microbial),
  S_obs = NA,
  S_chao1 = NA,
  completeness = NA
)

for (i in seq_len(nrow(count_matrix_microbial))) {
  s <- rownames(count_matrix_microbial)[i]
  est <- estimateR(round(count_matrix_microbial[s, ]))
  all_results$S_obs[i] <- est["S.obs"]
  all_results$S_chao1[i] <- round(est["S.chao1"], 1)
  all_results$completeness[i] <- round(est["S.obs"] / est["S.chao1"], 3)
}

all_results$host_species <- ifelse(grepl("fox", all_results$sample, ignore.case = TRUE), "fox",
                                   ifelse(grepl("badger", all_results$sample, ignore.case = TRUE), "badger", "unknown"))

# --- 3. IMPORTANT CAVEAT ---
# Chao1 completeness is unreliable at very low read counts (too few
# observations for the estimator to detect rare/singleton taxa properly -
# it can misleadingly show HIGH "completeness" simply because there's too
# little data to reveal how incomplete the sample really is). So exclusion
# decisions for very low-depth samples should rely primarily on RAW READ
# COUNT, not completeness score. Completeness is most informative and
# trustworthy as a discriminator in the low-to-moderate depth range (roughly
# a few thousand reads and up), where there's enough data for the estimator
# to work properly, but the sample still might be under-sequenced.

all_results <- all_results[order(all_results$total_microbial_reads), ]

write.csv(all_results, "final_depth_assessment.csv", row.names = FALSE)

cat("=== All 57 samples, sorted by microbial read depth (lowest first) ===\n")
print(all_results, row.names = FALSE)

cat("\n\n=== Same data, sorted by completeness (lowest first) ===\n")
cat("(remember: completeness is unreliable below ~1,000-2,000 reads - check\n")
cat(" total_microbial_reads alongside it, don't use completeness alone)\n\n")
print(all_results[order(all_results$completeness), ], row.names = FALSE)

## -----------------------------------------------------------------------------

## Figure 4.1

# Alpha diversity analysis: fox vs badger
# Uses the phylum count matrix built by rarefaction_curves.R
# (count_matrix: rows = samples, columns = phyla, values = read counts)

library(vegan)

# --- 1. Load the microbial-only matrix (Chordata already removed) ---
count_matrix_clean <- as.matrix(read.csv("phylum_count_matrix_microbial_only.csv",
                                         row.names = 1, check.names = FALSE))

# --- 2. Remove the 7 excluded low-depth samples (final list) ---
excluded_samples <- c(
  "0121002_bracken_phylum_badger", "0227016_bracken_phylum_badger",
  "0155007_bracken_phylum_badger", "0227007_bracken_phylum_badger",
  "0178008_bracken_phylum_badger", "0149008_bracken_phylum_badger",
  "0155005_bracken_phylum_badger"
)

count_matrix_clean <- count_matrix_clean[!(rownames(count_matrix_clean) %in% excluded_samples), ]
cat("Samples remaining after exclusion:", nrow(count_matrix_clean), "\n")

# --- 3. Calculate alpha diversity metrics per sample ---
richness   <- specnumber(count_matrix_clean)               # number of phyla detected
shannon    <- diversity(count_matrix_clean, index = "shannon")
simpson    <- diversity(count_matrix_clean, index = "simpson")

alpha_div <- data.frame(
  sample    = rownames(count_matrix_clean),
  richness  = richness,
  shannon   = shannon,
  simpson   = simpson
)

# --- 4. Add host species labels ---
alpha_div$host_species <- ifelse(grepl("fox", alpha_div$sample, ignore.case = TRUE), "fox",
                                 ifelse(grepl("badger", alpha_div$sample, ignore.case = TRUE), "badger", "unknown"))

write.csv(alpha_div, "alpha_diversity_results.csv", row.names = FALSE)

# --- 6. Compare fox vs badger for each metric ---
cat("\n--- Richness ---\n")
print(wilcox.test(richness ~ host_species, data = alpha_div, exact = FALSE))

cat("\n--- Shannon diversity ---\n")
print(wilcox.test(shannon ~ host_species, data = alpha_div, exact = FALSE))

cat("\n--- Simpson diversity ---\n")
print(wilcox.test(simpson ~ host_species, data = alpha_div, exact = FALSE))

# Summary stats per group
cat("\n--- Summary stats ---\n")
print(tapply(alpha_div$richness, alpha_div$host_species, summary))
print(tapply(alpha_div$shannon, alpha_div$host_species, summary))
print(tapply(alpha_div$simpson, alpha_div$host_species, summary))

# --- 7. Boxplots ---
png("alpha_diversity_boxplots.png", width = 1600, height = 600, res = 120)
par(mfrow = c(1, 3), cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.6, mar = c(5, 5, 4, 2))

boxplot(richness ~ host_species, data = alpha_div, main = "Richness",
        ylab = "Number of phyla", col = c("#4C6E5D", "#C97B3D"))
boxplot(shannon ~ host_species, data = alpha_div, main = "Shannon diversity",
        ylab = "Shannon index", col = c("#4C6E5D", "#C97B3D"))
boxplot(simpson ~ host_species, data = alpha_div, main = "Simpson diversity",
        ylab = "Simpson index", col = c("#4C6E5D", "#C97B3D"))

dev.off()
cat("\nSaved plot to alpha_diversity_boxplots.png\n")

## -----------------------------------------------------------------------------

## Figure 4.2


# Beta diversity analysis: fox vs badger
# CLR transformation + PERMANOVA + PCA visualisation
# Uses the same cleaned phylum count matrix as alpha_diversity.R


# --- STEP 1. Load the microbial-only matrix (Chordata already removed) and exclude low-depth samples ---
count_matrix_clean <- as.matrix(read.csv("phylum_count_matrix_microbial_only.csv",
                                         row.names = 1, check.names = FALSE))

excluded_samples <- c(
  "0121002_bracken_phylum_badger", "0227016_bracken_phylum_badger",
  "0155007_bracken_phylum_badger", "0227007_bracken_phylum_badger",
  "0178008_bracken_phylum_badger", "0149008_bracken_phylum_badger",
  "0155005_bracken_phylum_badger"
)
count_matrix_clean <- count_matrix_clean[!(rownames(count_matrix_clean) %in% excluded_samples), ]

# Remove any phyla with zero reads across ALL remaining samples (uninformative)
count_matrix_clean <- count_matrix_clean[, colSums(count_matrix_clean) > 0]

cat("Final matrix:", nrow(count_matrix_clean), "samples x", ncol(count_matrix_clean), "phyla\n")

# STEP 2. Handle zeros properly before CLR
# cmultRepl uses a Bayesian-multiplicative replacement method - standard
# approach for compositional microbiome data, better than adding an
# arbitrary small constant everywhere.
count_matrix_nozero <- cmultRepl(count_matrix_clean, method = "CZM", output = "p-counts")

# STEP 3. CLR transform ---
clr_data <- clr(count_matrix_nozero)
clr_matrix <- as.matrix(clr_data)

# STEP 4. Host species labels ---
host_species <- ifelse(grepl("fox", rownames(clr_matrix), ignore.case = TRUE), "fox",
                       ifelse(grepl("badger", rownames(clr_matrix), ignore.case = TRUE), "badger", "unknown"))

# STEP 5. PERMANOVA using Aitchison distance (= Euclidean distance on CLR data)
aitchison_dist <- dist(clr_matrix, method = "euclidean")

permanova_result <- adonis2(aitchison_dist ~ host_species, permutations = 999)
print(permanova_result)

# STEP 6. Check dispersion too
# PERMANOVA can be significant either because groups differ in true
# composition, OR because one group is just more variable than the other.
# This check tells you which is going on.
dispersion <- betadisper(aitchison_dist, host_species)
print(anova(dispersion))

# STEP 7. PCA plot for visualisation
pca_result <- prcomp(clr_matrix)
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$host_species <- host_species

var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)


# Colours matched to the alpha diversity boxplots (navy/dark green for
# badger, orange for fox) for visual consistency across figures
species_colors <- c("badger" = "#3B6E5B", "fox" = "#D97D26")

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = host_species, fill = host_species)) +
  stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, linewidth = 0.8) +
  geom_point(size = 3.5, shape = 21, color = "black", stroke = 0.4,
             aes(fill = host_species)) +
  scale_color_manual(values = species_colors, name = "Host species") +
  scale_fill_manual(values = species_colors, name = "Host species") +
  labs(
    title = "PCA of CLR-transformed phylum composition",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "right",
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

ggsave("beta_diversity_pca.png", p_pca, width = 8, height = 5.5, dpi = 300)
ggsave("beta_diversity_pca.pdf", p_pca, width = 8, height = 5.5)  # vector version for print

print(p_pca)

# --- 9. Save results for reference ---
write.csv(pca_scores, "beta_diversity_pca_scores.csv")

pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$host_species <- host_species
pca_scores$sample <- rownames(pca_scores)

# Top 3 by PC1 (the two on the right)
pca_scores[order(-pca_scores$PC1), ][1:3, ]

# Bottom 3 by PC2 (the one at the bottom)
pca_scores[order(pca_scores$PC2), ][1:3, ]



# Identify which phyla drive each outlier sample's PCA position, and check
# the raw read counts behind each driver phylum to assess whether it's a
# trustworthy detection or likely noise (lesson learned from 0155005/etc.)
# Run this AFTER beta_diversity.R in the same session
# (reuses clr_matrix and count_matrix_clean)

outlier_samples <- c(
  "0215005_bracken_phylum_badger", "0125005_bracken_phylum_badger",
  "0086002_bracken_phylum_badger", "0227010_bracken_phylum_badger",
  "0218003_bracken_phylum_fox", "0209006_bracken_phylum_badger"
)

clr_means <- colMeans(clr_matrix)

for (s in outlier_samples) {
  cat("\n\n============================================\n")
  cat("Sample:", s, "\n")
  cat("Total microbial reads:", sum(count_matrix_clean[s, ]), "\n")
  cat("============================================\n")
  
  # Which phyla are most unusual (CLR-based) for this sample
  sample_diff <- clr_matrix[s, ] - clr_means
  sample_diff_sorted <- sort(sample_diff, decreasing = TRUE)
  
  top_high <- head(sample_diff_sorted, 5)
  top_low <- tail(sample_diff_sorted, 5)
  
  cat("\nPhyla with UNUSUALLY HIGH abundance (vs. average sample):\n")
  for (p in names(top_high)) {
    raw_count <- count_matrix_clean[s, p]
    pct <- round(raw_count / sum(count_matrix_clean[s, ]) * 100, 3)
    cat(sprintf("  %-25s CLR diff = %6.2f | raw reads = %6d | %% of sample = %s%%\n",
                p, top_high[p], raw_count, pct))
  }
  
  cat("\nPhyla with UNUSUALLY LOW abundance (vs. average sample):\n")
  for (p in names(top_low)) {
    raw_count <- count_matrix_clean[s, p]
    pct <- round(raw_count / sum(count_matrix_clean[s, ]) * 100, 3)
    cat(sprintf("  %-25s CLR diff = %6.2f | raw reads = %6d | %% of sample = %s%%\n",
                p, top_low[p], raw_count, pct))
  }
}

cat("\n\nReminder: treat any driver phylum with only 1-5 raw reads with caution -\n")
cat("likely noise. Drivers with dozens+ reads and/or a meaningful % of the\n")
cat("sample are more likely to reflect genuine biological signal.\n")


## -----------------------------------------------------------------------------

## Figure 4.3

# Stacked bar plot: relative abundance of phyla per sample
# Publication-style version using ggplot2

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("tidyr")) install.packages("tidyr")
if (!require("dplyr")) install.packages("dplyr")
library(ggplot2)
library(tidyr)
library(dplyr)

# --- 1. Load the microbial-only matrix (Chordata already removed) ---
count_matrix_clean <- as.matrix(read.csv("phylum_count_matrix_microbial_only.csv",
                                         row.names = 1, check.names = FALSE))

excluded_samples <- c(
  "0121002_bracken_phylum_badger", "0227016_bracken_phylum_badger",
  "0155007_bracken_phylum_badger", "0227007_bracken_phylum_badger",
  "0178008_bracken_phylum_badger", "0149008_bracken_phylum_badger",
  "0155005_bracken_phylum_badger"
)
count_matrix_clean <- count_matrix_clean[!(rownames(count_matrix_clean) %in% excluded_samples), ]

# --- 3. Convert to relative abundance (%) per sample ---
rel_abund <- sweep(count_matrix_clean, 1, rowSums(count_matrix_clean), "/") * 100

# --- 4. Collapse rare phyla into "Other" ---
# A phylum is kept individually only if it exceeds this threshold in AT LEAST
# one sample; otherwise it's grouped into "Other". Lower = more phyla shown
# individually (more colours, more detail); higher = simpler plot, more
# grouped into Other.
threshold <- 0.2
keep_phyla <- colnames(rel_abund)[apply(rel_abund, 2, max) >= threshold]

rel_abund_grouped <- as.data.frame(rel_abund) %>%
  mutate(sample = rownames(rel_abund)) %>%
  pivot_longer(-sample, names_to = "phylum", values_to = "abundance") %>%
  mutate(phylum = ifelse(phylum %in% keep_phyla, phylum, "Other")) %>%
  group_by(sample, phylum) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# --- 5. Add host species, for ordering samples sensibly on the x-axis ---
rel_abund_grouped$host_species <- ifelse(
  grepl("fox", rel_abund_grouped$sample, ignore.case = TRUE), "Fox",
  ifelse(grepl("badger", rel_abund_grouped$sample, ignore.case = TRUE), "Badger", "Unknown")
)

# Order samples: badgers first then foxes, alphabetically within each
sample_order <- rel_abund_grouped %>%
  distinct(sample, host_species) %>%
  arrange(host_species, sample) %>%
  pull(sample)

rel_abund_grouped$sample <- factor(rel_abund_grouped$sample, levels = sample_order)

# --- 5b. Simplify sample labels to just the numeric ID ---
# e.g. "0097005_bracken_phylum_badger" -> "97005"
# This creates a separate label column so the underlying sample factor
# (used for ordering) is untouched, but the plot displays the short version.
simplify_id <- function(x) {
  as.character(as.numeric(sub("^0*([0-9]+)_.*", "\\1", x)))
}
label_lookup <- setNames(simplify_id(sample_order), sample_order)

# Order phyla by overall mean abundance (largest first), "Other" always last
phylum_order <- rel_abund_grouped %>%
  filter(phylum != "Other") %>%
  group_by(phylum) %>%
  summarise(mean_abund = mean(abundance)) %>%
  arrange(desc(mean_abund)) %>%
  pull(phylum)
phylum_order <- c(phylum_order, "Other")

rel_abund_grouped$phylum <- factor(rel_abund_grouped$phylum, levels = phylum_order)

# --- 6. Colour palette: hand-picked for maximum distinction, especially ---
#         among the most abundant phyla. "Other" always grey.
# This is a named lookup - phyla are matched by name if present, so the
# palette stays correct even if which phyla appear changes slightly (e.g.
# after re-running with a different threshold).
manual_colors <- c(
  "Pseudomonadota"          = "#64B5CD",  # dominant phylum #1
  "Bacillota"               = "#E68A58",  # dominant phylum #2
  "Bacteroidota"            = "#6FAE50",  # medium green
  "Fusobacteriota"          = "#2E7D32",  # darker green (distinct from Bacteroidota)
  "Actinomycetota"          = "#F2B6C1",  # soft pink
  "Uroviricota"             = "#C0143C",  # strong red (distinct from Actinomycetota)
  "Nucleocytoviricota"      = "#FFD34D",  # yellow-gold
  "Campylobacterota"        = "#7A4FA3",  # purple
  "Planctomycetota"         = "#4FB8D8",  # sky blue
  "Thermodesulfobacteriota" = "#8C564B",  # brown
  "Mycoplasmatota"          = "#E7CB94",  # tan
  "Euryarchaeota"           = "#5B5B5B",  # dark grey
  "Spirochaetota"           = "#B565A7",  # mauve
  "Myxococcota"             = "#F49AD8",  # bright pink
  "Verrucomicrobiota"       = "#9BC4E2",  # pale blue
  "Acidobacteriota"         = "#3CB371"   # sea green, in case it reappears
)

n_phyla <- length(phylum_order) - 1  # excluding "Other"
named_phyla <- phylum_order[phylum_order != "Other"]

# For any phylum not explicitly listed above (e.g. if your data has
# additional rare-but-above-threshold phyla), fall back to a distinct
# extra palette so nothing is left uncoloured.
missing_phyla <- setdiff(named_phyla, names(manual_colors))
if (length(missing_phyla) > 0) {
  if (!require("RColorBrewer")) install.packages("RColorBrewer")
  library(RColorBrewer)
  fallback_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Dark2"))(length(missing_phyla)),
    missing_phyla
  )
  manual_colors <- c(manual_colors, fallback_colors)
}

palette <- c(manual_colors[named_phyla], "grey80")
names(palette) <- c(named_phyla, "Other")

# --- 7. Build the plot ---
p <- ggplot(rel_abund_grouped, aes(x = sample, y = abundance, fill = phylum)) +
  geom_bar(stat = "identity", width = 0.85) +
  scale_x_discrete(labels = label_lookup) +
  scale_fill_manual(values = palette, name = "Phylum") +
  facet_grid(~ host_species, scales = "free_x", space = "free_x") +
  labs(
    title = "Relative abundance of phyla by sample",
    x = "Sample",
    y = "Relative abundance (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14)
  )

# --- 8. Save as high-resolution image, suitable for a thesis/paper ---
ggsave("phylum_relative_abundance.png", p, width = 12, height = 6, dpi = 300)
ggsave("phylum_relative_abundance.pdf", p, width = 12, height = 6)  # vector version, good for print

print(p)
cat("\nSaved phylum_relative_abundance.png and .pdf\n")


################################################################################


setwd("C:/Users/cedoy/Thesis/Chapter4/Reads/Kraken_Reads_Resis")

# Screen the READ-based Kraken2 report files (not the assembly-based ones)
# for any member of the Mycobacterium tuberculosis complex (MTBC).
#
# Same logic as the earlier assembly-based MTBC screen, just pointed at
# your original read-based Kraken2 report files.

library(stringr)
library(dplyr)

# --- UPDATE this to the folder containing your READ-based Kraken2 report
#     files (i.e. from your original per-sample Kraken2 runs on reads,
#     NOT the assembly-based ones in Kraken_files_assemblies) ---
report_dir <- "C:/Users/cedoy/Thesis/Chapter4/Reads/Kraken_Reads_Resis"

mtbc_pattern <- "Mycobacterium tuberculosis|Mycobacterium bovis|Mycobacterium africanum|Mycobacterium caprae|Mycobacterium microti|Mycobacterium pinnipedii|Mycobacterium orygis|Mycobacterium canettii"

report_files <- list.files(report_dir, pattern = "_kraken2_reads_report.txt$", full.names = TRUE)
cat("Found", length(report_files), "read-based report files\n")
cat("(should be 57 if pointed at the right folder)\n\n")

all_hits <- list()

for (r in report_files) {
  sample <- sub("_kraken2_reads_report.txt$", "", basename(r))
  
  report <- read.delim(r, header = FALSE, sep = "\t",
                       col.names = c("pct", "reads_clade", "reads_direct",
                                     "rank_code", "taxid", "name_indented"),
                       stringsAsFactors = FALSE)
  report$name <- str_trim(report$name_indented)
  
  hits <- report %>%
    dplyr::filter(str_detect(name, mtbc_pattern), reads_clade > 0)
  
  if (nrow(hits) > 0) {
    hits$sample <- sample
    all_hits[[sample]] <- hits %>%
      dplyr::select(sample, name, taxid, rank_code, reads_clade, reads_direct, pct)
  }
}

if (length(all_hits) == 0) {
  cat("No MTBC hits detected in ANY of the", length(report_files), "samples (read-based).\n")
} else {
  mtbc_results <- do.call(rbind, all_hits)
  rownames(mtbc_results) <- NULL
  
  cat("MTBC hits found in", length(unique(mtbc_results$sample)), "of",
      length(report_files), "samples (read-based):\n\n")
  print(mtbc_results)
  
  write.csv(mtbc_results, "MTBC_screening_results_READS.csv", row.names = FALSE)
  cat("\nSaved to MTBC_screening_results_READS.csv\n")
}







################################################################################
################################################################################
################################################################################

## In order to compare the rural and urban foxes, we need the bracken phyla files
## merged for the rural foxes


library(dplyr)
library(purrr)

# --- UPDATE this to the folder containing your 11 rural fox Bracken files ---
bracken_dir <- "C:/Users/cedoy/Thesis/Chapter4/bracken_reads_all"

cat("Current working directory:", getwd(), "\n")

files <- list.files(bracken_dir, pattern = "_bracken_phylum_fox.txt$", full.names = TRUE)
cat("Found", length(files), "files (should be 11)\n")
print(files)

if (length(files) == 0) {
  stop("No files found - check getwd() above matches the folder containing your Bracken files")
}

# --- Read each file and extract sample name + new_est_reads ---
# Standard Bracken columns: name, taxonomy_id, taxonomy_lvl,
# kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads
read_one_bracken <- function(filepath) {
  sample_name <- sub("_bracken_phylum_fox.txt$", "", basename(filepath))
  
  df <- read.delim(filepath, header = TRUE, stringsAsFactors = FALSE)
  
  df %>%
    select(name, taxonomy_id, taxonomy_lvl, new_est_reads) %>%
    rename(!!sample_name := new_est_reads)
}

all_samples <- map(files, read_one_bracken)

cat("\nNumber of samples successfully read:", length(all_samples), "\n")
if (length(all_samples) > 0) {
  cat("First sample preview:\n")
  print(head(all_samples[[1]]))
}

if (length(all_samples) == 0) {
  stop("all_samples is empty - something failed while reading the files, even though they were found")
}

# --- Merge all samples together on name/taxonomy_id/taxonomy_lvl ---
# Using full_join so a phylum present in some samples but not others still
# gets a row (filled with NA, converted to 0 below - i.e. genuinely absent)
merged <- reduce(all_samples, full_join, by = c("name", "taxonomy_id", "taxonomy_lvl"))

# Replace NA with 0 (a phylum not listed in that sample's Bracken report
# genuinely had zero reads, not missing data)
sample_cols <- setdiff(colnames(merged), c("name", "taxonomy_id", "taxonomy_lvl"))
merged[sample_cols] <- lapply(merged[sample_cols], function(x) ifelse(is.na(x), 0, round(x)))

cat("\nMerged table:", nrow(merged), "phyla x", length(sample_cols), "samples\n")

write.csv(merged, "bracken_phylum_merged_rural.csv", row.names = FALSE)
cat("Saved to bracken_phylum_merged_rural.csv\n")

################################################################################


## Figure 4.4 - Mapping phyla

# --- 1. Load coordinates ---
coords <- read.csv("sequenced_samples_x_y_low_removed.csv", stringsAsFactors = FALSE)
coords$xcoord <- suppressWarnings(as.numeric(coords$xcoord))
coords$ycoord <- suppressWarnings(as.numeric(coords$ycoord))
coords <- coords[!is.na(coords$xcoord) & !is.na(coords$ycoord), ]

# Standardise Sample ID to plain integer (strips any leading zeros)
coords$Sample_id <- as.integer(coords$Sample)
cat("Coordinates loaded:", nrow(coords), "samples with valid location\n")

# --- Add a small jitter to ALL sample coordinates, so overlapping or
#     closely-spaced points are easier to distinguish on the map. Jitter
#     is small (+/- up to 500m) relative to national map scale, so it
#     does not meaningfully misrepresent sample location. ---
set.seed(42)  # reproducible jitter
n <- nrow(coords)
coords$xcoord <- coords$xcoord + runif(n, -500, 500)
coords$ycoord <- coords$ycoord + runif(n, -500, 500)

# --- 2. Load phylum Bracken file ---
# Structure: one row per sample, "Sample" column + one column per phylum
bracken_phylum <- read.csv("phylum_count_matrix_microbial_only_low_removed.csv", check.names = FALSE, stringsAsFactors = FALSE)
cat("\nBracken phylum file columns (first 10):\n")
print(head(colnames(bracken_phylum), 10))
cat("Total columns:", ncol(bracken_phylum), "\n")
cat("\nRaw Sample column values (first 5):\n")
print(head(bracken_phylum$Sample, 5))
cat("Class of Sample column:", class(bracken_phylum$Sample), "\n")

# Extract just the digits from the Sample value, in case it contains extra
# text (e.g. "0086002_bracken_phylum_badger" rather than a plain number)
bracken_phylum$Sample_id <- as.integer(gsub("\\D", "", as.character(bracken_phylum$Sample)))
cat("Sample_id after extraction (first 5):\n")
print(head(bracken_phylum$Sample_id, 5))
cat("Number of NA Sample_id values:", sum(is.na(bracken_phylum$Sample_id)), "\n")

# --- 3. Reshape to long format and calculate relative abundance ---
# Exclude Chordata (host DNA), consistent with earlier analyses
phylum_cols <- setdiff(colnames(bracken_phylum), c("Sample", "Sample_id", "Chordata"))
long_data <- bracken_phylum %>%
  dplyr::select(Sample_id, all_of(phylum_cols)) %>%
  tidyr::pivot_longer(all_of(phylum_cols), names_to = "name", values_to = "reads")

long_data <- long_data %>%
  group_by(Sample_id) %>%
  mutate(rel_abund = reads / sum(reads)) %>%
  ungroup()

# --- 4. Find each sample's DOMINANT phylum ---
dominant_phylum <- long_data %>%
  group_by(Sample_id) %>%
  slice_max(rel_abund, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Sample_id, dominant_phylum = name, dominant_pct = rel_abund)

# --- 5. Merge with coordinates ---
merged <- coords %>%
  inner_join(dominant_phylum, by = "Sample_id")

cat("\nSuccessfully merged:", nrow(merged), "of", nrow(coords), "samples with coordinates\n")
if (nrow(merged) < nrow(coords)) {
  cat("Samples with coordinates but NO matching phylum data:\n")
  print(coords$Sample[!coords$Sample_id %in% merged$Sample_id])
}

write.csv(merged, "phylum_spatial_merged.csv", row.names = FALSE)
cat("Saved phylum_spatial_merged.csv\n")

# --- 6. Map: points coloured by dominant phylum ---
ireland_full <- ne_countries(country = "Ireland", scale = "medium", returnclass = "sf")
merged_sf <- st_as_sf(merged, coords = c("xcoord", "ycoord"), crs = 29903)
merged_sf <- st_transform(merged_sf, crs = 4326)

# Group rare dominant phyla into "Other" for a cleaner legend
top_dominant <- merged %>% count(dominant_phylum, sort = TRUE) %>% slice_head(n = 8) %>% pull(dominant_phylum)
merged_sf$dominant_grouped <- ifelse(merged_sf$dominant_phylum %in% top_dominant,
                                     merged_sf$dominant_phylum, "Other")

p <- ggplot() +
  geom_sf(data = ireland_full, fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_sf(data = merged_sf, aes(color = dominant_grouped), size = 4, alpha = 0.85) +
  labs(title = "Dominant phylum by sample location", color = "Dominant\nphylum") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 14)
  )
print(p)

ggsave("dominant_phylum_map.png", p, width = 9, height = 10, dpi = 300)
cat("\nSaved dominant_phylum_map.png\n")


# --- 1. Load and clean coordinates ---
coords <- read.csv("sequenced_samples_x_y_low_removed.csv", stringsAsFactors = FALSE)
coords$xcoord <- suppressWarnings(as.numeric(coords$xcoord))
coords$ycoord <- suppressWarnings(as.numeric(coords$ycoord))
coords <- coords[!is.na(coords$xcoord) & !is.na(coords$ycoord), ]
coords$Sample_id <- as.integer(coords$Sample)

# --- 2. Load phylum Bracken file and build the full sample x phylum matrix ---
bracken_phylum <- read.csv("phylum_count_matrix_microbial_only_low_removed.csv", check.names = FALSE, stringsAsFactors = FALSE)
bracken_phylum$Sample_id <- as.integer(gsub("\\D", "", as.character(bracken_phylum$Sample)))

phylum_cols <- setdiff(colnames(bracken_phylum), c("Sample", "Sample_id", "Chordata"))

# --- 3. Restrict to samples present in BOTH datasets, and align order exactly ---
common_ids <- intersect(coords$Sample_id, bracken_phylum$Sample_id)
cat("Samples with both coordinates and phylum data:", length(common_ids), "\n")

coords_matched <- coords[match(common_ids, coords$Sample_id), ]
bracken_matched <- bracken_phylum[match(common_ids, bracken_phylum$Sample_id), ]

# Double-check alignment before proceeding - this MUST be TRUE
cat("Sample order aligned correctly:", all(coords_matched$Sample_id == bracken_matched$Sample_id), "\n")

# --- 4. Build the phylum count matrix, remove all-zero phyla ---
phylum_matrix <- as.matrix(bracken_matched[, phylum_cols])
rownames(phylum_matrix) <- bracken_matched$Sample_id
phylum_matrix <- phylum_matrix[, colSums(phylum_matrix) > 0]

cat("Phylum matrix:", nrow(phylum_matrix), "samples x", ncol(phylum_matrix), "phyla\n")

# --- 5. CLR-transform (same compositional approach used for beta diversity
#         elsewhere in this analysis) ---
matrix_nozero <- cmultRepl(phylum_matrix, method = "CZM", output = "p-counts")
clr_matrix <- as.matrix(clr(matrix_nozero))

# --- 6. Build the two distance matrices ---
# Compositional (Aitchison) distance between samples' microbiomes
compositional_dist <- dist(clr_matrix, method = "euclidean")

# Geographic distance between samples (in metres, using the projected
# Irish Grid coordinates directly - no need to convert to lat/lon for this)
geo_coords <- as.matrix(coords_matched[, c("xcoord", "ycoord")])
rownames(geo_coords) <- coords_matched$Sample_id
geo_dist <- dist(geo_coords, method = "euclidean")

# --- 7. Run the Mantel test ---
mantel_result <- mantel(geo_dist, compositional_dist, method = "spearman", permutations = 999)

cat("\n=== Mantel test: geographic distance vs. phylum-level compositional distance ===\n")
print(mantel_result)

cat("\nInterpretation:\n")
cat("Mantel r =", round(mantel_result$statistic, 3), "\n")
cat("p-value =", mantel_result$signif, "\n")
if (mantel_result$signif < 0.05) {
  cat("=> Significant relationship: geographically closer samples tend to have\n")
  cat("   more similar phylum-level microbiome composition.\n")
} else {
  cat("=> No significant relationship detected between geographic distance and\n")
  cat("   phylum-level microbiome composition.\n")
}

##############################################################################

## Plasmid replicon type and geographic proximity

library(dplyr)
library(tidyr)
library(vegan)

# --- 1. Load coordinates ---
coords <- read.csv("sequenced_samples_x_y_low_removed.csv", stringsAsFactors = FALSE)
coords$xcoord <- suppressWarnings(as.numeric(coords$xcoord))
coords$ycoord <- suppressWarnings(as.numeric(coords$ycoord))
coords <- coords[!is.na(coords$xcoord) & !is.na(coords$ycoord), ]
coords$Sample_id <- as.integer(coords$Sample)

# --- 2. Load ALL plasmid replicon detections (regardless of ARG association) ---
all_plasmids <- read.csv("Plasmids_on_assemblies.csv", stringsAsFactors = FALSE)
# (update filename above to match your actual file)
all_plasmids$Sample_id <- as.integer(gsub("\\D", "", as.character(all_plasmids$Sample)))

cat("Total plasmid replicon detections:", nrow(all_plasmids), "\n")

# --- Check for and remove any rows with blank/missing GENE values, which
#     would otherwise break pivot_wider (can't create a column named "") ---
n_blank <- sum(all_plasmids$GENE == "" | is.na(all_plasmids$GENE))
if (n_blank > 0) {
  cat("WARNING:", n_blank, "row(s) have a blank/missing GENE value - excluding these\n")
  print(all_plasmids[all_plasmids$GENE == "" | is.na(all_plasmids$GENE),
                     c("Sample", "SEQUENCE", "GENE")])
  all_plasmids <- all_plasmids[all_plasmids$GENE != "" & !is.na(all_plasmids$GENE), ]
}

cat("Plasmid replicon detections after cleaning:", nrow(all_plasmids), "\n")
cat("Unique samples with at least one plasmid replicon detected:",
    n_distinct(all_plasmids$Sample_id), "\n")
cat("Unique plasmid replicon types detected:", n_distinct(all_plasmids$GENE), "\n")

# --- 3. Build sample x replicon-type presence/absence matrix ---
# (each row already represents one replicon detection - no need to split
# a semicolon-separated list, unlike the ARG-linked file)
presence_matrix <- all_plasmids %>%
  distinct(Sample_id, GENE) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = GENE, values_from = present, values_fill = 0)

sample_ids_matrix <- presence_matrix$Sample_id
presence_matrix <- as.matrix(presence_matrix[, -1])
rownames(presence_matrix) <- sample_ids_matrix

cat("\nPresence/absence matrix:", nrow(presence_matrix), "samples x",
    ncol(presence_matrix), "replicon types\n")

# --- 4. Restrict to samples present in BOTH datasets, align order exactly ---
common_ids <- intersect(coords$Sample_id, sample_ids_matrix)
cat("\nSamples with both coordinates and plasmid replicon data:", length(common_ids), "\n")

missing_from_coords <- setdiff(sample_ids_matrix, coords$Sample_id)
if (length(missing_from_coords) > 0) {
  cat("Plasmid-positive samples with NO coordinate data (excluded):\n")
  print(missing_from_coords)
}

if (length(common_ids) < 4) {
  stop("Fewer than 4 samples have both coordinates and plasmid data - a Mantel ",
       "test is not meaningful with this few samples. Check ID matching above.")
}

coords_matched <- coords[match(common_ids, coords$Sample_id), ]
presence_matched <- presence_matrix[match(common_ids, rownames(presence_matrix)), ]

cat("Sample order aligned correctly:",
    all(coords_matched$Sample_id == rownames(presence_matched)), "\n")

# --- 5. Build the two distance matrices ---
plasmid_dist <- vegdist(presence_matched, method = "jaccard", binary = TRUE)

geo_coords <- as.matrix(coords_matched[, c("xcoord", "ycoord")])
rownames(geo_coords) <- coords_matched$Sample_id
geo_dist <- dist(geo_coords, method = "euclidean")

# --- 6. Run the Mantel test ---
mantel_result <- mantel(geo_dist, plasmid_dist, method = "spearman", permutations = 999)

cat("\n=== Mantel test: geographic distance vs. ALL plasmid replicon type profiles ===\n")
print(mantel_result)

cat("\nInterpretation:\n")
cat("Mantel r =", round(mantel_result$statistic, 3), "\n")
cat("p-value =", mantel_result$signif, "\n")
if (mantel_result$signif < 0.05) {
  cat("=> Significant relationship: samples sharing similar plasmid replicon\n")
  cat("   type profiles tend to be geographically closer together.\n")
} else {

}

write.csv(as.data.frame(as.matrix(plasmid_dist)), "all_plasmid_jaccard_distance_matrix.csv")


# --- Coordinates for the four samples of interest ---
samples_of_interest <- data.frame(
  Sample = c("131002", "250004", "266001", "267005"),
  xcoord = c(242507, 201612.56, 107347.7, 239046),
  ycoord = c(240828, 329130.88, 173421.3, 302089)
)

# --- Calculate pairwise Euclidean distance (in metres, since Irish Grid
#     coordinates are already in metres) ---
coords_matrix <- as.matrix(samples_of_interest[, c("xcoord", "ycoord")])
rownames(coords_matrix) <- samples_of_interest$Sample

dist_matrix <- as.matrix(dist(coords_matrix, method = "euclidean"))

cat("=== Pairwise distances (metres) ===\n")
print(round(dist_matrix, 0))

cat("\n=== Pairwise distances (km), easier to interpret ===\n")
print(round(dist_matrix / 1000, 1))

# --- Print as a simple, readable list too ---
cat("\n=== Distance list ===\n")
pairs <- combn(samples_of_interest$Sample, 2)
for (i in seq_len(ncol(pairs))) {
  s1 <- pairs[1, i]
  s2 <- pairs[2, i]
  d_km <- round(dist_matrix[s1, s2] / 1000, 1)
  cat(s1, "-", s2, ":", d_km, "km\n")
}

# For context: what's a "close" vs "far" distance for Ireland?
# Ireland is roughly 480km at its longest (north-south) and 275km at its
# widest (east-west), for reference.
cat("\n(For context: Ireland is roughly 480 km north-south, 275 km east-west)\n")


################################################################################


# Link virulence gene hits (from combined ABRicate/vfdb file) to genus-level
# contig assignments (from the per-sample *_contig_genus_assignments.csv
# files produced by map_contigs_to_genus.R)

library(dplyr)

# --- 1. Load the combined virulence gene file ---
virulence <- read.csv("All_VGs_Assemblies.csv", header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)

cat("Virulence file:", nrow(virulence), "rows,", length(unique(virulence$Sample)), "unique samples\n")
cat("Sample ID format example:", head(unique(virulence$Sample), 3), "\n")

# --- 2. Load and combine all genus assignment CSVs ---
genus_files <- list.files(".", pattern = "_contig_genus_assignments.csv$", full.names = TRUE)
cat("\nFound", length(genus_files), "genus assignment files\n")

genus_all <- lapply(genus_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  # Extract sample ID from filename, e.g. "0131003_contig_genus_assignments.csv" -> "0131003"
  df$sample <- sub("_contig_genus_assignments.csv$", "", basename(f))
  df
})
genus_all <- do.call(rbind, genus_all)

cat("Combined genus table:", nrow(genus_all), "contigs across",
    length(unique(genus_all$sample)), "samples\n")
cat("Sample ID format example:", head(unique(genus_all$sample), 3), "\n")

# --- 2b. Standardise sample ID format ---
# virulence$Sample was read as plain integers (e.g. 97005), while genus_all$sample
# is text with leading zeros preserved from filenames (e.g. "0097005"). Strip
# leading zeros from genus_all$sample and convert both to plain integers so
# they match correctly.
genus_all$sample <- as.integer(genus_all$sample)
virulence$Sample <- as.integer(virulence$Sample)

# --- 3. IMPORTANT: check sample ID formats match before joining ---
# If these don't overlap, the join below will silently produce no matches -
# check the printed examples above (e.g. "97005" vs "0097005") and adjust
# with sub()/sprintf() as needed before proceeding.
sample_overlap <- intersect(unique(virulence$Sample), unique(genus_all$sample))
cat("\nSamples present in BOTH files:", length(sample_overlap), "\n")
if (length(sample_overlap) == 0) {
  cat("WARNING: no matching sample IDs found - check ID format (e.g. leading zeros)\n")
}



# --- 4. Join virulence genes to genus assignments ---
# Matches on sample AND contig ID (virulence "SEQUENCE" = genus table "contig_id")
linked <- virulence %>%
  inner_join(genus_all, by = c("Sample" = "sample", "SEQUENCE" = "contig_id"))

cat("\nLinked table:", nrow(linked), "virulence gene hits successfully matched to a contig\n")
cat("Of these,", sum(!is.na(linked$genus_name)), "have a genus-level taxonomic assignment\n")

write.csv(linked, "virulence_genes_with_genus.csv", row.names = FALSE)
cat("\nSaved to virulence_genes_with_genus.csv\n")

# --- 5. Per-sample virulence gene hit count, INCLUDING samples with zero hits ---
# This makes "no virulence genes detected" explicit for a sample, rather
# than that sample simply being absent from the linked table.
all_samples <- unique(genus_all$sample)  # your full set of 55 samples

hits_per_sample <- virulence %>%
  group_by(Sample) %>%
  summarise(n_virulence_hits = n(), n_unique_genes = n_distinct(GENE))

sample_virulence_summary <- data.frame(sample = all_samples) %>%
  left_join(hits_per_sample, by = c("sample" = "Sample")) %>%
  mutate(
    n_virulence_hits = ifelse(is.na(n_virulence_hits), 0, n_virulence_hits),
    n_unique_genes = ifelse(is.na(n_unique_genes), 0, n_unique_genes)
  ) %>%
  arrange(n_virulence_hits)

write.csv(sample_virulence_summary, "virulence_hits_per_sample_all55.csv", row.names = FALSE)

cat("\n=== Virulence gene detection across all", length(all_samples), "samples ===\n")
cat(sum(sample_virulence_summary$n_virulence_hits == 0), "samples had ZERO virulence gene hits\n")
cat(sum(sample_virulence_summary$n_virulence_hits > 0), "samples had at least one hit\n\n")
print(sample_virulence_summary)

# --- 6. Which genera carry the most virulence gene hits? ---
genus_summary <- linked %>%
  filter(!is.na(genus_name)) %>%
  group_by(genus_name) %>%
  summarise(n_virulence_hits = n(), n_unique_genes = n_distinct(GENE),
            n_samples = n_distinct(Sample)) %>%
  arrange(desc(n_virulence_hits))

write.csv(genus_summary, "virulence_by_genus_summary.csv", row.names = FALSE)
cat("\nTop genera by virulence gene hit count:\n")
print(head(genus_summary, 15))


mean(genus_mapping_summary$pct_with_genus)
range(genus_mapping_summary$pct_with_genus)
length(unique(genus_all$genus_name[!is.na(genus_all$genus_name)]))


median(genus_mapping_summary$total_contigs)


genus_prevalence <- genus_all %>%
  filter(!is.na(genus_name), genus_name != "Homo") %>%
  distinct(sample, genus_name) %>%
  count(genus_name, sort = TRUE, name = "n_samples")

head(genus_prevalence, 15)

sort(unique(genus_all$genus_name))


genus_prevalence$pct_samples <- round(genus_prevalence$n_samples / 55 * 100, 1)
head(genus_prevalence, 50)



genus_prevalence$pct_samples <- round(genus_prevalence$n_samples / 55 * 100, 1)
write.csv(head(genus_prevalence, 50), "genus_prevalence_top50.csv", row.names = FALSE)


## - Figure 4.5

# Heatmap: top virulence genes x top genera (by sample prevalence)

library(ggplot2)
library(dplyr)

virulence_linked <- read.csv("virulence_genes_with_genus.csv", stringsAsFactors = FALSE)

# --- 1. Select top genera by sample prevalence ---
n_top_genera <- 14  # includes essentially all genera with any virulence gene hit
n_genes_per_genus <- 8  # top genes to pull FROM EACH genus, not overall

top_genera <- virulence_linked %>%
  filter(!is.na(genus_name)) %>%
  distinct(Sample, genus_name) %>%
  count(genus_name, sort = TRUE, name = "n_samples") %>%
  slice_head(n = n_top_genera) %>%
  pull(genus_name)

# For each top genus, find ITS OWN top genes by sample prevalence -
# ensures every genus is represented, rather than the gene list being
# dominated by whichever genus happens to have the most gene diversity
top_genes <- virulence_linked %>%
  filter(genus_name %in% top_genera) %>%
  distinct(Sample, GENE, genus_name) %>%
  count(genus_name, GENE, name = "n_samples") %>%
  group_by(genus_name) %>%
  slice_max(n_samples, n = n_genes_per_genus, with_ties = FALSE) %>%
  ungroup() %>%
  pull(GENE) %>%
  unique()

# --- Force-include specific genes of known interest, regardless of
#     prevalence ranking - otherwise rare-but-important genes (e.g. stx2,
#     hly - only 1 sample) get excluded by very common Escherichia genes
#     (30+ samples) when only selecting "top N by prevalence" ---
must_include_genes <- c("stx2A", "stx2B", "hlyA", "hlyB", "hlyC", "hlyD")
top_genes <- union(top_genes, must_include_genes)


# --- 2. Build gene x genus matrix: number of samples where this
#         gene-genus combination occurs ---
heatmap_data <- virulence_linked %>%
  filter(GENE %in% top_genes, genus_name %in% top_genera) %>%
  distinct(Sample, GENE, genus_name) %>%
  count(GENE, genus_name, name = "n_samples")

# Fill in zero combinations (so the heatmap shows blanks, not gaps)
full_grid <- expand.grid(GENE = top_genes, genus_name = top_genera, stringsAsFactors = FALSE)
heatmap_data <- full_grid %>%
  left_join(heatmap_data, by = c("GENE", "genus_name")) %>%
  mutate(n_samples = ifelse(is.na(n_samples), 0, n_samples))

# Order genera by overall prevalence; order genes by which genus they were
# selected for (keeps each genus's genes grouped together visually)
gene_order <- virulence_linked %>%
  filter(GENE %in% top_genes, genus_name %in% top_genera) %>%
  distinct(GENE, genus_name) %>%
  count(GENE, genus_name, name = "n") %>%
  group_by(GENE) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(match(genus_name, top_genera), desc(n)) %>%
  pull(GENE)

heatmap_data$GENE <- factor(heatmap_data$GENE, levels = rev(gene_order))
heatmap_data$genus_name <- factor(heatmap_data$genus_name, levels = top_genera)

# --- 3. Plot ---
p <- ggplot(heatmap_data, aes(x = genus_name, y = GENE, fill = n_samples)) +
  geom_tile(color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_samples > 0, n_samples, "")),
            size = 3, color = ifelse(heatmap_data$n_samples > 15, "white", "black")) +
  scale_fill_gradient(low = "white", high = "#1B4D3E", name = "No. of\nsamples",
                      trans = "sqrt") +
  labs(
    title = "Virulence gene carriage by bacterial genus",
    x = "Genus", y = "Virulence gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
    axis.text.y = element_text(face = "italic", size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.title.position = "plot",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

ggsave("virulence_gene_genus_heatmap.png", p, width = 11, height = 13, dpi = 400)
ggsave("virulence_gene_genus_heatmap.pdf", p, width = 11, height = 13)
print(p)



# Compare virulence gene burden between badger and fox samples

library(dplyr)

virulence_linked <- read.csv("virulence_genes_with_genus.csv", stringsAsFactors = FALSE)

# --- 1. Per-sample summary: total hits, unique genes, host species ---
sample_summary <- virulence_linked %>%
  group_by(Sample, Species) %>%
  summarise(n_hits = n(), n_unique_genes = n_distinct(GENE), .groups = "drop")

cat("Samples with virulence gene data:", nrow(sample_summary), "\n")
print(table(sample_summary$Species))

# --- 2. Add back the zero-hit samples (not present in virulence_linked
#         at all, since it only contains rows where a hit occurred) ---
zero_hit_samples <- data.frame(
  Sample = c(86002, 97006, 121002, 149008, 155005, 164005, 178008,
             209005, 209006, 215005, 218003, 226003, 226007, 227007,
             227010, 238003),
  Species = c("Badger", "Badger", "Badger", "Badger", "Badger", "Badger",
              "Badger", "Badger", "Badger", "Badger", "Fox", "Fox", "Fox",
              "Badger", "Badger", "Fox"),
  n_hits = 0,
  n_unique_genes = 0
)

sample_summary$Sample <- as.integer(sample_summary$Sample)

sample_summary <- bind_rows(sample_summary, zero_hit_samples)

cat("\nFull sample set after adding zero-hit samples:", nrow(sample_summary), "\n")
print(table(sample_summary$Species))

# --- 3. Compare burden between species ---
cat("\n--- Total virulence gene hits ---\n")
print(wilcox.test(n_hits ~ Species, data = sample_summary, exact = FALSE))

cat("\n--- Unique virulence genes ---\n")
print(wilcox.test(n_unique_genes ~ Species, data = sample_summary, exact = FALSE))

cat("\n--- Summary stats ---\n")
print(tapply(sample_summary$n_hits, sample_summary$Species, summary))
print(tapply(sample_summary$n_unique_genes, sample_summary$Species, summary))

write.csv(sample_summary, "virulence_burden_by_species.csv", row.names = FALSE)



# Check the Yersinia and Salmonella virulence gene results more closely:
# 1. Which sample(s) do these hits come from?
# 2. Do the SAME genes also appear attributed to OTHER genera elsewhere in
#    the dataset? (VFDB genes are often shared/conserved across related
#    Enterobacteriaceae, so a gene could plausibly appear on contigs
#    classified to different genera in different samples)

library(dplyr)

virulence_linked <- read.csv("virulence_genes_with_genus.csv", stringsAsFactors = FALSE)

# --- 1. Which samples carry Yersinia / Salmonella virulence hits? ---
cat("=== Yersinia hits: which sample(s)? ===\n")
yersinia_hits <- virulence_linked %>% filter(genus_name == "Yersinia")
print(table(yersinia_hits$Sample))

cat("\n=== Salmonella hits: which sample(s)? ===\n")
salmonella_hits <- virulence_linked %>% filter(genus_name == "Salmonella")
print(table(salmonella_hits$Sample))

# --- 2. For each Yersinia/Salmonella gene, check if it ALSO appears under
#         a different genus anywhere else in the dataset ---
yersinia_genes <- unique(yersinia_hits$GENE)
salmonella_genes <- unique(salmonella_hits$GENE)

cat("\n\n=== For each Yersinia gene, which genera is it seen under overall? ===\n")
for (g in yersinia_genes) {
  genera_for_gene <- virulence_linked %>%
    filter(GENE == g) %>%
    distinct(genus_name) %>%
    pull(genus_name)
  cat(g, ":", paste(genera_for_gene, collapse = ", "), "\n")
}

cat("\n\n=== For each Salmonella gene, which genera is it seen under overall? ===\n")
for (g in salmonella_genes) {
  genera_for_gene <- virulence_linked %>%
    filter(GENE == g) %>%
    distinct(genus_name) %>%
    pull(genus_name)
  cat(g, ":", paste(genera_for_gene, collapse = ", "), "\n")
}

# --- 3. Summary: how many of these genes are genus-EXCLUSIVE vs shared ---
cat("\n\n=== Summary ===\n")
yersinia_exclusive <- sum(sapply(yersinia_genes, function(g) {
  n_distinct(virulence_linked$genus_name[virulence_linked$GENE == g]) == 1
}))
cat("Yersinia genes found ONLY under Yersinia (not shared with other genera):",
    yersinia_exclusive, "of", length(yersinia_genes), "\n")

salmonella_exclusive <- sum(sapply(salmonella_genes, function(g) {
  n_distinct(virulence_linked$genus_name[virulence_linked$GENE == g]) == 1
}))
cat("Salmonella genes found ONLY under Salmonella (not shared with other genera):",
    salmonella_exclusive, "of", length(salmonella_genes), "\n")




# Example: Extract a representative sample of contigs carrying Yersinia-associated
# virulence genes in sample X, for BLAST verification.

library(stringr)
library(dplyr)

sample <- "sampleX"

virulence_linked <- read.csv("virulence_genes_with_genus.csv", stringsAsFactors = FALSE)

yersinia_contigs <- virulence_linked %>%
  filter(genus_name == "Yersinia", Sample %in% c("sampleX", sampleX)) %>%
  distinct(SEQUENCE, GENE)

cat("Unique contigs carrying Yersinia genes in this sample:", n_distinct(yersinia_contigs$SEQUENCE), "\n")
print(yersinia_contigs)

# --- Pick a representative subset: a few different contigs, ideally
# covering different gene categories (flagellar vs O-antigen) ---
# Adjust this selection once you see the full list above
target_contigs <- unique(yersinia_contigs$SEQUENCE)[1:5]  # first 5 as a starting point
cat("\nSelected for extraction:", paste(target_contigs, collapse = ", "), "\n")

# --- Extract these contigs from the assembly FASTA ---
assembly_path <- "C:/PATH/TO/assembly.fasta"

if (!file.exists(assembly_path)) {
  cat("\nAssembly FASTA not found at:", assembly_path, "\n")
  cat("Update the 'assembly_path' variable to the correct location.\n")
} else {
  fasta_lines <- readLines(assembly_path)
  
  header_idx <- which(str_starts(fasta_lines, ">"))
  header_names <- str_remove(fasta_lines[header_idx], "^>") %>% str_split(" ", simplify = TRUE)
  header_names <- header_names[, 1]
  
  extracted <- character(0)
  for (cid in target_contigs) {
    match_idx <- which(header_names == cid)
    if (length(match_idx) == 0) {
      cat("WARNING: contig", cid, "not found in FASTA headers\n")
      next
    }
    start_line <- header_idx[match_idx]
    end_line <- ifelse(match_idx < length(header_idx), header_idx[match_idx + 1] - 1, length(fasta_lines))
    extracted <- c(extracted, fasta_lines[start_line:end_line])
  }
  
  writeLines(extracted, paste0(sample, "_Yersinia_contigs.fasta"))
}



## Figure S4.3 

# Presence/absence heatmap: top genera by prevalence, across all samples,
# grouped by host species (badger/fox)

library(ggplot2)
library(dplyr)

# --- 1. Build genus x sample presence/absence data ---
# Reuses genus_all (all contig-genus assignments) - adjust if loading fresh
genus_sample_list <- genus_all %>%
  filter(!is.na(genus_name), genus_name != "Homo") %>%
  distinct(sample, genus_name)

# --- 2. Select top N most prevalent genera ---
n_top <- 100

top_genera <- genus_sample_list %>%
  count(genus_name, sort = TRUE, name = "n_samples") %>%
  slice_head(n = n_top) %>%
  pull(genus_name)

# --- 3. Get the FULL list of 55 samples (so absent samples show as blank,
#         not just missing from the plot) ---
all_samples <- unique(genus_all$sample)

# --- 4. Build a complete grid: every sample x every top genus, filled in
#         with presence (1) or absence (0) ---
heatmap_data <- expand.grid(sample = all_samples, genus_name = top_genera,
                            stringsAsFactors = FALSE) %>%
  left_join(genus_sample_list %>% mutate(present = 1),
            by = c("sample", "genus_name")) %>%
  mutate(present = ifelse(is.na(present), 0, present))

# --- 5. Add host species, and order samples by species then ID ---
heatmap_data$host_species <- ifelse(grepl("fox", heatmap_data$sample, ignore.case = TRUE), "Fox",
                                    ifelse(grepl("badger", heatmap_data$sample, ignore.case = TRUE), "Badger",
                                           # If your sample IDs don't contain "fox"/"badger" text (likely, based on
                                           # IDs like "125003"), you'll need a proper lookup table instead - see note below
                                           "Unknown"))

# NOTE: your sample IDs (e.g. "125003") likely do NOT contain "fox"/"badger"
# text directly - unlike your bracken filenames did. You'll need to supply a
# sample -> host species mapping here instead. For example, if you have your
# original alpha_div or metadata table with sample + host_species columns:
#
# species_lookup <- alpha_div %>% distinct(sample, host_species)
# (adjust "sample" formatting/matching as needed - e.g. strip leading zeros)
#
# heatmap_data <- heatmap_data %>%
#   select(-host_species) %>%
#   left_join(species_lookup, by = "sample")

sample_order <- heatmap_data %>%
  distinct(sample, host_species) %>%
  arrange(host_species, sample) %>%
  pull(sample)

heatmap_data$sample <- factor(heatmap_data$sample, levels = sample_order)

# Order genera by overall prevalence (most prevalent at top)
heatmap_data$genus_name <- factor(heatmap_data$genus_name, levels = rev(top_genera))

# --- 6. Plot ---
p <- ggplot(heatmap_data, aes(x = sample, y = genus_name, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("0" = "grey90", "1" = "#3B6E5B"),
                    labels = c("Absent", "Present"), name = NULL) +
  facet_grid(~ host_species, scales = "free_x", space = "free_x") +
  labs(title = paste("Presence of top", n_top, "most prevalent genera across samples"),
       x = "Sample", y = "Genus") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
    axis.text.y = element_text(face = "italic", size = 16),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 13)
  )

p

ggsave("genus_presence_heatmap.png", p, width = 30, height = 24, dpi = 300)
ggsave("genus_presence_heatmap.pdf", p, width = 30, height = 24)

print(p)


## -----------------------------------------------------------------------------


# Survey ALL Mycobacterium species detected across all samples

library(stringr)
library(dplyr)

report_files <- list.files(".", pattern = "_kraken2_report.txt$", full.names = FALSE)
cat("Screening", length(report_files), "report files for all Mycobacterium species...\n\n")

all_hits <- list()

for (r in report_files) {
  sample <- sub("_kraken2_report.txt$", "", r)
  
  report <- read.delim(r, header = FALSE, sep = "\t",
                       col.names = c("pct", "reads_clade", "reads_direct",
                                     "rank_code", "taxid", "name_indented"),
                       stringsAsFactors = FALSE)
  report$name <- str_trim(report$name_indented)
  
  # Any row where the name starts with "Mycobacterium" and has an actual
  # hit (reads_clade > 0), at species level or finer (S, S1, S2, S3)
  hits <- report %>%
    dplyr::filter(str_starts(name, "Mycobacterium"),
                  rank_code %in% c("S", "S1", "S2", "S3"),
                  reads_clade > 0)
  
  if (nrow(hits) > 0) {
    hits$sample <- sample
    all_hits[[sample]] <- hits %>%
      dplyr::select(sample, name, taxid, rank_code, reads_clade, reads_direct, pct)
  }
}

if (length(all_hits) == 0) {
  cat("No species-level Mycobacterium hits detected in any sample.\n")
} else {
  myco_results <- do.call(rbind, all_hits)
  rownames(myco_results) <- NULL
  myco_results <- myco_results %>% arrange(name, sample)
  
  cat("Mycobacterium species-level hits found in", length(unique(myco_results$sample)),
      "of", length(report_files), "samples\n")
  cat("Unique species detected:", length(unique(myco_results$name)), "\n\n")
  
  print(myco_results)
  
  write.csv(myco_results, "Mycobacterium_all_species_results.csv", row.names = FALSE)
  cat("\nSaved to Mycobacterium_all_species_results.csv\n")
  
  # Quick summary: which Mycobacterium species are most common, and with
  # what typical read/contig support (helps spot single-read noise vs
  # genuinely well-supported detections)
  species_summary <- myco_results %>%
    group_by(name) %>%
    summarise(n_samples = n_distinct(sample),
              median_reads = median(reads_clade),
              max_reads = max(reads_clade)) %>%
    arrange(desc(n_samples))
  
  cat("\n=== Summary by species ===\n")
  print(species_summary, n = 50)
}




# Screen all 55 samples for Salmonella and Listeria (genera of public
# health importance), at both genus and species level, using the Kraken2
# REPORT files. Reports read/contig support alongside each hit so you can
# immediately judge confidence (same approach used for the MTBC screen).

library(stringr)
library(dplyr)

target_genera <- c("Salmonella", "Listeria")

report_files <- list.files(".", pattern = "_kraken2_report.txt$", full.names = FALSE)
cat("Screening", length(report_files), "report files for:", paste(target_genera, collapse = ", "), "\n\n")

all_hits <- list()

for (r in report_files) {
  sample <- sub("_kraken2_report.txt$", "", r)
  
  report <- read.delim(r, header = FALSE, sep = "\t",
                       col.names = c("pct", "reads_clade", "reads_direct",
                                     "rank_code", "taxid", "name_indented"),
                       stringsAsFactors = FALSE)
  report$name <- str_trim(report$name_indented)
  
  # Genus-level AND species-level (or finer) hits, for either target genus
  hits <- report %>%
    dplyr::filter(
      (str_starts(name, paste(target_genera, collapse = "|"))) &
        (rank_code %in% c("G", "S", "S1", "S2", "S3")) &
        reads_clade > 0
    )
  
  if (nrow(hits) > 0) {
    hits$sample <- sample
    all_hits[[sample]] <- hits %>%
      dplyr::select(sample, name, taxid, rank_code, reads_clade, reads_direct, pct)
  }
}

if (length(all_hits) == 0) {
  cat("No Salmonella or Listeria hits detected in any sample.\n")
} else {
  results <- do.call(rbind, all_hits)
  rownames(results) <- NULL
  results <- results %>% arrange(name, sample)
  
  write.csv(results, "Salmonella_Listeria_screening_results.csv", row.names = FALSE)
  cat("Saved full results to Salmonella_Listeria_screening_results.csv\n\n")
  
  cat("Hits found in", length(unique(results$sample)), "of", length(report_files), "samples\n")
  cat("Unique taxa detected:", length(unique(results$name)), "\n\n")
  
  tryCatch(print(results), error = function(e) {
    cat("(console print failed - open Salmonella_Listeria_screening_results.csv to view results)\n")
  })
  
  # Summary: how many samples per taxon, and typical read support
  summary_table <- results %>%
    group_by(name, rank_code) %>%
    summarise(n_samples = n_distinct(sample),
              median_reads = median(reads_clade),
              max_reads = max(reads_clade),
              .groups = "drop") %>%
    arrange(desc(n_samples))
  
  cat("\n=== Summary by taxon ===\n")
  tryCatch(print(summary_table), error = function(e) {
    cat("(console print failed - check the CSV file for full results)\n")
  })
}


# Extract the contig(s) behind the L. monocytogenes hits in sample X,
# for independent BLAST verification - same approach as the MTBC checks.

library(stringr)
library(dplyr)

sample <- "sampleX"
target_taxids <- c(1639)  # Listeria monocytogenes

# --- 1. Find which contig(s) in the OUTPUT file were assigned to this taxid ---
output_file <- paste0(sample, "_kraken2_output.txt")

contig_output <- read.delim(output_file, header = FALSE, sep = "\t",
                            col.names = c("classified", "contig_id", "taxon_field",
                                          "length", "lca_map"),
                            stringsAsFactors = FALSE)

contig_output$taxid <- as.integer(contig_output$taxon_field)

matching_contigs <- contig_output[contig_output$taxid %in% target_taxids, ]
cat("Contig(s) assigned to L. monocytogenes in sample", sample, ":\n")
print(matching_contigs)

if (nrow(matching_contigs) == 0) {
  stop("No matching contigs found - double check taxid/sample name")
}

# --- 2. Extract these contigs' sequences from the assembly FASTA ---
# UPDATE this path to point to the actual assembly.fasta location for this sample
assembly_path <- "C:/Users/cedoy/Susceptibles/0267005/metaflye/0267005_assembly_output/assembly.fasta"

if (!file.exists(assembly_path)) {
  cat("\nAssembly FASTA not found at:", assembly_path, "\n")
  cat("Update the 'assembly_path' variable to the correct location.\n")
} else {
  fasta_lines <- readLines(assembly_path)
  
  header_idx <- which(str_starts(fasta_lines, ">"))
  header_names <- str_remove(fasta_lines[header_idx], "^>") %>% str_split(" ", simplify = TRUE)
  header_names <- header_names[, 1]
  
  extracted <- character(0)
  for (cid in matching_contigs$contig_id) {
    match_idx <- which(header_names == cid)
    if (length(match_idx) == 0) {
      cat("WARNING: contig", cid, "not found in FASTA headers - check ID format\n")
      next
    }
    start_line <- header_idx[match_idx]
    end_line <- ifelse(match_idx < length(header_idx), header_idx[match_idx + 1] - 1, length(fasta_lines))
    extracted <- c(extracted, fasta_lines[start_line:end_line])
  }
  
  writeLines(extracted, paste0(sample, "_Lmonocytogenes_contigs.fasta"))
  cat("\nSaved extracted sequence(s) to", paste0(sample, "_Lmonocytogenes_contigs.fasta"), "\n")
  cat("This may contain MULTIPLE contigs (up to 8 expected) - BLAST each one\n")
  cat("individually via NCBI BLAST (blastn, nt database) to verify:\n")
  cat("https://blast.ncbi.nlm.nih.gov/Blast.cgi\n")
}

## -----------------------------------------------------------------------------

## PLASMIDS AND ARGs LINKAGE TO IDENTIFIED CONTIGS ##

# 1. Load ARG/plasmid data and combined genus assignments
args_data <- read.csv("arg_plasmid_context_full.csv", stringsAsFactors = FALSE)

genus_files <- list.files(".", pattern = "_contig_genus_assignments.csv$", full.names = TRUE)
genus_all <- lapply(genus_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$sample <- sub("_contig_genus_assignments.csv$", "", basename(f))
  df
})
genus_all <- do.call(rbind, genus_all)

# 2. Standardise sample ID format before joining 

args_data$Sample <- as.integer(args_data$Sample)
genus_all$sample <- as.integer(genus_all$sample)

# 3. Join on sample + contig ID
args_linked <- args_data %>%
  inner_join(genus_all, by = c("Sample" = "sample", "SEQUENCE" = "contig_id"))

cat("ARG hits successfully matched to a contig:", nrow(args_linked), "of", nrow(args_data), "\n")
cat("Of these, with a genus-level assignment:", sum(!is.na(args_linked$genus_name)), "\n\n")

write.csv(args_linked, "args_with_genus.csv", row.names = FALSE)

# 4. Overview numbers
print(table(args_linked$plasmid_associated))

# 5. Top ARG genes, resistance classes, and genera

cat("\n=== Top 15 PLASMID-ASSOCIATED ARG genes by sample prevalence ===\n")
top_args_plasmid <- args_linked %>%
  filter(plasmid_associated == TRUE) %>%
  distinct(Sample, Gene) %>%
  count(Gene, sort = TRUE, name = "n_samples")
print(head(top_args_plasmid, 15))

cat("\n=== Top 15 NON-PLASMID (chromosomal) ARG genes by sample prevalence ===\n")
top_args_nonplasmid <- args_linked %>%
  filter(plasmid_associated == FALSE) %>%
  distinct(Sample, Gene) %>%
  count(Gene, sort = TRUE, name = "n_samples")
print(head(top_args_nonplasmid, 15))

cat("\n=== Top 15 genera carrying PLASMID-ASSOCIATED ARGs, by sample prevalence ===\n")
top_genera_args_plasmid <- args_linked %>%
  filter(!is.na(genus_name), plasmid_associated == TRUE) %>%
  distinct(Sample, genus_name) %>%
  count(genus_name, sort = TRUE, name = "n_samples")
print(head(top_genera_args_plasmid, 15))

cat("\n=== Top 15 genera carrying NON-PLASMID ARGs, by sample prevalence ===\n")
top_genera_args_nonplasmid <- args_linked %>%
  filter(!is.na(genus_name), plasmid_associated == FALSE) %>%
  distinct(Sample, genus_name) %>%
  count(genus_name, sort = TRUE, name = "n_samples")
print(head(top_genera_args_nonplasmid, 15))

# --- 6. Flag anything of particular clinical significance ---
# Worth a manual look even in a lightweight treatment - e.g. carbapenemase
# genes are a standout finding regardless of chapter scope.
high_concern_patterns <- c("KPC", "NDM", "OXA-48", "VIM", "IMP", "mcr-",
                           "blaCTX-M", "blaOXA", "vanA", "vanB")
flagged <- args_linked %>%
  filter(grepl(paste(high_concern_patterns, collapse = "|"), Gene, ignore.case = TRUE))

if (nrow(flagged) > 0) {
  cat("\n=== FLAGGED: genes matching high-clinical-concern patterns ===\n")
  print(flagged %>% distinct(Sample, Gene, genus_name, RESISTANCE))
} else {
  cat("\nNo genes matched common high-clinical-concern patterns (KPC, NDM,",
      "OXA-48, VIM, IMP, mcr, ESBL blaCTX-M/blaOXA, vanA/B).\n")
}

write.csv(top_args_plasmid, "top_ARGs_plasmid_associated.csv", row.names = FALSE)
write.csv(top_args_nonplasmid, "top_ARGs_nonplasmid.csv", row.names = FALSE)
write.csv(top_genera_args_plasmid, "top_genera_ARGs_plasmid_associated.csv", row.names = FALSE)
write.csv(top_genera_args_nonplasmid, "top_genera_ARGs_nonplasmid.csv", row.names = FALSE)

# 7. Same tables, with the actual sample IDs listed
make_sample_list_table <- function(data, group_col) {
  data %>%
    distinct(Sample, .data[[group_col]]) %>%
    arrange(.data[[group_col]], Sample) %>%
    group_by(.data[[group_col]]) %>%
    summarise(n_samples = n(), samples = paste(Sample, collapse = ", "), .groups = "drop") %>%
    arrange(desc(n_samples))
}

args_by_sample_plasmid <- make_sample_list_table(
  args_linked %>% filter(plasmid_associated == TRUE), "Gene")
args_by_sample_nonplasmid <- make_sample_list_table(
  args_linked %>% filter(plasmid_associated == FALSE), "Gene")

genus_by_sample_args_plasmid <- make_sample_list_table(
  args_linked %>% filter(!is.na(genus_name), plasmid_associated == TRUE), "genus_name")
genus_by_sample_args_nonplasmid <- make_sample_list_table(
  args_linked %>% filter(!is.na(genus_name), plasmid_associated == FALSE), "genus_name")

write.csv(args_by_sample_plasmid, "ARGs_plasmid_associated_with_sample_list.csv", row.names = FALSE)
write.csv(args_by_sample_nonplasmid, "ARGs_nonplasmid_with_sample_list.csv", row.names = FALSE)
write.csv(genus_by_sample_args_plasmid, "genera_ARGs_plasmid_associated_with_sample_list.csv", row.names = FALSE)
write.csv(genus_by_sample_args_nonplasmid, "genera_ARGs_nonplasmid_with_sample_list.csv", row.names = FALSE)


## Figure 4.6

# Alluvial (Sankey-style) diagram: Genus -> Resistance class -> Plasmid
# replicon type, for plasmid-associated ARGs only.
# Shows the connections between which bacteria carry which resistance
# genes, and via which plasmid families - richer than a flat table alone.

library(ggalluvial)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1. Prepare data: one row per ARG hit, with genus, resistance class,
#         and (simplified) plasmid replicon type ---
alluvial_data <- args_linked %>%
  filter(plasmid_associated == TRUE, !is.na(genus_name)) %>%
  # RESISTANCE can list multiple classes separated by ";" - take the first
  # for simplicity in this lightweight visualisation
  mutate(resistance_class = sub(";.*", "", RESISTANCE)) %>%
  mutate(resistance_class = paste0(toupper(substr(resistance_class, 1, 1)),
                                   substr(resistance_class, 2, nchar(resistance_class)))) %>%
  # replicon_types_on_contig can list multiple types separated by "; " -
  # take the first for the same reason
  mutate(plasmid_type = sub(";.*", "", replicon_types_on_contig)) %>%
  filter(!is.na(plasmid_type), plasmid_type != "NA", plasmid_type != "") %>%
  count(genus_name, resistance_class, plasmid_type, name = "freq")

cat("Alluvial data built:", nrow(alluvial_data), "unique genus-resistance-plasmid combinations\n")
print(alluvial_data)

# 2. Plot
n_genera <- n_distinct(alluvial_data$genus_name)

# Use a strong, maximally distinct palette rather than ggplot's default hue
# scale, which produces similar muted colours once you have this many levels
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)
genus_palette <- colorRampPalette(brewer.pal(12, "Paired"))(n_genera)
names(genus_palette) <- sort(unique(alluvial_data$genus_name))


p <- ggplot(alluvial_data,
            aes(axis1 = genus_name, axis2 = resistance_class, axis3 = plasmid_type,
                y = freq)) +
  geom_alluvium(aes(fill = genus_name), width = 0.28, alpha = 0.85, color = "grey20", linewidth = 0.1) +
  geom_stratum(width = 0.28, fill = "grey92", color = "grey20", linewidth = 0.4) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 7, fontface = "bold") +
  scale_fill_manual(values = genus_palette) +
  scale_x_discrete(limits = c("Genus", "Resistance class", "Plasmid replicon type"),
                   expand = c(0.12, 0.12)) +
  labs(
    title = "Plasmid-associated ARGs: genus, resistance class, and plasmid type",
    x = NULL, y = "Number of ARG hits"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )

ggsave("ARG_plasmid_alluvial.png", p, width = 15, height = 16, dpi = 300)
ggsave("ARG_plasmid_alluvial.pdf", p, width = 34, height = 22)
print(p)
cat("\nSaved ARG_plasmid_alluvial.png and .pdf\n")



# Genus-by-genus correlation heatmap, based on ARG hit counts per sample
# (full, unfiltered ARG dataset - not just plasmid-associated).

library(corrplot)
library(dplyr)
library(tidyr)

# --- 1. Build a genus x sample matrix of ARG hit counts ---
#

genus_sample_matrix <- args_linked %>%
  filter(!is.na(genus_name)) %>%
  count(Sample, genus_name, name = "n_hits") %>%
  pivot_wider(names_from = genus_name, values_from = n_hits, values_fill = 0)

sample_ids <- genus_sample_matrix$Sample
genus_sample_matrix <- as.matrix(genus_sample_matrix[, -1])
rownames(genus_sample_matrix) <- sample_ids

cat("Matrix built:", nrow(genus_sample_matrix), "samples x",
    ncol(genus_sample_matrix), "genera\n")

# --- 2. Filter out genera present in too few samples to give a stable
#         correlation estimate (a common-sense minimum, e.g. >=3 samples) ---
min_samples <- 2
genus_prevalence <- colSums(genus_sample_matrix > 0)
keep_genera <- names(genus_prevalence[genus_prevalence >= min_samples])

cat("Keeping", length(keep_genera), "of", ncol(genus_sample_matrix),
    "genera with ARG hits in at least", min_samples, "samples\n")

genus_sample_matrix <- genus_sample_matrix[, keep_genera]

# --- 3. Calculate genus-by-genus correlation matrix ---
cor_matrix <- cor(genus_sample_matrix, method = "spearman")

# --- 4. Plot using corrplot, hierarchically clustered for readability ---
png("ARG_genus_correlation.png", width = 2400, height = 2400, res = 300)
corrplot(cor_matrix, method = "color", order = "hclust",
         tl.col = "black", tl.cex = 0.7, tl.srt = 45,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         title = "Genus-genus correlation in ARG carriage across samples",
         mar = c(0, 0, 2, 0))
dev.off()



# Build a detailed genus-level ARG summary table: prevalence, gene
# diversity, total hits, and dominant resistance class(es)

library(dplyr)

# Reload if not already in session:
# args_linked <- read.csv("args_with_genus.csv")

# --- Simplify RESISTANCE to its first-listed class, for a clean "dominant
#     class" summary (same simplification used in the alluvial plot) ---
args_linked <- args_linked %>%
  mutate(resistance_class = sub(";.*", "", RESISTANCE)) %>%
  mutate(resistance_class = paste0(toupper(substr(resistance_class, 1, 1)),
                                   substr(resistance_class, 2, nchar(resistance_class))))

genus_arg_summary <- args_linked %>%
  filter(!is.na(genus_name)) %>%
  group_by(genus_name) %>%
  summarise(
    n_samples = n_distinct(Sample),
    n_unique_ARGs = n_distinct(Gene),
    total_ARG_hits = n(),
    dominant_resistance_class = names(sort(table(resistance_class), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples))

cat("Genus-level ARG summary (", nrow(genus_arg_summary), "genera):\n\n")
print(genus_arg_summary, n = 100)

write.csv(genus_arg_summary, "genus_ARG_summary_table.csv", row.names = FALSE)
cat("\nSaved to genus_ARG_summary_table.csv\n")


# Check whether any contigs carry BOTH a plasmid-associated ARG AND a
# virulence gene - i.e. do these two datasets share any contigs in common?

library(dplyr)

# --- 1. Load both datasets ---
virulence_linked <- read.csv("virulence_genes_with_genus.csv", stringsAsFactors = FALSE)

# Reuse args_linked if already in session, otherwise reload:
# args_linked <- read.csv("args_with_genus.csv")

plasmid_args <- args_linked %>% filter(plasmid_associated == TRUE)

# --- 2. Standardise sample ID format on both sides (same fix as before) ---
plasmid_args$Sample <- as.integer(plasmid_args$Sample)
virulence_linked$Sample <- as.integer(virulence_linked$Sample)

# --- 3. Find contigs present in BOTH datasets (same Sample + same contig) ---
shared_contigs <- plasmid_args %>%
  distinct(Sample, SEQUENCE) %>%
  inner_join(virulence_linked %>% distinct(Sample, SEQUENCE), by = c("Sample", "SEQUENCE"))

cat("Contigs carrying BOTH a plasmid-associated ARG AND a virulence gene:",
    nrow(shared_contigs), "\n\n")

if (nrow(shared_contigs) > 0) {
  print(shared_contigs)
  
  # --- 4. Pull the full detail for these specific shared contigs, from
  #         both datasets, so you can see exactly what's co-occurring ---
  cat("\n=== ARG detail for shared contigs ===\n")
  arg_detail <- plasmid_args %>%
    inner_join(shared_contigs, by = c("Sample", "SEQUENCE")) %>%
    select(Sample, SEQUENCE, Gene, RESISTANCE, replicon_types_on_contig, genus_name)
  print(arg_detail)
  
  cat("\n=== Virulence gene detail for shared contigs ===\n")
  vir_detail <- virulence_linked %>%
    inner_join(shared_contigs, by = c("Sample", "SEQUENCE")) %>%
    select(Sample, SEQUENCE, GENE, PRODUCT, genus_name)
  print(vir_detail)
  
  write.csv(arg_detail, "shared_contigs_ARG_detail.csv", row.names = FALSE)
  write.csv(vir_detail, "shared_contigs_virulence_detail.csv", row.names = FALSE)
  cat("\nSaved shared_contigs_ARG_detail.csv and shared_contigs_virulence_detail.csv\n")
} else {
  cat("No overlap found - no contig carries both a plasmid-associated ARG\n")
  cat("and a virulence gene in this dataset.\n")
}

## -----------------------------------------------------------------------------

# Extract the FULL contig_64 sequence from sampleX, to check via
# BLASTn whether this contig (carrying the normally-chromosomal ampC gene,
# but flagged as plasmid-associated via Col440I/ColRNAI replicon markers)
# genuinely resembles a known plasmid.


sample <- "sampleX"
target_contig <- "contig_X"

# UPDATE this path to point to the actual assembly.fasta location
assembly_path <- "C:/PATH/TO/assembly.fasta"

if (!file.exists(assembly_path)) {
  stop(paste("Assembly FASTA not found at:", assembly_path,
             "- update the 'assembly_path' variable at the top of the script"))
}

fasta_lines <- readLines(assembly_path)
header_idx <- which(str_starts(fasta_lines, ">"))
header_names <- str_remove(fasta_lines[header_idx], "^>") %>% str_split(" ", simplify = TRUE)
header_names <- header_names[, 1]

match_idx <- which(header_names == target_contig)
if (length(match_idx) == 0) {
  stop(paste("Contig", target_contig, "not found in FASTA headers - check ID format"))
}

start_line <- header_idx[match_idx]
end_line <- ifelse(match_idx < length(header_idx), header_idx[match_idx + 1] - 1, length(fasta_lines))
seq_lines <- fasta_lines[(start_line + 1):end_line]
full_seq <- paste(seq_lines, collapse = "")

cat("Contig", target_contig, "length:", nchar(full_seq), "bp\n")

# --- Save the full contig ---
out_file <- paste0(sample, "_", target_contig, "_full.fasta")
writeLines(c(paste0(">", sample, "_", target_contig), full_seq), out_file)
cat("Saved full contig to", out_file, "\n")



