library(tidyverse)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(vegan)
library(dplyr)
install.packages("rstatix", type = "binary", dependencies = TRUE)
library(rstatix)
library(RColorBrewer)

## WGS validation check

validation <-read.csv("WGS_isolates.csv")
validation$Total_reads <- as.numeric(gsub(",", "", validation$Total_reads))
validation$Depth <- as.numeric(validation$Depth)
validation$Breadth_coverage <- as.numeric(validation$Breadth_coverage)
str(validation)


## Reads vs breadth coverage
cor.test(validation$Total_reads,
         validation$Breadth_coverage,
         method = "spearman")

cor.test(validation$Total_reads,
         validation$Depth,
         method = "spearman")

## -----------------------------------------------------------------------------

## Figure 3.3 - Map

# 1. Load an outline map of the Republic of Ireland
ireland_full <- ne_countries(country = "Ireland", scale = "medium", returnclass = "sf")

# 2. Load sample coordinates
samples <- read.csv("sequenced_samples_x_y.csv", stringsAsFactors = FALSE)

cat("Loaded", nrow(samples), "samples\n")

# Check for missing/non-numeric coordinates BEFORE converting
samples$xcoord_numeric <- suppressWarnings(as.numeric(samples$xcoord))
samples$ycoord_numeric <- suppressWarnings(as.numeric(samples$ycoord))

bad_rows <- samples[is.na(samples$xcoord_numeric) | is.na(samples$ycoord_numeric), ]

if (nrow(bad_rows) > 0) {
  cat("\nWARNING:", nrow(bad_rows), "row(s) have missing or non-numeric coordinates:\n")
  print(bad_rows[, c("Sample", "xcoord", "ycoord", "SPECIES")])
  cat("\nSpecies breakdown of EXCLUDED samples:\n")
  print(table(bad_rows$SPECIES))
  cat("\nThese rows will be excluded from the map. Fix the source CSV if this\n")
  cat("is unexpected (e.g. check for blank cells, stray text, or commas used\n")
  cat("as decimal separators instead of periods).\n\n")
}

# Remove bad rows before proceeding
samples <- samples[!is.na(samples$xcoord_numeric) & !is.na(samples$ycoord_numeric), ]
samples$xcoord <- samples$xcoord_numeric
samples$ycoord <- samples$ycoord_numeric
samples$xcoord_numeric <- NULL
samples$ycoord_numeric <- NULL

cat("Proceeding with", nrow(samples), "samples with valid coordinates\n")
cat("Species breakdown of samples being plotted:\n")
print(table(samples$SPECIES))
cat("Coordinate range - X:", range(samples$xcoord), "| Y:", range(samples$ycoord), "\n")

#    Add a small jitter to ALL sample coordinates, so overlapping or
#     closely-spaced points are easier to distinguish on the map. Jitter
#     is small (+/- up to 500m) relative to national map scale, so it
#     does not meaningfully misrepresent sample location. ---
set.seed(42)  # reproducible jitter
n <- nrow(samples)
samples$xcoord <- samples$xcoord + runif(n, -500, 500)
samples$ycoord <- samples$ycoord + runif(n, -500, 500)

# 3. Convert to a spatial object using the Irish Grid CRS (EPSG:29903)
samples_sf <- st_as_sf(samples, coords = c("xcoord", "ycoord"), crs = 29903)

# Reproject to WGS84 (standard lat/lon) so it lines up with the Ireland
# outline map, which uses that system by default
samples_sf <- st_transform(samples_sf, crs = 4326)

# 4. Plot

# Main map: full extent
p_main <- ggplot() +
  geom_sf(data = ireland_full, fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_sf(data = samples_sf, aes(color = SPECIES), size = 3, alpha = 1) +
  scale_color_manual(values = c("Badger" = "#3B6E5B", "Fox" = "#D97D26")) +
  labs(title = "Sample locations", color = "Species") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 18),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin(5, 2, 5, 5)
  )

# --- Identify the tight cluster EXPLICITLY (the automatic distance-based
#     approach chained together too many points across a wider area) ---
cluster_samples <- c("218003", "226002", "226003")

cat("\nZooming in on these specific samples:\n")
print(samples[samples$Sample %in% cluster_samples, c("Sample", "SPECIES")])

# Zoomed bounding box around just these points, with a small buffer
cluster_sf <- samples_sf[samples$Sample %in% cluster_samples, ]
cluster_bbox <- st_bbox(st_buffer(st_transform(cluster_sf, crs = 29903), dist = 1500))
cluster_bbox_wgs84 <- st_bbox(st_transform(st_as_sfc(cluster_bbox, crs = 29903), crs = 4326))

p_inset <- ggplot() +
  geom_sf(data = ireland_full, fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_sf(data = samples_sf, aes(color = SPECIES), size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Badger" = "#3B6E5B", "Fox" = "#D97D26")) +
  coord_sf(
    xlim = c(cluster_bbox_wgs84["xmin"], cluster_bbox_wgs84["xmax"]),
    ylim = c(cluster_bbox_wgs84["ymin"], cluster_bbox_wgs84["ymax"])
  ) +
  labs(title = "Zoomed: clustered samples") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(5, 5, 5, 2)
  )

# Combine
p <- p_main + p_inset + plot_layout(widths = c(2.2, 1))

print(p)

## -----------------------------------------------------------------------------

## Figure 3.6 Comparing filtering workflows - unfiltered v filtered ######

## 1.Load and label each workflow 

raw       <- read_csv("raw_unfiltered.csv") %>% mutate(Workflow = "Raw unfiltered")
raw_kma   <- read_csv("raw_filtered.csv")   %>% mutate(Workflow = "Raw filtered")
bp300     <- read_csv("300_unfiltered.csv") %>% mutate(Workflow = "300 bp unfiltered")
bp300_kma <- read_csv("300_filtered.csv")   %>% mutate(Workflow = "300 bp filtered")

## 2.Combine into one long-format table 

workflow_summary <- bind_rows(raw, raw_kma, bp300, bp300_kma) %>%
  mutate(Workflow = factor(Workflow,
                           levels = c("Raw unfiltered", "Raw filtered",
                                      "300 bp unfiltered", "300 bp filtered")))

## 3.Sanity check - every sample should appear 4 times

workflow_summary %>% count(Sample) %>% filter(n != 4)
# If this returns any rows, those samples are missing from one or more workflow
# CSVs - check before running paired statistics, as they require complete cases

## 4.Friedman test (non-parametric repeated-measures ANOVA equivalent)
##    Appropriate here since: (a) same samples measured under 4 conditions,
##    (b) count data is unlikely to be normally distributed

friedman_richness <- workflow_summary %>% friedman_test(Gene_richness ~ Workflow | Sample)
friedman_count    <- workflow_summary %>% friedman_test(Total_count   ~ Workflow | Sample)

friedman_richness
friedman_count

## 5.Post-hoc pairwise comparisons (if Friedman is significant)
##    Wilcoxon signed-rank test for paired data, with Bonferroni correction
##    for multiple comparisons across the 6 possible workflow pairs

posthoc_richness <- workflow_summary %>%
  wilcox_test(Gene_richness ~ Workflow, paired = TRUE, p.adjust.method = "bonferroni")

posthoc_count <- workflow_summary %>%
  wilcox_test(Total_count ~ Workflow, paired = TRUE, p.adjust.method = "bonferroni")

posthoc_richness
posthoc_count

## 6.Panel A — Richness

p1 <- ggplot(workflow_summary, aes(x = Workflow, y = Gene_richness)) +
  geom_boxplot(aes(fill = Workflow), alpha = 0.25, outlier.shape = NA) +
  geom_point(aes(colour = Workflow), position = position_jitter(width = 0.08),
             size = 2.5, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set2") +
  labs(x = NULL, y = "ARG richness") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14))

## 7.Panel B — Total count

p2 <- ggplot(workflow_summary, aes(x = Workflow, y = Total_count)) +
  geom_boxplot(aes(fill = Workflow), alpha = 0.25, outlier.shape = NA) +
  geom_point(aes(colour = Workflow), position = position_jitter(width = 0.08),
             size = 2.5, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Set2") +
  labs(x = NULL, y = "Total ARG hits") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
        axis.title.y = element_text(size = 14))

## 8.Combine

p1 + p2 + plot_annotation(tag_levels = "A")


## -----------------------------------------------------------------------------

## Figure 3.7 - Sequencing depth bias


## 1. Build per-sample diversity table

alpha_div <- kma_data %>%
  group_by(Sample) %>%
  summarise(
    Richness        = n_distinct(Gene),
    Total_ARG_reads = sum(Depth * Template_length),
    .groups = "drop"
  ) %>%
  left_join(sample_metadata, by = "Sample") %>%
  mutate(
    log_reads     = log10(Total_reads),
    log_ARG_reads = log10(Total_ARG_reads)
  )

## 2. Spearman label helper

spearman_label <- function(x, y) {
  res   <- cor.test(x, y, method = "spearman", exact = FALSE)
  p_str <- if (res$p.value < 0.001) "< 0.001" else sprintf("= %.3f", res$p.value)
  sprintf("Spearman rho = %.2f, p %s", res$estimate, p_str)
}

## 3. Panel A — ARG richness vs sequencing depth

f2a <- ggplot(alpha_div, aes(x = log_reads, y = Richness)) +
  geom_point(aes(colour = Species, shape = Species), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = "lm", se = TRUE,
    colour = "grey30", fill = "grey80", linewidth = 1
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = spearman_label(alpha_div$log_reads, alpha_div$Richness),
    hjust = -0.05, vjust = 1.5, size = 3.2, colour = "grey30"
  ) +
  scale_colour_manual(values = palette_species) +
  scale_shape_manual(values = c(Badger = 16, Fox = 17)) +
  labs(
    x = expression(Log[10]~"total reads"),
    y = "ARG richness (no. unique genes detected)"
  )


## 4. Panel B — Total ARG-mapped reads vs sequencing depth

f2b <- ggplot(alpha_div, aes(x = log_reads, y = log_ARG_reads)) +
  geom_point(aes(colour = Species, shape = Species), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = "lm", se = TRUE,
    colour = "grey30", fill = "grey80", linewidth = 1
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = spearman_label(alpha_div$log_reads, alpha_div$log_ARG_reads),
    hjust = -0.05, vjust = 1.5, size = 3.2, colour = "grey30"
  ) +
  scale_colour_manual(values = palette_species) +
  scale_shape_manual(values = c(Badger = 16, Fox = 17)) +
  labs(
    x = expression(Log[10]~"total reads"),
    y = expression(Log[10]~"total ARG-mapped reads")
  )

## 5. Combine

f2a + f2b +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")


## -----------------------------------------------------------------------------



## Figure 3.8 - HEATMAP resistants

# ---- Load data ----
args_kma_resis <- read_csv("resis_kma_output.csv")

# ---- Get all sample names including blanks ----
all_samples <- args_kma_resis %>%
  distinct(Sample) %>%
  pull(Sample)

# ---- Get gene frequency across resistant samples with ARGs only ----
gene_freq <- args_kma_resis %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Sample, Gene) %>%
  count(Gene, name = "n_samples") %>%
  arrange(desc(n_samples))

# ---- Take top 30 most commonly occurring genes ----
top30_genes <- gene_freq %>%
  slice_max(n_samples, n = 30) %>%
  pull(Gene)

# ---- Build presence/absence matrix for samples WITH ARGs ----
args_resis_mat <- args_kma_resis %>%
  filter(!is.na(Gene), Gene != "") %>%
  filter(Gene %in% top30_genes) %>%
  distinct(Sample, Gene) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Gene,
              values_from = present,
              values_fill = 0) %>%
  column_to_rownames("Sample")

# ---- Add zero rows for samples with no ARGs ----
missing_samples <- setdiff(all_samples, rownames(args_resis_mat))

if (length(missing_samples) > 0) {
  zero_rows <- matrix(0,
                      nrow = length(missing_samples),
                      ncol = ncol(args_resis_mat),
                      dimnames = list(missing_samples, colnames(args_resis_mat)))
  args_resis_mat <- rbind(args_resis_mat, zero_rows)
}

# ---- Reorder columns alphabetically ----
args_resis_mat <- args_resis_mat[, sort(colnames(args_resis_mat))]

# ---- Transpose so genes are rows and samples are columns ----
args_resis_mat_t <- t(args_resis_mat)

# ---- Create column annotation ----
col_annotation <- data.frame(
  ARG_status = ifelse(colnames(args_resis_mat_t) %in% missing_samples,
                      "No ARGs detected", "ARGs detected"),
  row.names = colnames(args_resis_mat_t)
)

annotation_colours <- list(
  ARG_status = c("No ARGs detected" = "#999999",
                 "ARGs detected" = "black")
)

# ---- Cluster ARG-positive samples only ----
positive_samples <- setdiff(colnames(args_resis_mat_t), missing_samples)
positive_mat <- args_resis_mat_t[, positive_samples]

# euclidean distance clustering
col_clust <- hclust(dist(t(positive_mat), method = "euclidean"),
                    method = "complete")

# get clustered order of positive samples
clustered_order <- positive_samples[col_clust$order]

# ---- Combine: clustered positive samples + zero samples at end ----
final_col_order <- c(clustered_order, missing_samples)

# reorder matrix and annotation
args_resis_mat_t_ordered <- args_resis_mat_t[, final_col_order]
col_annotation_ordered <- col_annotation[final_col_order, , drop = FALSE]


# ---- Plot heatmap ----
args_resis_heatmap <- pheatmap(args_resis_mat_t_ordered,
                               color = c("white", "#E07B6A"),
                               annotation_col = col_annotation_ordered,
                               annotation_colors = annotation_colours,
                               legend = FALSE,
                               cluster_rows = FALSE,
                               cluster_cols = FALSE,
                               border_color = "grey80",
                               fontsize_row = 12,
                               fontsize_col = 12)



## Figure 3.8 HEATMAP Susceptibles 

args_kma_susc <- read_csv("susc_kma_output.csv")


# Get all sample names including blanks
all_susc_samples <- args_kma_susc %>%
  distinct(Sample) %>%
  pull(Sample)

# Get gene frequency across susceptible samples with ARGs only
gene_freq_susc <- args_kma_susc %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Sample, Gene) %>%
  count(Gene, name = "n_samples") %>%
  arrange(desc(n_samples))

# Take top 30 most commonly occurring genes
top30_genes_susc <- gene_freq_susc %>%
  slice_max(n_samples, n = 30) %>%
  pull(Gene)

# Build presence/absence matrix for samples WITH ARGs
args_susc_mat <- args_kma_susc %>%
  filter(!is.na(Gene), Gene != "") %>%
  filter(Gene %in% top30_genes) %>%
  distinct(Sample, Gene) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Gene,
              values_from = present,
              values_fill = 0) %>%
  column_to_rownames("Sample")

# Add zero rows for samples with no ARGs
missing_samples <- setdiff(all_susc_samples, rownames(args_susc_mat))

if (length(missing_samples) > 0) {
  zero_rows <- matrix(0,
                      nrow = length(missing_samples),
                      ncol = ncol(args_susc_mat),
                      dimnames = list(missing_samples, colnames(args_susc_mat)))
  args_susc_mat <- rbind(args_susc_mat, zero_rows)
}

# Reorder columns alphabetically
args_susc_mat <- args_susc_mat[, sort(colnames(args_susc_mat))]

# Transpose so genes are rows and samples are columns
args_susc_mat_t <- t(args_susc_mat)

# Create column annotation
col_annotation <- data.frame(
  ARG_status = ifelse(colnames(args_susc_mat_t) %in% missing_samples,
                      "No ARGs detected", "ARGs detected"),
  row.names = colnames(args_susc_mat_t)
)

annotation_colours <- list(
  ARG_status = c("No ARGs detected" = "#999999",
                 "ARGs detected" = "black")
)

# ---- Cluster ARG-positive samples only ----
positive_samples <- setdiff(colnames(args_susc_mat_t), missing_samples)
positive_mat <- args_susc_mat_t[, positive_samples]

# euclidean distance clustering
col_clust <- hclust(dist(t(positive_mat), method = "euclidean"),
                    method = "complete")

# get clustered order of positive samples
clustered_order <- positive_samples[col_clust$order]

# ---- Combine: clustered positive samples + zero samples at end ----
final_col_order <- c(clustered_order, missing_samples)

# reorder matrix and annotation
args_susc_mat_t_ordered <- args_susc_mat_t[, final_col_order]
col_annotation_ordered <- col_annotation[final_col_order, , drop = FALSE]


# ---- Plot heatmap ----
args_susc_heatmap <- pheatmap(args_susc_mat_t_ordered,
                              color = c("white", "steelblue"),
                              annotation_col = col_annotation_ordered,
                              annotation_colors = annotation_colours,
                              legend = FALSE,
                              cluster_rows = FALSE,
                              cluster_cols = FALSE,
                              border_color = "grey80",
                              fontsize_row = 12,
                              fontsize_col = 12)

## Figure 3.9a and b


# ---- Load data ----
df_resistant <- read.csv("ARGs_resistant_class.csv", stringsAsFactors = FALSE)
df_susceptible <- read.csv("ARGs_susceptible_class.csv", stringsAsFactors = FALSE)

# ---- Clean class names ----
df_resistant <- df_resistant %>% mutate(Class = str_trim(Class))
df_susceptible <- df_susceptible %>% mutate(Class = str_trim(Class))

# ---- Collapse to distinct gene families per sample per class ----
df_resistant_collapsed <- df_resistant %>%
  distinct(Sample, Gene, Class) %>%
  group_by(Sample, Class) %>%
  summarise(Count = n_distinct(Gene), .groups = "drop")

df_susceptible_collapsed <- df_susceptible %>%
  distinct(Sample, Gene, Class) %>%
  group_by(Sample, Class) %>%
  summarise(Count = n_distinct(Gene), .groups = "drop")

# ---- Build consistent class levels from raw data ----
all_classes <- union(unique(df_resistant$Class),
                     unique(df_susceptible$Class)) %>%
  na.omit() %>%
  str_trim()

# ---- Add missing classes to resistant collapsed data with zero counts ----
# This forces all classes into p1 legend without plotting anything visible
missing_classes <- setdiff(all_classes,
                           unique(df_resistant_collapsed$Class))

# get real sample names for resistant
resistant_samples <- unique(df_resistant_collapsed$Sample)

if (length(missing_classes) > 0) {
  dummy_rows <- expand.grid(
    Sample = resistant_samples[1],
    Class = missing_classes,
    Count = 0,
    stringsAsFactors = FALSE
  )
  df_resistant_collapsed <- bind_rows(df_resistant_collapsed, dummy_rows)
}

# ---- Set factor levels ----
df_resistant_collapsed$Class <- factor(df_resistant_collapsed$Class,
                                       levels = all_classes)
df_susceptible_collapsed$Class <- factor(df_susceptible_collapsed$Class,
                                         levels = all_classes)

# ---- Define colour palette ----
class_colours <- c(
  "Tetracyclines" = "#8DD3C7",
  "Aminoglycosides" = "#FFFFB3",
  "Lincosamides" = "#BEBADA",
  "Phenicols" = "#FB8072",
  "Sulfonamides" = "#80B1D3",
  "Trimethoprim" = "#FDB462",
  "beta-lactams" = "#B3DE69",
  "MLSB" = "#FCCDE5",
  "Macrolides" = "#D9D9D9",
  "Oxazolidinones" = "#BC80BD",
  "MDR" = "#CCEBC5",
  "Glycopeptides" = "#FFED6F",
  "Fosfomycin" = "#E41A1C",
  "Ionophore resistance" = "#377EB8",
  "Streptogramins (A)" = "#4DAF4A",
  "Fluoroquinolones/Quinolones" = "gold"
)

# ---- Plot 1: Resistant ----
p1 <- ggplot(df_resistant_collapsed %>% filter(Sample != "DUMMY"),
             aes(x = factor(Sample), y = Count, fill = Class)) +
  geom_bar(stat = "identity") +
  # invisible dummy layer to force all classes into legend
  geom_bar(data = df_resistant_collapsed %>% filter(Count == 0),
           aes(x = factor(Sample), y = Count, fill = Class),
           stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = class_colours, drop = FALSE) +
  labs(
    x = "Sample",
    y = "Number of distinct ARG families",
    fill = "Antibiotic Class",
    title = "Phenotypically Resistant"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

# ---- Plot 2: Susceptible ----
p2 <- ggplot(df_susceptible_collapsed, aes(x = factor(Sample), y = Count, fill = Class)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = class_colours, drop = FALSE) +
  labs(
    x = "Sample",
    y = "Number of distinct ARG families",
    fill = "Antibiotic Class",
    title = "Phenotypically Susceptible"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )

# ---- Combine with patchwork ----
ARGfamilydiversity_byclass <- p1 / p2 +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "ARG family diversity by antibiotic class",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
  )

ARGfamilydiversity_byclass


## Figure 3.9 c and d

df <- read.csv("ARGs_res_susc.csv", stringsAsFactors = FALSE)
# ---- Define group sample sizes ----
n_resistant <- 30
n_susceptible <- 16
n_badger <- 35
n_fox <- 11


# ---- Phenotype counts ----
df_counts_char <- df %>%
  distinct(Sample, Gene, Characteristic, Class) %>%
  group_by(Characteristic, Class) %>%
  summarise(n_samples = n_distinct(Sample), .groups = "drop") %>%
  mutate(Percentage = case_when(
    Characteristic == "Resistant" ~ (n_samples / n_resistant) * 100,
    Characteristic == "Susceptible" ~ (n_samples / n_susceptible) * 100
  ))

# Species counts
df_counts_species <- df %>%
  distinct(Sample, Gene, Species, Class) %>%
  group_by(Species, Class) %>%
  summarise(n_samples = n_distinct(Sample), .groups = "drop") %>%
  mutate(Percentage = case_when(
    Species == "Badger" ~ (n_samples / n_badger) * 100,
    Species == "Fox" ~ (n_samples / n_fox) * 100
  ))


# Plot phenotype
px <- ggplot(df_counts_char, aes(x = Class, y = Percentage, fill = Characteristic)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  theme_minimal() +
  scale_fill_manual(values = c("Resistant" = "#E07B6A", "Susceptible" = "#6ABCE0")) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 20)) +  # gridlines at 0,20,40,60,80,100
  labs(
    x = "Antibiotic Class",
    y = "Percentage of samples (%)",
    fill = "Phenotype",
    title = "ARG family diversity by antibiotic class and phenotype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )


# Plot species
py <- ggplot(df_counts_species, aes(x = Class, y = Percentage, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  theme_minimal() +
  scale_fill_manual(values = c("Badger" = "#4C6E5D", "Fox" = "#C97B3D")) +
  labs(
    x = "Antibiotic Class",
    y = "Percentage of samples (%)",
    fill = "Species",
    title = "ARG family diversity by antibiotic class and host species"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

# Combine with patchwork
ARGfamilydiversity_pheno_host <- px / py +
  plot_annotation(
    title = "ARG family diversity by antibiotic class",
    tag_levels = list(c("C", "D")),
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
  )

print(ARGfamilydiversity_pheno_host)


## -----------------------------------------------------------------------------


# Sequencing depth checks


metadata_total <- read_csv("Metadata_total.csv")

depth_check <- full_metadata %>%
  left_join(metadata_total %>% select(Sample, Total_reads), by = "Sample")

sum(is.na(depth_check$Total_reads))   # should be 0

depth_check %>%
  group_by(Characteristic) %>%
  summarise(
    n            = n(),
    min_reads    = min(Total_reads),
    max_reads    = max(Total_reads),
    mean_reads   = mean(Total_reads),
    median_reads = median(Total_reads)
  )

wilcox.test(Total_reads ~ Characteristic, data = depth_check)

# ---- Load data ----
metadata <- read.csv("Metadata_total.csv", stringsAsFactors = FALSE)
df <- read.csv("kma_all_samples.csv", stringsAsFactors = FALSE)

## 1. Check sequencing depths across cohorts

metadata %>%
  group_by(Characteristic) %>%
  summarise(n = n(), median_reads = median(Total_reads),
            min_reads = min(Total_reads), max_reads = max(Total_reads))

wilcox.test(Total_reads ~ Characteristic, data = metadata)

metadata %>%
  group_by(Species) %>%
  summarise(n = n(), median_reads = median(Total_reads),
            min_reads = min(Total_reads), max_reads = max(Total_reads))

wilcox.test(Total_reads ~ Species, data = metadata)

# Diversity metrics calculations

df_clean <- df %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::group_by(Sample, Gene, Characteristic) %>%
  dplyr::summarise(Depth = sum(Depth, na.rm = TRUE),
                   Template_length = mean(Template_length, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::mutate(abundance = Depth * Template_length)

abund_matrix <- df_clean %>%
  dplyr::select(Sample, Gene, abundance) %>%
  pivot_wider(names_from = Gene, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("Sample")

diversity_df <- data.frame(
  Sample = rownames(abund_matrix),
  Richness = rowSums(abund_matrix > 0),
  Shannon = diversity(abund_matrix, index = "shannon"),
  Simpson = diversity(abund_matrix, index = "simpson")
)

diversity_df <- diversity_df %>%
  mutate(Sample = as.character(Sample)) %>%
  left_join(metadata %>% mutate(Sample = as.character(Sample)), by = "Sample")

zero_samples <- metadata %>%
  filter(!Sample %in% diversity_df$Sample) %>%
  mutate(Sample = as.character(Sample), Richness = 0, Shannon = 0, Simpson = 0)

diversity_df <- bind_rows(diversity_df, zero_samples)

dim(diversity_df)

diversity_df %>%
  group_by(Characteristic) %>%
  summarise(n = n(), median_Shannon = median(Shannon),
            median_Simpson = median(Simpson), median_Richness = median(Richness))

diversity_df %>%
  group_by(Species) %>%
  summarise(n = n(), median_Shannon = median(Shannon),
            median_Simpson = median(Simpson), median_Richness = median(Richness))



# 1. Reusable boxplot builder 

make_diversity_boxplot <- function(data, x_var, y_var, fill_values,
                                   x_label, y_label, title) {
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2) +
    theme_bw() +
    scale_fill_manual(values = fill_values) +
    labs(x = x_label, y = y_label, title = title) +
    theme(legend.position = "none")
}

# 2. Reusable depth-vs-metric scatter plot builder 
make_depth_scatter <- function(data, colour_var, colour_values, title) {
  ggplot(data, aes(x = log10(Total_reads), y = ARG_count,
                   colour = .data[[colour_var]])) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_bw() +
    scale_colour_manual(values = colour_values) +
    labs(x = "log10(Total reads)", y = "ARG count", title = title)
}

# ---- 3. Reusable additive + interaction lm() pair ----

run_lm_pair <- function(response, group_var, data) {
  f_add <- as.formula(paste(response, "~ log10(Total_reads) +", group_var))
  f_int <- as.formula(paste(response, "~ log10(Total_reads) *", group_var))
  list(
    additive    = lm(f_add, data = data),
    interaction = lm(f_int, data = data)
  )
}

# Colour palettes used throughout - defined once
char_colours <- c("Resistant" = "#E07B6A", "Susceptible" = "#6ABCE0")
sp_colours   <- c("Badger" = "#4C6E5D", "Fox" = "#C97B3D")


# Depth vs ARG count correlation - Figure 3.10


cor.test(log10(metadata$Total_reads), metadata$ARG_count, method = "spearman")

depth_ARG_char <- make_depth_scatter(metadata, "Characteristic", char_colours,
                                     "Sequencing depth vs ARG count by phenotype")

depth_ARG_species <- make_depth_scatter(metadata, "Species", sp_colours,
                                        "Sequencing depth vs ARG count by host species")

depth_plots <- depth_ARG_char + depth_ARG_species + plot_annotation(tag_levels = "A")
depth_plots


# Linear models - phenotype and host species

# Phenotype
models_richness_char <- run_lm_pair("ARG_count", "Characteristic", metadata)
summary(models_richness_char$additive)
summary(models_richness_char$interaction)

models_shannon_char <- run_lm_pair("Shannon", "Characteristic", diversity_df)
summary(models_shannon_char$additive)
summary(models_shannon_char$interaction)

models_simpson_char <- run_lm_pair("Simpson", "Characteristic", diversity_df)
summary(models_simpson_char$additive)
summary(models_simpson_char$interaction)

# ---- Host species ----
models_richness_sp <- run_lm_pair("ARG_count", "Species", metadata)
summary(models_richness_sp$additive)
summary(models_richness_sp$interaction)

models_shannon_sp <- run_lm_pair("Shannon", "Species", diversity_df)
summary(models_shannon_sp$additive)
summary(models_shannon_sp$interaction)

models_simpson_sp <- run_lm_pair("Simpson", "Species", diversity_df)
summary(models_simpson_sp$additive)
summary(models_simpson_sp$interaction)

## Figures 3.11 and 3.12

# Phenotype boxplots
p_richness_char <- make_diversity_boxplot(diversity_df, "Characteristic", "Richness",
                                          char_colours, "Phenotype", "ARG richness", "ARG richness")
p_shannon_char  <- make_diversity_boxplot(diversity_df, "Characteristic", "Shannon",
                                          char_colours, "Phenotype", "Shannon diversity", "Shannon diversity")
p_simpson_char  <- make_diversity_boxplot(diversity_df, "Characteristic", "Simpson",
                                          char_colours, "Phenotype", "Simpson diversity", "Simpson diversity")

# Species boxplots
p_richness_sp <- make_diversity_boxplot(diversity_df, "Species", "Richness",
                                        sp_colours, "Host species", "ARG richness", "ARG richness")
p_shannon_sp  <- make_diversity_boxplot(diversity_df, "Species", "Shannon",
                                        sp_colours, "Host species", "Shannon diversity", "Shannon diversity")
p_simpson_sp  <- make_diversity_boxplot(diversity_df, "Species", "Simpson",
                                        sp_colours, "Host species", "Simpson diversity", "Simpson diversity")

# Combine
combined_char <- p_richness_char + p_shannon_char + p_simpson_char +
  plot_annotation(
    title = "ARG diversity metrics by phenotype",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
  )

combined_sp <- p_richness_sp + p_shannon_sp + p_simpson_sp +
  plot_annotation(
    title = "ARG diversity metrics by host species",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
  )

combined_char
combined_sp

## -----------------------------------------------------------------------------


# Step 1: Load data
df <- read.csv("kma_all_samples.csv", stringsAsFactors = FALSE)
metadata <- read.csv("Metadata_total.csv", stringsAsFactors = FALSE)

# Step 2: Clean data
df_clean <- df %>%
  filter(!is.na(Gene), Gene != "") %>%
  group_by(Sample, Gene, Characteristic) %>%
  summarise(Depth = sum(Depth, na.rm = TRUE),
            Template_length = mean(Template_length, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(abundance = Depth * Template_length)

# ---- Step 3: Build presence/absence matrix ----
pa_matrix <- df_clean %>%
  distinct(Sample, Gene) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Gene,
              values_from = present,
              values_fill = 0)

# add zero rows for the 11 empty samples
all_samples <- df %>% distinct(Sample, Characteristic)
missing_samples <- all_samples %>%
  filter(!Sample %in% pa_matrix$Sample) %>%
  select(Sample)

pa_matrix <- bind_rows(pa_matrix,
                       missing_samples %>%
                         mutate(across(everything(), ~replace_na(., 0)))) %>%
  replace(is.na(.), 0)

# ---- Step 4: Join metadata ----
meta <- metadata %>%
  select(Sample, Species, Characteristic, Total_reads) %>%
  mutate(Sample = as.character(Sample))

pa_matrix <- pa_matrix %>%
  mutate(Sample = as.character(Sample))

pa_joined <- pa_matrix %>%
  left_join(meta, by = "Sample")

# Step 5: Extract matrix and metadata
meta_ordered <- pa_joined %>%
  select(Sample, Species, Characteristic, Total_reads)

pa_mat <- pa_joined %>%
  select(-Sample, -Species, -Characteristic, -Total_reads) %>%
  as.matrix()

rownames(pa_mat) <- meta_ordered$Sample

# Step 6: Remove zero samples for beta diversity
zero_samples <- rowSums(pa_mat) == 0
pa_mat_beta <- pa_mat[!zero_samples, ]
meta_beta <- meta_ordered[!zero_samples, ]

# Sanity checks
dim(pa_mat_beta)        # should be 46 x 165
nrow(meta_beta)         # should be 46
table(meta_beta$Species)
table(meta_beta$Characteristic)

# Step 7: Jaccard distance matrix
jaccard_dist <- vegdist(pa_mat_beta, method = "jaccard", binary = TRUE)

# Step 8: NMDS
set.seed(123)
nmds <- metaMDS(jaccard_dist, k = 2, trymax = 100)
nmds$stress

# ---- Step 9: Extract NMDS scores ----
nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  left_join(meta_beta, by = "Sample")

# ---- Step 10: Figure 3.13 - NMDS plot

library(patchwork)

# Plot by phenotype
p_char <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2,
                                  colour = Characteristic,
                                  shape = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Characteristic),
               type = "t", level = 0.95, linewidth = 0.8) +
  theme_bw() +
  labs(subtitle = paste0("Stress = ", round(nmds$stress, 3)),
       colour = "Phenotype",
       shape = "Host species") +
  scale_colour_manual(values = c("Resistant" = "#E07B6A",
                                 "Susceptible" = "#6ABCE0")) +
  theme(legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))

# Plot by species
p_species <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2,
                                     colour = Species,
                                     shape = Characteristic)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Species),
               type = "t", level = 0.95, linewidth = 0.8) +
  theme_bw() +
  labs(subtitle = paste0("Stress = ", round(nmds$stress, 3)),
       colour = "Host species",
       shape = "Phenotype") +
  scale_colour_manual(values = c("Badger" = "#4C6E5D",
                                 "Fox" = "#C97B3D")) +
  theme(legend.title = element_text(size = 11),    # this was missing
        legend.text = element_text(size = 10))      # closing out p_species

# Combine
p_char + p_species +
  plot_annotation(
    title = "NMDS of ARG profiles (Jaccard distance)",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
  )
 theme(plot.title = element_text(size = 14, hjust = 0.5))
  
 
print(nmds_combined)



# Step 11a: PERMANOVA
set.seed(123)

permanova_char <- adonis2(jaccard_dist ~ Characteristic,
                          data = meta_beta,
                          permutations = 999)

permanova_species <- adonis2(jaccard_dist ~ Species,
                             data = meta_beta,
                             permutations = 999)

permanova_both <- adonis2(jaccard_dist ~ Characteristic + Species,
                          data = meta_beta,
                          permutations = 999)

permanova_char
permanova_species
permanova_both

permanova_both <- adonis2(jaccard_dist ~ Characteristic + Species,
                          data = meta_beta,
                          permutations = 999,
                          by = "term")  # tests each term individually

permanova_both


# 11b PERMANOVA with sequencing depth as covariate
set.seed(123)

# depth + phenotype + species together
permanova_depth_both <- adonis2(jaccard_dist ~ log10(Total_reads) + Characteristic + Species,
                                data = meta_beta,
                                permutations = 999,
                                by = "term")

permanova_depth_both

# depth + phenotype only
permanova_depth_char <- adonis2(jaccard_dist ~ log10(Total_reads) + Characteristic,
                                data = meta_beta,
                                permutations = 999,
                                by = "term")

permanova_depth_char

# depth + species only
permanova_depth_species <- adonis2(jaccard_dist ~ log10(Total_reads) + Species,
                                   data = meta_beta,
                                   permutations = 999,
                                   by = "term")

permanova_depth_species


# ---- Step 12: betadisper ----

disp_char <- betadisper(jaccard_dist, meta_beta$Characteristic)
permutest_char <- permutest(disp_char, permutations = 999)
print(permutest_char)
print(disp_char)

disp_species <- betadisper(jaccard_dist, meta_beta$Species)
permutest_species <- permutest(disp_species, permutations = 999)
print(permutest_species)
print(disp_species)

# ---- Step 13: Dispersion vs sequencing depth ----
disp_df <- data.frame(
  Sample = meta_beta$Sample,
  Characteristic = meta_beta$Characteristic,
  Species = meta_beta$Species,
  Total_reads = meta_beta$Total_reads,
  Distance_to_centroid_char = disp_char$distances,
  Distance_to_centroid_species = disp_species$distances
)

# ---- Correlation tests ----
cor_char <- cor.test(disp_df$Distance_to_centroid_char,
                     log10(disp_df$Total_reads),
                     method = "spearman")
print(cor_char)

cor_species <- cor.test(disp_df$Distance_to_centroid_species,
                        log10(disp_df$Total_reads),
                        method = "spearman")
print(cor_species)




## Figure 3.14
# Plot by phenotype
disperse_plot1 <- ggplot(disp_df, aes(x = log10(Total_reads),
                                      y = Distance_to_centroid_char,
                                      colour = Characteristic)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  scale_colour_manual(values = c("Resistant" = "#E07B6A",
                                 "Susceptible" = "#6ABCE0")) +
  labs(x = "log10(Total reads)",
       y = "Distance to group centroid",
       title = "By phenotype")

# ---- Plot by species ----
disperse_plot2 <- ggplot(disp_df, aes(x = log10(Total_reads),
                                      y = Distance_to_centroid_species,
                                      colour = Species)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  scale_colour_manual(values = c("Badger" = "#4C6E5D",
                                 "Fox" = "#C97B3D")) +
  labs(x = "log10(Total reads)",
       y = "Distance to group centroid",
       title = "By species")

# Combine and print
dispersion_combined <- disperse_plot1 + disperse_plot2 +
  plot_annotation(
    title = "Dispersion versus sequencing depth",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
  )

print(dispersion_combined)

## -----------------------------------------------------------------------------

# Step 14: Wilcoxon test - ARG richness by species - Figure 3.15
richness_df <- data.frame(
  Sample = rownames(pa_mat_beta),
  Richness = rowSums(pa_mat_beta)) %>%
  left_join(meta_beta, by = "Sample")

wilcox.test(Richness ~ Species, data = richness_df)

ARGrichnessbyspecies <- ggplot(richness_df, aes(x = Species, y = Richness, fill = Species)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 4) +
  theme_bw() +
  labs(x = "Host species",
       y = "ARG richness (number of distinct genes)",
       title = "ARG richness by host species") +
  scale_fill_manual(values = c("Badger" = "#4C6E5D", "Fox" = "#C97B3D")) +
  theme(legend.position = "none",
        axis.text = element_text(size = 24),   # tick label size
        axis.title = element_text(size = 24))  # axis title size


# check actual mean dispersions per group for each betadisper
print(disp_char)
print(disp_species)

## -----------------------------------------------------------------------------

## ARGs in Assembled Contigs: Heatmap, Summary Stats, and Assembly-Quality Check
## Badger vs Fox

# READ AND CLEAN DATA

df <- read_csv("ARG_FILE.csv", show_col_types = FALSE) %>%
  mutate(Sample = as.character(Sample)) %>%   # Sample IDs look numeric (e.g. 86002),
  # so read_csv imports as double by default.
  # column_to_rownames() later forces character
  # anyway, so fixing this once here avoids a
  # silent type mismatch in joins downstream.
  filter(!is.na(Sample)) %>%                   # remove any completely blank rows
  filter(!Sample %in% c("155007", "227016"))   # exclude the 2 failed-assembly samples

# Quick check: each Sample should map to exactly one Species
dup_species_check <- df %>% distinct(Sample, Species) %>% count(Sample) %>% filter(n > 1)
if (nrow(dup_species_check) > 0) {
  warning("Some samples have more than one Species label - check these before continuing:")
  print(dup_species_check)
}


message("Sample counts after cleaning:")
df %>% distinct(Sample, Species) %>% count(Species) %>% print()


# 2. HEATMAP: ARG hit counts, Badger vs Fox

# ---- 2a. Select top genes ----
# Exclude blank/NA Gene rows (= zero-hit samples) before counting, so a
# zero-hit sample's blank row can't accidentally be counted as a "gene"
top_genes <- df %>%
  filter(!is.na(Gene)) %>%
  count(Gene, sort = TRUE) %>%
  slice_max(n, n = 90) %>%
  pull(Gene)

df_top <- df %>% filter(Gene %in% top_genes)

# ---- 2b. Build Sample x Gene hit-count matrix ----
count_df <- df_top %>% count(Sample, Gene, name = "n")

mat <- count_df %>%
  pivot_wider(names_from = Gene, values_from = n, values_fill = 0) %>%
  column_to_rownames("Sample") %>%
  as.matrix()

# ---- 2c. Build sample metadata (Species) ----
# distinct() uses the FULL df (not df_top), so zero-hit samples - which
# still have a blank/NA Gene row but a real Sample/Species value - are
# captured too
sample_meta <- df %>%
  distinct(Sample, Species) %>%
  column_to_rownames("Sample")

# ---- 2d. Add back any samples missing from mat entirely ----
# "Missing from mat" can mean two different things:
#   (a) genuinely zero ARG hits anywhere in the full dataset, or
#   (b) the sample DOES have ARG hits, but none happen to be among the
#       top 90 genes selected above
# These need to be told apart, or (b)-type samples get mislabelled as
# having no ARGs at all.

full_samples <- rownames(sample_meta)
missing_samples <- setdiff(full_samples, rownames(mat))

sample_total_hits <- df %>%
  filter(!is.na(Gene)) %>%
  count(Sample, name = "total_hits_all_genes")

hit_check <- tibble(Sample = missing_samples) %>%
  left_join(sample_total_hits, by = "Sample") %>%
  mutate(total_hits_all_genes = replace_na(total_hits_all_genes, 0))

genuinely_zero <- hit_check %>% filter(total_hits_all_genes == 0) %>% pull(Sample)
excluded_by_filter <- hit_check %>% filter(total_hits_all_genes > 0) %>% pull(Sample)

if (length(genuinely_zero) > 0) {
  message(length(genuinely_zero), " sample(s) have NO ARG hits at all: ",
          paste(genuinely_zero, collapse = ", "))
}

if (length(excluded_by_filter) > 0) {
  message(length(excluded_by_filter),
          " sample(s) DO have ARG hits, but none among the top 90 genes shown - ",
          "these will still appear as blank rows below, which could be misleading: ",
          paste(excluded_by_filter, collapse = ", "))
}

if (length(missing_samples) > 0) {
  zero_mat <- matrix(0,
                     nrow = length(missing_samples),
                     ncol = ncol(mat),
                     dimnames = list(missing_samples, colnames(mat)))
  mat <- rbind(mat, zero_mat)
}

# ---- 2e. Order rows: Badger block first, then Fox; by total richness within each ----
sample_order <- sample_meta %>%
  rownames_to_column("Sample") %>%
  left_join(
    tibble(Sample = rownames(mat), total = rowSums(mat)),
    by = "Sample"
  ) %>%
  arrange(Species, desc(total)) %>%
  pull(Sample)

mat <- mat[sample_order, ]
sample_meta <- sample_meta[sample_order, , drop = FALSE]

gap_after <- sum(sample_meta$Species == "Badger")

# ---- 2f. Drop genuinely zero-hit samples for THIS figure only ----
# (purely visual decluttering - doesn't affect any stats elsewhere)
n_before <- nrow(mat)
mat <- mat[rowSums(mat) > 0, , drop = FALSE]
sample_meta <- sample_meta[rownames(mat), , drop = FALSE]
n_after <- nrow(mat)

message(n_before - n_after, " zero-hit sample(s) dropped from this figure. ",
        "N shown = ", n_after, " of ", n_before, " total. ",
        "Report this reduced N in the figure caption.")

gap_after <- sum(sample_meta$Species == "Badger")  # recompute after dropping rows

# ---- 2g. Transpose: Samples = columns (x-axis), Genes = rows (y-axis) ----
mat_t <- t(mat)

# ---- 2h. Color scale: one distinct colour per integer count value ----
max_val <- max(mat_t)
break_points <- seq(-0.5, max_val + 0.5, by = 1)
n_colors <- length(break_points) - 1

orange_red_palette <- c("grey92",
                        colorRampPalette(brewer.pal(9, "YlOrRd"))(n_colors - 1))

annotation_colors <- list(
  Species = c(Badger = "#4C6E5D", Fox = "#C97B3D")
)

# ---- 2i. Heatmap - Figure 3.17
pheatmap(
  mat_t,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = sample_meta,
  annotation_colors = annotation_colors,
  gaps_col = gap_after,
  color = orange_red_palette,
  breaks = break_points,
  legend_breaks = 0:max_val,
  legend_labels = as.character(0:max_val),
  border_color = "grey85",
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "ARG hit counts on assembled contigs (Badger vs Fox)",
  treeheight_row = 40,
  legend = TRUE
)


# TOTAL UNIQUE GENES DETECTED (full dataset, not just top 90)

n_unique_genes <- df %>%
  filter(!is.na(Gene)) %>%
  distinct(Gene) %>%
  nrow()

message("Total unique genes detected (CARD, all samples): ", n_unique_genes)


# DOES ARG DETECTION CORRELATE WITH ASSEMBLY QUALITY?

# Uses the SAME cleaned `df` from Section 1 - no second read of the CSV,
# so the failed-assembly exclusion applies here automatically too.

quast <- read_csv("Quast_Results.csv", show_col_types = FALSE) %>%
  mutate(Sample = as.character(Sample)) %>%
  select(Sample, total_length = `Total length`, N50)

df_hits <- df %>% filter(!is.na(Gene))

# ---- 4a. Per-sample ARG richness and total hit count ----
arg_summary <- df_hits %>%
  group_by(Sample, Species) %>%
  summarise(
    arg_richness = n_distinct(Gene),
    arg_total_hits = n(),
    .groups = "drop"
  )

# Add back zero-ARG samples as genuine zeros
full_samples_df <- df %>% distinct(Sample, Species)
arg_summary <- full_samples_df %>%
  left_join(arg_summary, by = c("Sample", "Species")) %>%
  mutate(across(c(arg_richness, arg_total_hits), ~replace_na(., 0)))

diag_df <- arg_summary %>% left_join(quast, by = "Sample")

# Sanity check - no failed-assembly samples, no missing QUAST data
message("diag_df sample count by species:")
diag_df %>% count(Species) %>% print()

missing_quast <- diag_df %>% filter(is.na(total_length) | is.na(N50))
if (nrow(missing_quast) > 0) {
  warning("Some samples are missing QUAST data - check these before trusting correlations:")
  print(missing_quast)
}

# ---- 4b. Spearman correlations: overall
cat("=== OVERALL ===\n")

r1 <- cor.test(diag_df$N50, diag_df$arg_richness, method = "spearman", exact = FALSE)
cat("Richness vs N50: rho =", round(r1$estimate, 3), ", p =", signif(r1$p.value, 3), "\n")

r2 <- cor.test(diag_df$total_length, diag_df$arg_richness, method = "spearman", exact = FALSE)
cat("Richness vs Total length: rho =", round(r2$estimate, 3), ", p =", signif(r2$p.value, 3), "\n")

r3 <- cor.test(diag_df$N50, diag_df$arg_total_hits, method = "spearman", exact = FALSE)
cat("Total hits vs N50: rho =", round(r3$estimate, 3), ", p =", signif(r3$p.value, 3), "\n")

r4 <- cor.test(diag_df$total_length, diag_df$arg_total_hits, method = "spearman", exact = FALSE)
cat("Total hits vs Total length: rho =", round(r4$estimate, 3), ", p =", signif(r4$p.value, 3), "\n\n")

# ---- 4c. Spearman correlations: by species ----
cat("=== BY SPECIES ===\n")

diag_df %>%
  group_by(Species) %>%
  summarise(
    n = n(),
    rho_richness_N50 = cor.test(N50, arg_richness, method = "spearman", exact = FALSE)$estimate,
    p_richness_N50 = cor.test(N50, arg_richness, method = "spearman", exact = FALSE)$p.value,
    rho_richness_length = cor.test(total_length, arg_richness, method = "spearman", exact = FALSE)$estimate,
    p_richness_length = cor.test(total_length, arg_richness, method = "spearman", exact = FALSE)$p.value,
    rho_hits_N50 = cor.test(N50, arg_total_hits, method = "spearman", exact = FALSE)$estimate,
    p_hits_N50 = cor.test(N50, arg_total_hits, method = "spearman", exact = FALSE)$p.value,
    rho_hits_length = cor.test(total_length, arg_total_hits, method = "spearman", exact = FALSE)$estimate,
    p_hits_length = cor.test(total_length, arg_total_hits, method = "spearman", exact = FALSE)$p.value
  ) %>%
  print(width = Inf)




## -----------------------------------------------------------------------------

## Plasmids

# 1. Read and clean PlasmidFinder data

pf <- read_csv("Plasmids_on_Assemblies.csv", show_col_types = FALSE) %>%
  mutate(Sample = as.character(Sample)) %>%
  filter(!is.na(Sample))  # drop any stray blank rows

# NOTE: GENE names are assumed to already be manually cleaned in the CSV
# (allele numbers like "_1"/"_2" removed, e.g. ColRNAI_1 -> ColRNAI),
# so no further name-cleaning is done here.

# Sample -> Species lookup (assumes zero-hit samples still appear as a
# row with Sample/Species filled but GENE blank - verify this holds)
sample_meta <- pf %>%
  distinct(Sample, Species) %>%
  column_to_rownames("Sample")

full_samples <- rownames(sample_meta)

# 2. Per-sample summary metrics

pf_hits <- pf %>% filter(!is.na(GENE))

# Metric 1: total distinct PLASMID CONTIGS per sample
# (collapses multi-replicon contigs like contig_304 down to 1 contig)
contig_counts <- pf_hits %>%
  distinct(Sample, SEQUENCE) %>%
  count(Sample, name = "plasmid_contigs")

# Metric 2: replicon-type RICHNESS per sample
# (counts distinct replicon types, e.g. IncFIA + IncFIB + IncFIC = 3)
replicon_richness <- pf_hits %>%
  distinct(Sample, GENE) %>%
  count(Sample, name = "replicon_richness")

summary_df <- tibble(Sample = full_samples) %>%
  left_join(sample_meta %>% rownames_to_column("Sample"), by = "Sample") %>%
  left_join(contig_counts, by = "Sample") %>%
  left_join(replicon_richness, by = "Sample") %>%
  mutate(
    plasmid_contigs = replace_na(plasmid_contigs, 0),
    replicon_richness = replace_na(replicon_richness, 0)
  )


# 3. Descriptive summary: Badger vs Fox

species_summary <- summary_df %>%
  group_by(Species) %>%
  summarise(
    n_samples = n(),
    n_with_any_plasmid = sum(plasmid_contigs > 0),
    pct_with_any_plasmid = round(100 * mean(plasmid_contigs > 0), 1),
    median_plasmid_contigs = median(plasmid_contigs),
    mean_plasmid_contigs = round(mean(plasmid_contigs), 2),
    median_replicon_richness = median(replicon_richness),
    mean_replicon_richness = round(mean(replicon_richness), 2)
  )

print(species_summary)
print(species_summary, width = Inf)


# 4. Total plasmid contig count, Badger vs Fox - Figure 3.16

p1 <- ggplot(summary_df, aes(x = Species, y = plasmid_contigs, fill = Species)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(x = NULL, y = "Number of plasmid-associated contigs",
       title = "Plasmid contig count per sample: Badger vs Fox") +
  scale_fill_manual(values = c(Badger = "#4C6E5D", Fox = "#C97B3D")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

# replicon-type richness, Badger vs Fox

p2 <- ggplot(summary_df, aes(x = Species, y = replicon_richness, fill = Species)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(x = NULL, y = "Number of distinct replicon types detected",
       title = "Plasmid replicon-type richness per sample: Badger vs Fox") +
  scale_fill_manual(values = c(Badger = "#4C6E5D", Fox = "#C97B3D")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )


# ---- Combine side by side ----
boxes <- p1 + p2 + plot_annotation(tag_levels = "A")

boxes

## -----------------------------------------------------------------------------
## Figure 3.18

# 6. Presence/absence heatmap of replicon types across samples

presence_df <- pf_hits %>%
  distinct(Sample, GENE) %>%
  mutate(present = 1)

pres_mat <- presence_df %>%
  pivot_wider(names_from = GENE, values_from = present, values_fill = 0) %>%
  column_to_rownames("Sample") %>%
  as.matrix()

# Add back samples with zero plasmid hits as all-zero rows
missing_samples <- setdiff(full_samples, rownames(pres_mat))
if (length(missing_samples) > 0) {
  zero_mat <- matrix(0, nrow = length(missing_samples), ncol = ncol(pres_mat),
                     dimnames = list(missing_samples, colnames(pres_mat)))
  pres_mat <- rbind(pres_mat, zero_mat)
}

# Order by species, then by total replicon richness within species
sample_order <- summary_df %>%
  arrange(Species, desc(replicon_richness)) %>%
  pull(Sample)

pres_mat <- pres_mat[sample_order, ]
sample_meta_ordered <- sample_meta[sample_order, , drop = FALSE]
gap_after <- sum(sample_meta_ordered$Species == "Badger")

annotation_colors <- list(Species = c(Badger = "#4C6E5D", Fox = "#C97B3D"))

pheatmap(
  t(pres_mat),  # transpose: samples on x-axis, replicon types on y-axis
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = sample_meta_ordered,
  annotation_colors = annotation_colors,
  gaps_col = gap_after,
  color = c("grey92", "#B2182B"),  # binary: absent / present
  legend_breaks = c(0, 1),
  legend_labels = c("Absent", "Present"),
  border_color = "grey85",
  fontsize = 10,
  fontsize_row = 7,
  fontsize_col = 8,
  main = "Plasmid replicon type presence/absence (Badger vs Fox)",
  treeheight_row = 40
)


# 6b. Heatmap of HOW MANY TIMES each replicon type occurs per sample

# Same structure as the presence/absence heatmap above - just counts
# instead of 1/0, and a color gradient instead of two colors.

count_df <- pf_hits %>%
  count(Sample, GENE, name = "n")

count_mat <- count_df %>%
  pivot_wider(names_from = GENE, values_from = n, values_fill = 0) %>%
  column_to_rownames("Sample") %>%
  as.matrix()

# Add back samples with zero plasmid hits as all-zero rows
missing_samples_count <- setdiff(full_samples, rownames(count_mat))
if (length(missing_samples_count) > 0) {
  zero_mat_count <- matrix(0, nrow = length(missing_samples_count), ncol = ncol(count_mat),
                           dimnames = list(missing_samples_count, colnames(count_mat)))
  count_mat <- rbind(count_mat, zero_mat_count)
}

# Use the SAME sample order as the presence/absence heatmap above
count_mat <- count_mat[sample_order, ]

# Log-scale the color mapping only (real counts still shown in the legend).
# Linear scale washes out most values here since some replicon types
# (e.g. ColRNAI) reach 20-40+ while most others sit at 1-2.
count_mat_log <- log1p(count_mat)
real_ticks <- unique(c(0, 1, 5, 10, 25, max(count_mat)))
real_ticks <- real_ticks[real_ticks <= max(count_mat)]
log_tick_positions <- log1p(real_ticks)

pheatmap(
  t(count_mat_log),  # transpose: samples on x-axis, replicon types on y-axis
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = sample_meta_ordered,
  annotation_colors = annotation_colors,
  gaps_col = gap_after,
  color = colorRampPalette(c("grey92", brewer.pal(9, "YlOrRd")))(100),
  legend_breaks = log_tick_positions,
  legend_labels = real_ticks,
  border_color = "grey85",
  fontsize = 10,
  fontsize_row = 7,
  fontsize_col = 8,
  main = "Plasmid replicon contig counts per sample (Badger vs Fox)",
  treeheight_row = 40
)

# 7. Diagnostic: is this an assembly-quality confound?

quast <- read_csv("Quast_Results.csv", show_col_types = FALSE) %>%
  mutate(Sample = as.character(Sample)) %>%
  select(Sample, total_length = `Total length`, n_contigs = `#contigs`, N50)

diag_df <- summary_df %>% left_join(quast, by = "Sample")

ggplot(diag_df, aes(x = total_length, y = plasmid_contigs, color = Species)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(x = "Total assembly length (bp)", y = "Plasmid contig count",
       title = "Plasmid contigs vs assembly size, by host species") +
  scale_color_manual(values = c(Badger = "#4C6E5D", Fox = "#C97B3D")) +
  theme_minimal()

ggplot(diag_df, aes(x = N50, y = plasmid_contigs, color = Species)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(x = "Assembly N50", y = "Plasmid contig count",
       title = "Plasmid contigs vs N50, by host species") +
  scale_color_manual(values = c(Badger = "#4C6E5D", Fox = "#C97B3D")) +
  theme_minimal()

diag_df %>%
  group_by(Species) %>%
  summarise(
    median_total_length = median(total_length, na.rm = TRUE),
    median_N50 = median(N50, na.rm = TRUE),
    median_n_contigs = median(n_contigs, na.rm = TRUE)
  )

write_csv(summary_df, "plasmid_summary_per_sample.csv")



# Spearman correlation: plasmid contig count vs assembly quality

# Tests whether plasmid detection tracks assembly contiguity/size

# ---- Overall (all samples pooled) ----
cor_n50_overall <- cor.test(diag_df$N50, diag_df$plasmid_contigs,
                            method = "spearman", exact = FALSE)

cor_length_overall <- cor.test(diag_df$total_length, diag_df$plasmid_contigs,
                               method = "spearman", exact = FALSE)

cat("Overall: Plasmid contigs vs N50\n")
cat("  rho =", round(cor_n50_overall$estimate, 3),
    ", p =", signif(cor_n50_overall$p.value, 3), "\n\n")

cat("Overall: Plasmid contigs vs Total assembly length\n")
cat("  rho =", round(cor_length_overall$estimate, 3),
    ", p =", signif(cor_length_overall$p.value, 3), "\n\n")

# ---- Within each species separately ----

spearman_by_species <- diag_df %>%
  group_by(Species) %>%
  summarise(
    n = n(),
    rho_N50 = cor.test(N50, plasmid_contigs, method = "spearman", exact = FALSE)$estimate,
    p_N50 = cor.test(N50, plasmid_contigs, method = "spearman", exact = FALSE)$p.value,
    rho_length = cor.test(total_length, plasmid_contigs, method = "spearman", exact = FALSE)$estimate,
    p_length = cor.test(total_length, plasmid_contigs, method = "spearman", exact = FALSE)$p.value
  )

print(spearman_by_species)


diag_df %>%
  group_by(Species) %>%
  summarise(
    median_total_length = median(total_length, na.rm = TRUE),
    median_N50 = median(N50, na.rm = TRUE),
    median_n_contigs = median(n_contigs, na.rm = TRUE)
  )


## -----------------------------------------------------------------------------

## Genomic context of ARGs: association with plasmid-associated contigs

# 1. Read in data

card <- read_csv("CARD_Ass_ALL.csv", show_col_types = FALSE) %>%
  mutate(Sample = as.character(Sample)) %>%
  filter(!is.na(Sample))

pf <- read_csv("Plasmids_on_Assemblies.csv", show_col_types = FALSE) %>%
  mutate(Sample = as.character(Sample)) %>%
  filter(!is.na(Sample))

card_hits <- card %>% filter(!is.na(Gene))
pf_hits <- pf %>% filter(!is.na(GENE))


# 2. Identify which contigs carry a plasmid replicon marker


# One row per plasmid-positive contig, with all replicon types on it
# combined into a single string (a contig can carry more than one,
# e.g. a multi-replicon plasmid like IncFIA + IncFIB + IncFIC)
plasmid_contigs <- pf_hits %>%
  group_by(Sample, SEQUENCE) %>%
  summarise(
    replicon_types_on_contig = paste(sort(unique(GENE)), collapse = "; "),
    n_replicons_on_contig = n_distinct(GENE),
    .groups = "drop"
  )


# 3. Join ARG hits to plasmid contig info


arg_plasmid_context <- card_hits %>%
  left_join(plasmid_contigs, by = c("Sample", "SEQUENCE")) %>%
  mutate(
    plasmid_associated = !is.na(replicon_types_on_contig),
    contig_length = END - START  # rough proxy; replace with real contig
  )

# 4. Summary: The proportion of ARG hits that are plasmid-associated

overall_summary <- arg_plasmid_context %>%
  summarise(
    n_arg_hits = n(),
    n_plasmid_associated = sum(plasmid_associated),
    pct_plasmid_associated = round(100 * mean(plasmid_associated), 1)
  )

print(overall_summary)

species_summary <- arg_plasmid_context %>%
  group_by(Species) %>%
  summarise(
    n_arg_hits = n(),
    n_plasmid_associated = sum(plasmid_associated),
    pct_plasmid_associated = round(100 * mean(plasmid_associated), 1)
  )

print(species_summary)

# Per-sample breakdown, useful for the bar chart below
per_sample_summary <- arg_plasmid_context %>%
  group_by(Sample, Species) %>%
  summarise(
    n_plasmid_associated = sum(plasmid_associated),
    n_chromosomal = sum(!plasmid_associated),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(n_plasmid_associated, n_chromosomal),
               names_to = "location", values_to = "count") %>%
  mutate(location = recode(location,
                           n_plasmid_associated = "Plasmid-associated",
                           n_chromosomal = "Chromosomal/unidentified"))

# 5. Which specific ARGs are most often plasmid-linked?


arg_by_location <- arg_plasmid_context %>%
  group_by(Gene) %>%
  summarise(
    n_total = n(),
    n_plasmid_associated = sum(plasmid_associated),
    pct_plasmid_associated = round(100 * mean(plasmid_associated), 1)
  ) %>%
  arrange(desc(pct_plasmid_associated))

print(arg_by_location, n = Inf)


# 6. ARG x replicon-type co-occurrence table

# Which specific ARGs turn up linked to which specific plasmid types -
# often one of the more interesting findings in this kind of analysis

arg_replicon_pairs <- arg_plasmid_context %>%
  filter(plasmid_associated) %>%
  separate_rows(replicon_types_on_contig, sep = "; ") %>%
  count(Gene, replicon_types_on_contig, sort = TRUE, name = "n_contigs")

print(arg_replicon_pairs, n = Inf)


# 7. Figure 3.19: plasmid-associated vs chromosomal ARG counts per sample

ggplot(per_sample_summary, aes(x = Sample, y = count, fill = location)) +
  geom_col() +
  facet_grid(. ~ Species, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Number of ARG hits",
       fill = NULL,
       title = "ARG genomic context: plasmid-associated vs chromosomal/unidentified") +
  scale_fill_manual(values = c("Plasmid-associated" = "#B2182B",
                               "Chromosomal/unidentified" = "grey70")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

write_csv(arg_plasmid_context, "arg_plasmid_context_full.csv")
write_csv(arg_by_location, "arg_pct_plasmid_associated.csv")
write_csv(arg_replicon_pairs, "arg_replicon_cooccurrence.csv")


# 4. Do the same plasmid-ARG pairings occur in both host species?

arg_replicon_pairs_by_species <- arg_plasmid_context %>%
  filter(plasmid_associated) %>%
  separate_rows(replicon_types_on_contig, sep = "; ") %>%
  count(Species, Gene, replicon_types_on_contig, sort = TRUE, name = "n_contigs")

print(arg_replicon_pairs_by_species, n = Inf)

# ---- Which pairings are shared vs species-specific? ----
pairing_species_summary <- arg_replicon_pairs_by_species %>%
  group_by(Gene, replicon_types_on_contig) %>%
  summarise(
    species_present = paste(sort(unique(Species)), collapse = ", "),
    n_species = n_distinct(Species),
    total_contigs = sum(n_contigs),
    .groups = "drop"
  ) %>%
  arrange(desc(n_species), desc(total_contigs))

print(pairing_species_summary, n = Inf)

# Quick counts: how many pairings are shared vs unique to one species
pairing_species_summary %>%
  count(species_present, name = "n_pairings") %>%
  arrange(desc(n_pairings))


# 5. Are plasmid-associated ARGs widespread or concentrated in a few samples?


# Number of plasmid-associated ARG hits per sample
plasmid_arg_per_sample <- arg_plasmid_context %>%
  filter(plasmid_associated) %>%
  count(Sample, Species, sort = TRUE, name = "n_plasmid_arg_hits")

print(plasmid_arg_per_sample, n = Inf)

# How many samples actually contribute at least one plasmid-associated hit,
# vs how many total samples exist per species (use full sample list, not
# just those with hits, so "0" samples are visible too)
full_sample_species <- arg_plasmid_context %>% distinct(Sample, Species)

concentration_summary <- full_sample_species %>%
  left_join(plasmid_arg_per_sample, by = c("Sample", "Species")) %>%
  mutate(n_plasmid_arg_hits = replace_na(n_plasmid_arg_hits, 0)) %>%
  group_by(Species) %>%
  summarise(
    n_samples_total = n(),
    n_samples_with_hit = sum(n_plasmid_arg_hits > 0),
    pct_samples_with_hit = round(100 * mean(n_plasmid_arg_hits > 0), 1),
    max_hits_single_sample = max(n_plasmid_arg_hits),
    pct_of_total_from_top_sample = round(100 * max(n_plasmid_arg_hits) / sum(n_plasmid_arg_hits), 1),
    pct_of_total_from_top2_samples = round(100 * sum(sort(n_plasmid_arg_hits, decreasing = TRUE)[1:2]) /
                                             sum(n_plasmid_arg_hits), 1)
  )

print(concentration_summary)


write_csv(arg_replicon_pairs_by_species, "arg_replicon_pairs_by_species.csv")
write_csv(pairing_species_summary, "arg_replicon_pairing_shared_vs_unique.csv")
write_csv(plasmid_arg_per_sample, "plasmid_arg_hits_per_sample.csv")













