library(tidyverse)
library(pheatmap)
library(patchwork)
library(ggplot2)
library(vegan)
library(ggrepel)
library(viridis)
library(stringr)
library(dplyr)


## 1. E. coli AMR profiles - ordination, PERMANOVA & spatial (Mantel) test

# STEP 1. Load data

data     <- read.csv("E.coli_heatmap_updated.csv", check.names = FALSE)
metadata <- read.csv("fox_metadata_all.csv", check.names = FALSE)

# STEP 2. Build IDs/grouping variables, then the AMR matrix

data <- data %>%
  mutate(
    BaseSample  = gsub("-ESBL", "", Sample),
    ESBL_status = ifelse(grepl("ESBL", Sample), "ESBL", "Non-ESBL"),
    Isolate_ID  = make.unique(as.character(Sample))
  )

mat <- data %>%
  column_to_rownames("Isolate_ID") %>%
  select(where(is.numeric)) %>%
  as.matrix()

dim(mat)

# STEP 3. Sanity checks before running any distance-based stats

zero_rows <- rownames(mat)[rowSums(mat) == 0]
if (length(zero_rows) > 0) {
  message("Isolates with resistance to nothing tested (check before trusting Bray-Curtis): ",
          paste(zero_rows, collapse = ", "))
}

missing_meta <- setdiff(data$BaseSample, metadata$Sample)
if (length(missing_meta) > 0) {
  message("These samples have no metadata match: ", paste(missing_meta, collapse = ", "))
}

table(data$BaseSample)[table(data$BaseSample) > 1]  # confirms only the 2 known foxes are duplicated


# STEP 4. Ordination

set.seed(123)
nmds <- metaMDS(mat, distance = "bray", k = 2, trymax = 100)

print(nmds)
nmds$stress
stressplot(nmds)

nmds_points <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("Isolate_ID")

plot_data <- nmds_points %>%
  left_join(data %>% select(Isolate_ID, BaseSample, ESBL_status), by = "Isolate_ID") %>%
  left_join(metadata, by = c("BaseSample" = "Sample"))

stopifnot(nrow(plot_data) == nrow(nmds_points))
plot_data$Apparent_health <- factor(plot_data$Apparent_health, levels = c("Healthy", "Sick"))

dist_mat <- vegdist(mat, method = "bray")  # isolate-level distance matrix (n = 16)


# STEP 5. One isolate per fox - for the HOST-ATTRIBUTE tests only (Sex, Health)

mat_by_fox       <- mat[!grepl("-ESBL", rownames(mat)), ]
plot_data_by_fox <- plot_data %>% filter(!grepl("-ESBL", Isolate_ID))

stopifnot(nrow(mat_by_fox) == nrow(plot_data_by_fox))
stopifnot(!anyDuplicated(plot_data_by_fox$BaseSample))

dim(mat_by_fox)  # 14 foxes
dist_by_fox <- vegdist(mat_by_fox, method = "bray")

table(plot_data_by_fox$Apparent_health)
table(plot_data_by_fox$Sex)

bd_health <- betadisper(dist_by_fox, plot_data_by_fox$Apparent_health)
anova(bd_health)

bd_sex <- betadisper(dist_by_fox, plot_data_by_fox$Sex)
anova(bd_sex)

adonis_health <- adonis2(mat_by_fox ~ Apparent_health, data = plot_data_by_fox, method = "bray", permutations = 9999)
adonis_sex    <- adonis2(mat_by_fox ~ Sex,             data = plot_data_by_fox, method = "bray", permutations = 9999)

adonis_health
adonis_sex


# STEP 6. Location tests - isolate level

table(plot_data$Location)  # now legitimately includes both isolates per duplicated fox

bd_location <- betadisper(dist_mat, plot_data$Location)
anova(bd_location)

adonis_location <- adonis2(mat ~ Location, data = plot_data, method = "bray", permutations = 9999)
adonis_location 


# STEP 7. Spatial (Mantel) test - isolate level

park_coords <- tribble(
  ~Location,                  ~lat,     ~lon,
  "Bushy Park",                53.3013, -6.2900,
  "Clonskeagh",                53.3083, -6.2403,
  "Fairview Park",                  53.3616, -6.2332,
  "Irishtown Nature Reserve",  53.3361, -6.1969,
  "Saint Annes Park",          53.3721, -6.1796,
  "Tolka Valley Park",         53.3770, -6.3010
)

plot_data <- plot_data %>% left_join(park_coords, by = "Location")
stopifnot(!any(is.na(plot_data$lat)))

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  to_rad <- function(x) x * pi / 180
  dlat <- to_rad(lat2 - lat1)
  dlon <- to_rad(lon2 - lon1)
  a <- sin(dlat / 2)^2 + cos(to_rad(lat1)) * cos(to_rad(lat2)) * sin(dlon / 2)^2
  R * 2 * atan2(sqrt(a), sqrt(1 - a))
}

n_iso <- nrow(plot_data)
geo_dist_mat <- matrix(0, n_iso, n_iso)
for (i in seq_len(n_iso)) {
  for (j in seq_len(n_iso)) {
    geo_dist_mat[i, j] <- haversine_km(plot_data$lat[i], plot_data$lon[i],
                                       plot_data$lat[j], plot_data$lon[j])
  }
}
rownames(geo_dist_mat) <- colnames(geo_dist_mat) <- plot_data$Isolate_ID
geo_dist <- as.dist(geo_dist_mat)

set.seed(123)
mantel_result <- mantel(dist_mat, geo_dist, method = "spearman", permutations = 9999)
mantel_result


# STEP 8. ESBL comparison - descriptive only

esbl_pairs <- data %>%
  filter(BaseSample %in% BaseSample[ESBL_status == "ESBL"]) %>%
  arrange(BaseSample, ESBL_status)

esbl_pairs

for (id in unique(esbl_pairs$BaseSample)) {
  rows <- esbl_pairs %>% filter(BaseSample == id)
  d <- vegdist(as.matrix(rows %>% column_to_rownames("Isolate_ID") %>% select(where(is.numeric))),
               method = "bray")
  message("Fox ", id, ": Bray-Curtis distance between paired isolates = ", round(d, 3))
}

summary(dist_mat)  # context for how large those within-fox distances are, population-wide


# STEP 9. Which antibiotics drive the ordination?

mat_filtered <- mat[, apply(mat, 2, sd) > 0]
fit <- envfit(nmds, mat_filtered, permutations = 9999)
fit

plot_data$Year <- format(as.Date(plot_data$Collection_date, format = "%d/%m/%Y"), "%Y")
table(plot_data$Year, plot_data$Location)

bushy <- plot_data %>% filter(Location == "Bushy Park")
mat_bushy <- mat[bushy$Isolate_ID, ]
adonis2(mat_bushy ~ Year, data = bushy, method = "bray", permutations = 9999)


# STEP 10. Figure 5.3 - NMDS ordination

stress_lab <- paste0("Stress = ", round(nmds$stress, 3))

p_nmds <- ggplot(plot_data, aes(NMDS1, NMDS2, colour = Location, shape = ESBL_status)) +
  stat_ellipse(aes(group = Location), type = "t", level = 0.68,
               linetype = 2, linewidth = 0.4, show.legend = FALSE) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = Isolate_ID), size = 3, max.overlaps = 20,
                  box.padding = 0.4, segment.color = "grey60", show.legend = FALSE) +
  coord_equal() +
  # NOTE: viridis "turbo" is a rainbow-style palette, not one of viridis's
  # actual colourblind-vetted options - with 6 categories, adjacent hues
  # (e.g. two different greens) can be hard to tell apart at a glance.
  # scale_colour_brewer's Dark2 is a genuinely qualitative, distinguishable
  # palette for up to 8 groups.
  scale_colour_brewer(palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    axis.title       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, colour = "grey40")
  ) +
  labs(
    title    = "NMDS of E. coli AMR profiles from fox isolates",
    subtitle = stress_lab,
    x = "NMDS1", y = "NMDS2",
    colour = "Park", shape = "ESBL status"
  )

p_nmds
ggsave("nmds_fox_AMR.png", p_nmds, width = 8, height = 6, dpi = 300)



## Figure 5.4 - Scatter across geographic distance

p_mantel <- ggplot(mantel_df, aes(geo_km, amr_bray)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 2, colour = "#1f78b4") +
  geom_smooth(method = "loess", method.args = list(degree = 1),
              se = TRUE, colour = "firebrick", linewidth = 0.6, span = 1) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40")
  ) +
  labs(
    title    = "AMR dissimilarity vs geographic distance between isolates",
    subtitle = paste0("Mantel r (Spearman) = ", round(mantel_result$statistic, 3),
                      ", p = ", round(mantel_result$signif, 3)),
    x = "Geographic distance between parks (km)",
    y = "Bray-Curtis dissimilarity (AMR profile)"
  )

p_mantel
ggsave("mantel_fox_AMR.png", p_mantel, width = 7, height = 5.5, dpi = 300)



## -----------------------------------------------------------------------------

metadata <- read.csv("fox_metadata_seq.csv", check.names = FALSE)
args     <- read.csv("ARGs_on_reads.csv", check.names = FALSE)

## Figure 5.5

# ---- Load data ----
args_urban <- read_csv("ARGs_on_reads.csv")

# ---- Get gene frequency across samples with ARGs only ----
gene_freq <- args_urban %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Sample, Gene) %>%
  count(Gene, name = "n_samples") %>%
  arrange(desc(n_samples))

# ---- Take top 75 most commonly occurring genes ----
top75_genes <- gene_freq %>%
  slice_max(n_samples, n = 75) %>%
  pull(Gene)

# ---- Build presence/absence matrix for samples WITH ARGs ----
args_urban_mat <- args_urban %>%
  filter(!is.na(Gene), Gene != "") %>%
  filter(Gene %in% top75_genes) %>%
  distinct(Sample, Gene) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Gene,
              values_from = present,
              values_fill = 0) %>%
  column_to_rownames("Sample")

# ---- Reorder columns (genes) alphabetically ----
args_urban_mat <- args_urban_mat[, sort(colnames(args_urban_mat))]

# ---- Transpose so genes are rows and samples are columns ----
args_urban_mat_t <- t(args_urban_mat)


# ---- Compute Jaccard-equivalent distance between samples ----
# R's dist(method = "binary") is the Jaccard distance for 0/1 data - it
# ignores cases where both samples lack a gene, so two "gene-poor" samples
# aren't treated as similar just because they're both missing a lot.
# Computed on args_urban_mat (samples as rows) before transposing.
sample_dist <- dist(args_urban_mat, method = "binary")

# ---- Build the Location annotation strip ----
metadata <- read_csv("fox_metadata_seq.csv")

sample_annotation <- metadata %>%
  filter(Sample %in% rownames(args_urban_mat)) %>%
  select(Sample, Location) %>%
  column_to_rownames("Sample")

# ---- Plot heatmap with Jaccard-based clustering + Location annotation ----
args_urban_heatmap <- pheatmap(args_urban_mat_t,
                               color = c("white", "#1B5E3F"),
                               legend = FALSE,
                               cluster_rows = FALSE,
                               cluster_cols = TRUE,
                               clustering_distance_cols = sample_dist,
                               clustering_method = "complete",
                               annotation_col = sample_annotation,
                               border_color = "grey80",
                               fontsize_row = 14,
                               fontsize_col = 18)


## -----------------------------------------------------------------------------

## Figure 5.6

# ---- Load data ----
args_classes <- read.csv("ARGs_classes.csv", stringsAsFactors = FALSE) %>%
  mutate(Class = str_trim(Class))

metadata <- read.csv("fox_metadata_seq.csv", check.names = FALSE)

# ---- Collapse to distinct gene families per sample per class ----
class_counts <- args_classes %>%
  distinct(Sample, Gene, Class) %>%
  group_by(Sample, Class) %>%
  summarise(Count = n_distinct(Gene), .groups = "drop") %>%
  left_join(metadata %>% select(Sample, Location), by = "Sample")

# ---- Set consistent class factor levels (fixes legend order) ----
all_classes <- sort(unique(class_counts$Class))
class_counts$Class <- factor(class_counts$Class, levels = all_classes)

# ---- Order samples by park, then by total ARG family count within park ----
sample_order <- class_counts %>%
  group_by(Sample, Location) %>%
  summarise(total = sum(Count), .groups = "drop") %>%
  arrange(Location, desc(total)) %>%
  pull(Sample)

class_counts$Sample <- factor(class_counts$Sample, levels = sample_order)

# ---- Colour palette - only classes actually present, consistent with earlier figures ----
class_colours <- c(
  "Tetracyclines"               = "#8DD3C7",
  "Aminoglycosides"              = "#FFFFB3",
  "Lincosamides"                  = "#BEBADA",
  "Phenicols"                     = "#FB8072",
  "Sulfonamides"                  = "#80B1D3",
  "Trimethoprim"                  = "#FDB462",
  "beta-lactams"                  = "#B3DE69",
  "MLSB"                          = "#FCCDE5",
  "Macrolides"                    = "#D9D9D9",
  "Oxazolidinones"                = "#BC80BD",
  "MDR"                           = "#CCEBC5",
  "Glycopeptides"                 = "#FFED6F",
  "Fosfomycin"                    = "#E41A1C",
  "Ionophore resistance"          = "#377EB8",
  "Streptogramins (A)"            = "#4DAF4A",
  "Fluoroquinolones/Quinolones"   = "gold"
)

# ---- Plot ----
p1 <- ggplot(class_counts, aes(x = Sample, y = Count, fill = Class)) +
  geom_col(colour = "grey30", linewidth = 0.2, width = 0.75) +
  facet_grid(~ Location, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = class_colours, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 15),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "grey90", colour = "grey30"),
    strip.text          = element_text(face = "bold"),
    legend.title         = element_text(face = "bold"),
    plot.title            = element_text(face = "bold", hjust = 0.5),
    panel.border           = element_rect(colour = "grey30", linewidth = 0.4)
  ) +
  labs(
    x = "Sample",
    y = "Number of distinct ARG families",
    fill = "Antibiotic class",
    title = "ARG classes detected per sample"
  )

p1
ggsave("ARG_classes_per_sample.png", p1, width = 9, height = 6, dpi = 300)



##------------------------------------------------------------------------------

## Figure 5.7a

# 1. Read in data
bracken_phylum <- read.csv("bracken_phylum_merged.csv", check.names = FALSE)

# 2. Filter out Chordata (host DNA) BEFORE calculating relative abundance
bracken_phylum <- bracken_phylum %>%
  filter(name != "Chordata")

# 3. Reshape to long format
phylum_long <- bracken_phylum %>%
  select(-taxonomy_id, -taxonomy_lvl) %>%
  pivot_longer(-name, names_to = "sample", values_to = "reads")

# 4. Clean up sample names
phylum_long <- phylum_long %>%
  mutate(sample = str_remove(sample, "_bracken_phylum"))

# 5. Assign parkland groups
parkland_map <- c(
  "14_53"   = "Fairview Park",
  "12_62"   = "Fairview Park",
  "15_54"   = "Fairview Park",
  "10_60"   = "Fairview Park",
  "11_5003" = "Bushy Park",
  "13_63"   = "Bushy Park",
  "9_57"    = "Bushy Park"
)

phylum_long <- phylum_long %>%
  mutate(parkland = parkland_map[sample])

# 6. Recalculate relative abundance (now Chordata-free)
phylum_long <- phylum_long %>%
  group_by(sample) %>%
  mutate(rel_abund = reads / sum(reads)) %>%
  ungroup()

# 7. Select top N phyla by overall abundance across all samples
N <- 10 

top_taxa <- phylum_long %>%
  group_by(name) %>%
  summarise(total_rel = sum(rel_abund)) %>%
  arrange(desc(total_rel)) %>%
  slice_head(n = N) %>%
  pull(name)

phylum_long <- phylum_long %>%
  mutate(name_grouped = if_else(name %in% top_taxa, name, "Other"))

# Order by abundance, "Other" always last
taxon_order <- phylum_long %>%
  group_by(name_grouped) %>%
  summarise(total = sum(rel_abund)) %>%
  arrange(desc(total)) %>%
  pull(name_grouped)

taxon_order <- c(setdiff(taxon_order, "Other"), "Other")

phylum_long <- phylum_long %>%
  mutate(name_grouped = factor(name_grouped, levels = taxon_order))

# ---- 8. Build a publication-quality color palette ----
n_taxa <- length(levels(phylum_long$name_grouped))
# Combine two qualitative palettes to get enough distinct colors, grey for "Other"
palette <- colorRampPalette(brewer.pal(12, "Paired"))(n_taxa - 1)
palette <- c(palette, "grey70")
names(palette) <- taxon_order

# ---- 9. Plot ----
phyla <- ggplot(phylum_long, aes(x = sample, y = rel_abund, fill = name_grouped)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  facet_wrap(~parkland, scales = "free_x") +
  scale_fill_manual(values = palette) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 10),  # phylum names in italics per convention
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.4)
  )


phyla


# 10. Save
ggsave("phylum_composition_barplot.png", width = 10, height = 6, dpi = 300)
ggsave("phylum_composition_barplot.pdf", width = 10, height = 6)

## -----------------------------------------------------------------------------

## Figure 5.7b

## GENUS LEVEL BAR PLOT ##

library(dplyr)
library(tidyr)
library(ggplot2)

# --- 1. Load and reshape the merged genus Bracken file ---
bracken_genus <- read.csv("bracken_genus_merged.csv", stringsAsFactors = FALSE)

cat("Column names found in file:\n")
print(colnames(bracken_genus))

sample_cols <- setdiff(colnames(bracken_genus), c("name", "taxonomy_id", "taxonomy_lvl"))
cat("\nUsing these as sample columns:\n")
print(sample_cols)

# Filter out Homo
if ("Homo" %in% bracken_genus$name) {
  homo_reads <- bracken_genus %>% filter(name == "Homo") %>% select(all_of(sample_cols))
  cat("\nHomo reads found and removed - counts per sample before removal:\n")
  print(homo_reads)
}
bracken_genus <- bracken_genus %>% filter(name != "Homo")

long_data <- bracken_genus %>%
  select(name, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "reads") %>%
  mutate(sample = gsub("^X|_bracken_genus$", "", sample))  # clean up sample names

# --- Assign parkland groups (same mapping as the phylum script) ---
parkland_map <- c(
  "14_53"   = "Fairview Park",
  "12_62"   = "Fairview Park",
  "15_54"   = "Fairview Park",
  "10_60"   = "Fairview Park",
  "11_5003" = "Bushy Park",
  "13_63"   = "Bushy Park",
  "9_57"    = "Bushy Park"
)
long_data <- long_data %>%
  mutate(parkland = parkland_map[sample])

# --- 2. Convert to relative abundance (proportion 0-1, matching the
#         phylum script, so both plots use the same scale/percent labels)
long_data <- long_data %>%
  group_by(sample) %>%
  mutate(rel_abund = reads / sum(reads)) %>%
  ungroup()

# --- 3. Most common genera overview (by mean relative abundance) ---
genus_summary <- long_data %>%
  group_by(name) %>%
  summarise(
    mean_rel_abund_pct = round(mean(rel_abund) * 100, 2),
    max_rel_abund_pct = round(max(rel_abund) * 100, 2),
    n_samples_present = sum(reads > 0),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_rel_abund_pct))

cat("\n=== Top 20 most common genera (by mean relative abundance) ===\n")
for (i in 1:20) {
  r <- genus_summary[i, ]
  cat(r$name, ": mean =", r$mean_rel_abund_pct, "%, max =", r$max_rel_abund_pct,
      "%, present in", r$n_samples_present, "of 7 samples\n")
}

write.csv(genus_summary, "fox_genus_abundance_summary.csv", row.names = FALSE)
cat("\nSaved fox_genus_abundance_summary.csv\n")

# --- 4. Select top N genera by overall abundance, rest grouped as "Other" ---
N <- 15  # showing top 15 genera individually
top_taxa <- long_data %>%
  group_by(name) %>%
  summarise(total_rel = sum(rel_abund), .groups = "drop") %>%
  arrange(desc(total_rel)) %>%
  slice_head(n = N) %>%
  pull(name)

long_data <- long_data %>%
  mutate(name_grouped = if_else(name %in% top_taxa, name, "Other"))

# --- Explicitly combine all rows into ONE row per sample x genus category
#     (summing rel_abund) before plotting. Without this, every individual
#     rare genus lumped into "Other" stays as its own separate row, and
#     ggplot has to stack potentially dozens of hairline-thin segments -
#     which can visually fragment/vanish rather than rendering as one
#     clean block. ---
long_data <- long_data %>%
  group_by(sample, parkland, name_grouped) %>%
  summarise(rel_abund = sum(rel_abund), .groups = "drop")

cat("\nRows after aggregation (should be samples x categories, e.g. 7 x 16 = 112):",
    nrow(long_data), "\n")

# Order by abundance, "Other" always last
taxon_order <- long_data %>%
  group_by(name_grouped) %>%
  summarise(total = sum(rel_abund), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(name_grouped)
taxon_order <- c(setdiff(taxon_order, "Other"), "Other")

long_data <- long_data %>%
  mutate(name_grouped = factor(name_grouped, levels = taxon_order))

# --- 5. Colour palette ---

if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)
n_taxa <- length(levels(long_data$name_grouped)) - 1  # excluding "Other"

if (n_taxa <= 12) {
  base_colors <- brewer.pal(max(n_taxa, 3), "Paired")[1:n_taxa]
} else {
  # Combine Paired + Set2 + Dark2 for enough genuinely distinct colours.
  # NOTE: Paired's default position 11 is a very pale, near-white yellow
  # (#FFFF99) which disappears against a white background, especially in
  # thin low-abundance segments - swapped out for a more visible gold.
  combined <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))
  combined[combined == "#FFFF99"] <- "#D4A017"  # replace pale yellow with a visible gold
  base_colors <- combined[1:n_taxa]
}

palette <- c(base_colors, "grey70")
names(palette) <- taxon_order

# --- 6. Plot (matching the phylum script's style, including parkland facets) ---
genera <- ggplot(long_data, aes(x = sample, y = rel_abund, fill = name_grouped)) +
  geom_col(width = 0.75) +
  facet_wrap(~parkland, scales = "free_x") +
  scale_fill_manual(values = palette) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Relative Abundance",
    fill = "Genus"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 10, face = "italic"),  # genus names in italics per convention
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.4)
  )

genera

ggsave("genus_composition_barplot.png", genera, width = 10, height = 6, dpi = 300)
ggsave("genus_composition_barplot.pdf", genera, width = 10, height = 6)
cat("\nSaved genus_composition_barplot.png and .pdf\n")




## PATCHWORK to merge ##

library(patchwork)

# --- Stacked vertically (phyla on top, genera below), with A/B panel labels ---
combined_plot <- (phyla / genera) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right", plot.tag = element_text(size = 16, face = "bold"))

combined_plot

ggsave("phylum_genus_combined.png", combined_plot, width = 10, height = 12, dpi = 300)
ggsave("phylum_genus_combined.pdf", combined_plot, width = 12, height = 12)


##------------------------------------------------------------------------------



## Species-level bracken statistical analyses ##

library(tidyverse)
library(vegan)

# ---- 1. Read in species-level data ----
bracken_species <- read.csv("bracken_species_merged.csv", check.names = FALSE)

# ---- 2. Filter out Chordata (host DNA) again at species level ----
# (also worth checking for other host-associated species, e.g. Vulpes vulpes if present)
bracken_species <- bracken_species %>%
  filter(!name %in% c("Chordata", "Vulpes vulpes"))

# ---- 3. Build a sample x species matrix (rows = samples, cols = species) ----
species_matrix <- bracken_species %>%
  select(-taxonomy_id, -taxonomy_lvl) %>%
  column_to_rownames("name") %>%
  t() %>%
  as.data.frame()

rownames(species_matrix) <- str_remove(rownames(species_matrix), "_bracken_species")

# ---- 4. Metadata: parkland assignment ----
parkland_map <- c(
  "14_53"   = "Fairview Park",
  "12_62"   = "Fairview Park",
  "15_54"   = "Fairview Park",
  "10_60"   = "Fairview Park",
  "11_5003" = "Bushy Park",
  "13_63"   = "Bushy Park",
  "9_57"    = "Bushy Park"
)

metadata <- data.frame(
  sample = rownames(species_matrix),
  parkland = parkland_map[rownames(species_matrix)]
)

# ---- 5. Check sequencing depth per sample (important before diversity calcs) ----
depth <- rowSums(species_matrix)
print(depth)
# If depth varies a lot (e.g. >5-10x difference), consider rarefying — see step 5b


# ---- 6. Alpha diversity ----
alpha_div <- data.frame(
  sample = rownames(species_matrix),
  richness = specnumber(species_matrix),
  shannon = diversity(species_matrix, index = "shannon"),
  simpson = diversity(species_matrix, index = "simpson")
) %>%
  left_join(metadata, by = "sample")

print(alpha_div)
alpha_div


# ---- 7. Compare alpha diversity between parklands ----
# Wilcoxon rank-sum test (non-parametric, appropriate for small n)
wilcox_shannon <- wilcox.test(shannon ~ parkland, data = alpha_div)
wilcox_simpson <- wilcox.test(simpson ~ parkland, data = alpha_div)
wilcox_richness <- wilcox.test(richness ~ parkland, data = alpha_div)

print(wilcox_shannon)
print(wilcox_simpson)
print(wilcox_richness)

# STEP 8. Figure S5.3 - Visualize alpha diversity
alpha_long <- alpha_div %>%
  pivot_longer(cols = c(richness, shannon, simpson),
               names_to = "metric", values_to = "value")

ggplot(alpha_long, aes(x = parkland, y = value, fill = parkland)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c("Fairview Park" = "#F4A9A0", "Bushy Park" = "#7FDBDA")) +
  labs(x = NULL, y = "Diversity Index Value", title = "Alpha Diversity by Parkland") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave("alpha_diversity_boxplot.png", width = 9, height = 6, dpi = 300)

# STEP 9. Beta diversity: Bray-Curtis dissimilarity ----
# Convert to relative abundance first
species_rel <- species_matrix / rowSums(species_matrix)

bray_dist <- vegdist(species_rel, method = "bray")

# ---- 10. PCoA ordination ----
pcoa <- cmdscale(bray_dist, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  sample = rownames(species_rel),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  left_join(metadata, by = "sample")

# % variance explained by each axis
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]), 1)

## Figure S5.4

ggplot(pcoa_df, aes(x = PC1, y = PC2, color = parkland)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.4) +
  geom_point(size = 4.5, alpha = 0.9) +
  scale_color_manual(values = c("Fairview Park" = "#F4A9A0", "Bushy Park" = "#7FDBDA")) +
  labs(
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)"),
    title = "Community Composition by Parkland",
    subtitle = "PCoA of Bray-Curtis dissimilarity, species-level",
    color = "Parkland"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.position = "right",
    axis.title = element_text(size = 12)
  )

ggsave("pcoa_beta_diversity_final.png", width = 8, height = 6, dpi = 300)
ggsave("pcoa_beta_diversity_final.pdf", width = 8, height = 6)



# ---- 11. PERMANOVA — test whether parklands differ significantly ----
permanova_result <- adonis2(bray_dist ~ parkland, data = metadata, permutations = 999)
print(permanova_result)

# ---- 12. PERMDISP — check homogeneity of within-group dispersion ----
# (important: PERMANOVA can be confounded by unequal within-group variance)
dispersion <- betadisper(bray_dist, metadata$parkland)
anova(dispersion)
permutest(dispersion, permutations = 999)




##############################################################################

## effects of metadata


library(tidyverse)
library(vegan)

# ---- 1. Read in metadata ----
metadata_full <- read.csv("fox_metadata_seq.csv", check.names = FALSE)
# Expected columns: Sample, Weight_Kg, Apparent_health, Sex, Location

# ---- 2. Merge with alpha diversity results ----
alpha_div <- alpha_div %>%
  select(-parkland) %>%  # drop old parkland column, we'll get it fresh from metadata
  left_join(metadata_full, by = c("sample" = "Sample"))

# ---- 3. Alpha diversity vs WEIGHT (continuous) — correlation ----
cor_shannon_weight <- cor.test(alpha_div$shannon, alpha_div$Weight_Kg, method = "spearman")
cor_simpson_weight <- cor.test(alpha_div$simpson, alpha_div$Weight_Kg, method = "spearman")
cor_richness_weight <- cor.test(alpha_div$richness, alpha_div$Weight_Kg, method = "spearman")

print(cor_shannon_weight)
print(cor_simpson_weight)
print(cor_richness_weight)

# Quick scatter plot
ggplot(alpha_div, aes(x = Weight_Kg, y = shannon)) +
  geom_point(size = 3, color = "#7FDBDA") +
  geom_smooth(method = "lm", se = TRUE, color = "grey40") +
  labs(x = "Weight (kg)", y = "Shannon Diversity", title = "Shannon Diversity vs Weight") +
  theme_minimal(base_size = 13)

# ---- 4. Alpha diversity vs HEALTH STATUS (categorical) ----
table(alpha_div$Apparent_health)

kruskal_shannon_health <- kruskal.test(shannon ~ Apparent_health, data = alpha_div)
kruskal_simpson_health <- kruskal.test(simpson ~ Apparent_health, data = alpha_div)
print(kruskal_shannon_health)
print(kruskal_simpson_health)


wilcox_shannon_health <- wilcox.test(shannon ~ Apparent_health, data = alpha_div)
wilcox_simpson_health <- wilcox.test(simpson ~ Apparent_health, data = alpha_div)
print(wilcox_shannon_health)
print(wilcox_simpson_health)



# ---- 5. Alpha diversity vs SEX (categorical, 2 levels) ----
wilcox_shannon_sex <- wilcox.test(shannon ~ Sex, data = alpha_div)
wilcox_simpson_sex <- wilcox.test(simpson ~ Sex, data = alpha_div)
print(wilcox_shannon_sex)
print(wilcox_simpson_sex)



table(alpha_div$Sex)

table(alpha_div$Location, alpha_div$Apparent_health)
table(alpha_div$Location, alpha_div$Sex)
table(alpha_div$Sex, alpha_div$Apparent_health)


ggplot(alpha_div, aes(x = Apparent_health, y = simpson, fill = Apparent_health)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2) +
  labs(x = NULL, y = "Simpson Diversity", title = "Simpson Diversity by Health Status") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggplot(alpha_div, aes(x = Sex, y = simpson, fill = Sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2) +
  labs(x = NULL, y = "Simpson Diversity", title = "Simpson Diversity by Sex") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

alpha_div %>% select(sample, Location, Apparent_health, Sex, simpson) %>% arrange(simpson)


# ---- 6. Beta diversity vs each metadata variable — individual PERMANOVA models ----
# IMPORTANT: test each variable in its own model given small n — do not combine into one model

metadata_beta <- metadata_full %>%
  filter(Sample %in% rownames(species_matrix)) %>%
  arrange(match(Sample, rownames(species_matrix)))  # ensure same order as species_matrix/bray_dist

permanova_weight <- adonis2(bray_dist ~ Weight_Kg, data = metadata_beta, permutations = 999)
permanova_health <- adonis2(bray_dist ~ Apparent_health, data = metadata_beta, permutations = 999)
permanova_sex <- adonis2(bray_dist ~ Sex, data = metadata_beta, permutations = 999)

print(permanova_weight)
print(permanova_health)
print(permanova_sex)


# ---- 7. Visualize PCoA colored by each variable ----
pcoa_df <- pcoa_df %>%
  left_join(metadata_full, by = c("sample" = "Sample"))

# By health status
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Apparent_health)) +
  geom_point(size = 4) +
  labs(title = "PCoA colored by Apparent Health") +
  theme_minimal(base_size = 13)

# By sex
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 4) +
  labs(title = "PCoA colored by Sex") +
  theme_minimal(base_size = 13)

# By weight (continuous — use color gradient)
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Weight_Kg)) +
  geom_point(size = 4) +
  scale_color_viridis_c() +
  labs(title = "PCoA colored by Weight") +
  theme_minimal(base_size = 13)


## -----------------------------------------------------------------------------


## Figure S5.5 - Dendrogram


install.packages("dendextend")
library(dendextend)

# 1. Hierarchical clustering on Bray-Curtis distance matrix

hc <- hclust(bray_dist, method = "average")

# ---- 2. Convert to a dendrogram object ----
dend <- as.dendrogram(hc)

# ---- 3. Match label order to metadata, then assign colors by parkland ----
labels_order <- labels(dend)
parkland_ordered <- metadata$parkland[match(labels_order, metadata$sample)]

color_map <- c("Fairview Park" = "#F4A9A0", "Bushy Park" = "#7FDBDA")
label_colors <- color_map[parkland_ordered]

# ---- 4. Style the dendrogram ----
dend <- dend %>%
  set("labels_colors", label_colors) %>%   # corrected from labels_col
  set("labels_cex", 1.3) %>%
  set("branches_lwd", 1.5)

# ---- 5. Plot
par(mar = c(10, 4, 4, 2), font = 2)  # font = 2 makes text bold
plot(dend, main = "Hierarchical Clustering of Fox Microbiome (Bray-Curtis, UPGMA)",
     ylab = "Height (Bray-Curtis dissimilarity)")
legend("topright", legend = names(color_map), fill = color_map, border = NA, bty = "n")

# 6. Save 
png("dendrogram_bray_curtis_upgma.png", width = 10, height = 6, units = "in", res = 300)
par(mar = c(10, 4, 4, 2), font = 2)
plot(dend, main = "Hierarchical Clustering of Fox Microbiome (Bray-Curtis, UPGMA)",
     ylab = "Height (Bray-Curtis dissimilarity)")
legend("topright", legend = names(color_map), fill = color_map, border = NA, bty = "n")
dev.off()



## -----------------------------------------------------------------------------

## Figure 5.8

library(dplyr)
library(ggplot2)

vgs_fox <- read.csv("vgs_fox_with_genus.csv", stringsAsFactors = FALSE)

# --- 1. Select top genera by sample prevalence ---
n_top_genera <- 17
n_genes_per_genus <- 6  # top genes to pull FROM EACH genus, not overall

top_genera <- vgs_fox %>%
  filter(!is.na(genus_name)) %>%
  distinct(Sample, genus_name) %>%
  count(genus_name, sort = TRUE, name = "n_samples") %>%
  slice_head(n = n_top_genera) %>%
  pull(genus_name)

# For each top genus, find ITS OWN top genes by sample prevalence
top_genes <- vgs_fox %>%
  filter(genus_name %in% top_genera) %>%
  distinct(Sample, GENE, genus_name) %>%
  count(genus_name, GENE, name = "n_samples") %>%
  group_by(genus_name) %>%
  slice_max(n_samples, n = n_genes_per_genus, with_ties = FALSE) %>%
  ungroup() %>%
  pull(GENE) %>%
  unique()

cat("Selected", length(top_genes), "genes across", n_top_genera, "genera\n")

# 2. Build gene x genus matrix: number of SAMPLES with each combination
heatmap_data <- vgs_fox %>%
  filter(GENE %in% top_genes, genus_name %in% top_genera) %>%
  distinct(Sample, GENE, genus_name) %>%
  count(GENE, genus_name, name = "n_samples")

full_grid <- expand.grid(GENE = top_genes, genus_name = top_genera, stringsAsFactors = FALSE)
heatmap_data <- full_grid %>%
  left_join(heatmap_data, by = c("GENE", "genus_name")) %>%
  mutate(n_samples = ifelse(is.na(n_samples), 0, n_samples))

# Order genera by overall prevalence; order genes by which genus they were

gene_order <- vgs_fox %>%
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

# 3. Plot
p <- ggplot(heatmap_data, aes(x = genus_name, y = GENE, fill = n_samples)) +
  geom_tile(color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_samples > 0, n_samples, "")),
            size = 3, color = ifelse(heatmap_data$n_samples > 5, "white", "black")) +
  scale_fill_gradient(low = "white", high = "#1B4D3E", name = "No. of\nsamples",
                      trans = "sqrt") +
  labs(
    title = "Virulence gene carriage by bacterial genus (fox samples)",
    x = "Genus", y = "Virulence gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 13),
    axis.text.y = element_text(face = "italic", size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.title.position = "plot",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

ggsave("fox_VG_genus_heatmap.png", p, width = 10, height = 12, dpi = 400)
ggsave("fox_VG_genus_heatmap_17.pdf", p, width = 10, height = 12)
print(p)
cat("\nSaved fox_VG_genus_heatmap.png and .pdf\n")

## -----------------------------------------------------------------------------


# Fox-only dataset (7 samples): 
# 1. Map contigs to genus using Kraken2 report/output files
# 2. Independently link ARGs, VGs, and Plasmids to genus


# Reuse the genus-mapping function from the original analysis
map_contigs_to_genus <- function(report_path, output_path) {
  
  report <- read.delim(report_path, header = FALSE, sep = "\t",
                       col.names = c("pct", "reads_clade", "reads_direct",
                                     "rank_code", "taxid", "name_indented"),
                       stringsAsFactors = FALSE)
  report$indent_level <- str_count(str_extract(report$name_indented, "^\\s*"), "  ")
  report$name <- str_trim(report$name_indented)
  
  n <- nrow(report)
  genus_taxid_vec <- integer(n)
  genus_name_vec <- character(n)
  current_genus_taxid <- NA_integer_
  current_genus_name <- NA_character_
  last_genus_indent <- -1
  
  for (i in seq_len(n)) {
    lvl <- report$indent_level[i]
    rank <- report$rank_code[i]
    tid <- report$taxid[i]
    nm <- report$name[i]
    
    if (lvl <= last_genus_indent) {
      current_genus_taxid <- NA_integer_
      current_genus_name <- NA_character_
      last_genus_indent <- -1
    }
    if (rank == "G") {
      current_genus_taxid <- tid
      current_genus_name <- nm
      last_genus_indent <- lvl
    }
    genus_taxid_vec[i] <- current_genus_taxid
    genus_name_vec[i] <- current_genus_name
  }
  
  taxid_to_genus <- data.frame(taxid = report$taxid, genus_taxid = genus_taxid_vec,
                               genus_name = genus_name_vec)
  
  contig_output <- read.delim(output_path, header = FALSE, sep = "\t",
                              col.names = c("classified", "contig_id", "taxon_field",
                                            "length", "lca_map"),
                              stringsAsFactors = FALSE)
  
  if (all(grepl("^\\d+$", contig_output$taxon_field[1:min(20, nrow(contig_output))]))) {
    contig_output$taxid <- as.integer(contig_output$taxon_field)
  } else {
    contig_output$taxid <- as.integer(str_extract(contig_output$taxon_field, "(?<=taxid )\\d+"))
  }
  
  result <- contig_output %>%
    left_join(taxid_to_genus, by = "taxid") %>%
    dplyr::select(contig_id, classified, length, assigned_taxid = taxid,
                  assigned_taxon = taxon_field, genus_taxid, genus_name)
  
  return(result)
}

# 1. Run genus mapping across all 7 fox samples

kraken_dir <- "C:/PATH/TO/KRAKEN_FILES"

reports <- list.files(kraken_dir, pattern = "_kraken2_report.txt$", full.names = FALSE)
cat("Found", length(reports), "report files (should be 7)\n\n")

for (r in reports) {
  sample <- sub("_kraken2_report.txt$", "", r)
  o <- paste0(sample, "_kraken2_output.txt")
  if (file.exists(o)) {
    cat("Processing", sample, "...\n")
    tbl <- map_contigs_to_genus(r, o)
    tbl$Sample <- sample
    tbl <- tbl[, c("Sample", setdiff(names(tbl), "Sample"))]  # Sample as first column
    write.csv(tbl, paste0(sample, "_contig_genus_assignments.csv"), row.names = FALSE)
  } else {
    cat("Skipping", sample, "- output file not found\n")
  }
}

# --- 2. Combine all 7 samples' genus assignments into one table
genus_files <- list.files(kraken_dir, pattern = "_contig_genus_assignments.csv$", full.names = TRUE)
genus_all_fox <- lapply(genus_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  # Sample column should already be present from step 1 above, but this
  # re-derives it from the filename as a safety net in case any file
  # lacks it (e.g. if step 1 was skipped and files were sourced elsewhere)
  if (!"Sample" %in% names(df)) {
    df$Sample <- sub("_contig_genus_assignments.csv$", "", basename(f))
  }
  df
})
genus_all_fox <- do.call(rbind, genus_all_fox)
genus_all_fox <- genus_all_fox[, c("Sample", setdiff(names(genus_all_fox), "Sample"))]

cat("\nCombined genus table:", nrow(genus_all_fox), "contigs across",
    n_distinct(genus_all_fox$Sample), "samples\n")

write.csv(genus_all_fox, "genus_all_fox.csv", row.names = FALSE)

# 3. Sample IDs in this dataset are strings (e.g. "9_57"), not plain
#         numbers - keep as character, no integer conversion needed here
#         (unlike the earlier badger/fox combined dataset) ---

# 4. Link ARGs, VGs, and Plasmids to genus - INDEPENDENTLY, not merged ---
# NOTE: ARGs_on_contigs.csv uses "Contig" as its contig ID column name,
# while VGs_on_contigs.csv and Plasmids_on_contigs.csv use "SEQUENCE" -
# the contig_col argument below handles this difference.
link_to_genus <- function(filename, label, contig_col = "SEQUENCE") {
  data <- read.csv(filename, stringsAsFactors = FALSE)
  # Sample kept as character (e.g. "9_57") - no integer conversion needed
  
  join_spec <- c("Sample", "contig_col_placeholder" = "contig_id")
  names(join_spec)[2] <- contig_col
  
  linked <- data %>%
    inner_join(genus_all_fox, by = join_spec)
  
  cat("\n===", label, "===\n")
  cat("Hits matched to a contig:", nrow(linked), "of", nrow(data), "\n")
  cat("Of these, with a genus assignment:", sum(!is.na(linked$genus_name)), "\n")
  cat("Samples with at least one hit:", n_distinct(linked$Sample), "\n")
  cat("Unique genera involved:", n_distinct(linked$genus_name, na.rm = TRUE), "\n")
  
  out_file <- paste0(tolower(label), "_fox_with_genus.csv")
  write.csv(linked, out_file, row.names = FALSE)
  cat("Saved to", out_file, "\n")
  
  return(linked)
}

args_fox <- link_to_genus("ARGs_on_contigs.csv", "ARGs", contig_col = "Contig")
vgs_fox <- link_to_genus("VGs_on_contigs.csv", "VGs", contig_col = "SEQUENCE")
plasmids_fox <- link_to_genus("Plasmids_on_contigs.csv", "Plasmids", contig_col = "SEQUENCE")

cat("\n\nAll done. Three independently-linked files saved:\n")
cat(" - args_fox_with_genus.csv\n")
cat(" - vgs_fox_with_genus.csv\n")
cat(" - plasmids_fox_with_genus.csv\n")
cat("\n(Plasmid-ARG merging held off for now, as planned.)\n")



## COUNT 

# Summarise contig classification rates per sample: total contigs,
# % unclassified, % classified but above genus level, % classified to
# genus level - same approach used for the badger/fox combined dataset.

library(dplyr)

classification_summary <- genus_all_fox %>%
  group_by(Sample) %>%
  summarise(
    total_contigs = n(),
    n_unclassified = sum(classified == "U"),
    n_classified_above_genus = sum(classified == "C" & is.na(genus_name)),
    n_classified_to_genus = sum(classified == "C" & !is.na(genus_name)),
    .groups = "drop"
  ) %>%
  mutate(
    pct_unclassified = round(n_unclassified / total_contigs * 100, 1),
    pct_classified_above_genus = round(n_classified_above_genus / total_contigs * 100, 1),
    pct_classified_to_genus = round(n_classified_to_genus / total_contigs * 100, 1)
  )

print(classification_summary)

cat("\n=== Averages across all", nrow(classification_summary), "fox samples ===\n")
cat("Mean total contigs per sample:", round(mean(classification_summary$total_contigs), 1), "\n")
cat("Mean % unclassified:", round(mean(classification_summary$pct_unclassified), 1), "%\n")
cat("Mean % classified but above genus level:",
    round(mean(classification_summary$pct_classified_above_genus), 1), "%\n")
cat("Mean % classified to genus level:",
    round(mean(classification_summary$pct_classified_to_genus), 1), "%\n")
cat("Range of % classified to genus level:",
    round(min(classification_summary$pct_classified_to_genus), 1), "-",
    round(max(classification_summary$pct_classified_to_genus), 1), "%\n")

write.csv(classification_summary, "fox_classification_summary.csv", row.names = FALSE)



## Count of genera occurring across the samples

# Count genus occurrence across the 7 fox samples: how many samples each
# genus appears in, plus total contig count for context.

library(dplyr)

genus_prevalence_fox <- genus_all_fox %>%
  filter(!is.na(genus_name)) %>%
  group_by(genus_name) %>%
  summarise(
    n_samples = n_distinct(Sample),
    total_contigs = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples), desc(total_contigs))

cat("Total unique genera detected across the 7 fox samples:",
    nrow(genus_prevalence_fox), "\n\n")

print(genus_prevalence_fox, n = 50)

write.csv(genus_prevalence_fox, "genus_prevalence_fox.csv", row.names = FALSE)
cat("\nSaved to genus_prevalence_fox.csv\n")



############################################################################

## Count of the virulence genes as associated with genera

# Summarise virulence gene (VG) counts across the 7 fox samples

library(dplyr)

# Reuse vgs_fox if already in session, otherwise reload:
# vgs_fox <- read.csv("vgs_fox_with_genus.csv", stringsAsFactors = FALSE)

cat("=== Overview ===\n")
cat("Total VG hits:", nrow(vgs_fox), "\n")
cat("Unique VG genes:", n_distinct(vgs_fox$GENE), "\n")
cat("Samples with at least one VG hit:", n_distinct(vgs_fox$Sample), "of 7\n")
cat("Unique genera carrying VGs:", n_distinct(vgs_fox$genus_name, na.rm = TRUE), "\n\n")

# --- Per-sample VG counts ---
vg_per_sample <- vgs_fox %>%
  group_by(Sample) %>%
  summarise(n_hits = n(), n_unique_genes = n_distinct(GENE), .groups = "drop") %>%
  arrange(desc(n_hits))

cat("=== VG hits per sample ===\n")
print(vg_per_sample)

# --- Top VG genes by sample prevalence ---
top_vgs_fox <- vgs_fox %>%
  distinct(Sample, GENE) %>%
  count(GENE, sort = TRUE, name = "n_samples")

cat("\n=== Top VG genes by sample prevalence ===\n")
print(head(top_vgs_fox, 20))

# --- Top genera carrying VGs, by sample prevalence ---
top_genera_vgs_fox <- vgs_fox %>%
  filter(!is.na(genus_name)) %>%
  distinct(Sample, genus_name) %>%
  count(genus_name, sort = TRUE, name = "n_samples")

cat("\n=== Top genera carrying VGs, by sample prevalence ===\n")
print(head(top_genera_vgs_fox, 20))

write.csv(vg_per_sample, "vg_per_sample_fox.csv", row.names = FALSE)
write.csv(top_vgs_fox, "top_vgs_fox_prevalence.csv", row.names = FALSE)
write.csv(top_genera_vgs_fox, "top_genera_vgs_fox_prevalence.csv", row.names = FALSE)
cat("\nSaved vg_per_sample_fox.csv, top_vgs_fox_prevalence.csv, top_genera_vgs_fox_prevalence.csv\n")



# Count unique virulence genes (VGs) across the entire VGs_on_contigs.csv
# dataset (all 7 fox samples combined)

vgs_raw <- read.csv("VGs_on_contigs.csv", stringsAsFactors = FALSE)

cat("Total VG hits (rows) in file:", nrow(vgs_raw), "\n")
cat("Unique virulence genes across the whole dataset:", length(unique(vgs_raw$GENE)), "\n")






# Extract a single contig from an assembly FASTA - rebuilt from scratch
# with diagnostics at every step to catch any errors before they produce
# a bad output file.

# Example contig check
Sample <- "SAMPLE_NAME"
target_contig <- "contig_X"
assembly_path <- "/PATH/TO/assembly.fasta"
# -----------------------------

# --- Step 1: confirm the file exists and can be read ---
if (!file.exists(assembly_path)) {
  stop("File not found at: ", assembly_path)
}
cat("File found. Size:", round(file.info(assembly_path)$size / 1e6, 2), "MB\n")

lines <- readLines(assembly_path, warn = FALSE)
cat("Total lines read:", length(lines), "\n")

# --- Step 2: strip any hidden carriage returns (Windows line endings) ---
lines <- gsub("\r", "", lines)

# --- Step 3: find ALL header lines, and check for duplicates ---
header_positions <- grep("^>", lines)
cat("Total sequences (headers) in this file:", length(header_positions), "\n")

header_ids <- sub("^>(\\S+).*", "\\1", lines[header_positions])

# Check specifically how many times our target contig ID appears
target_matches <- which(header_ids == target_contig)
cat("\nNumber of headers matching '", target_contig, "': ", length(target_matches), "\n", sep = "")

if (length(target_matches) == 0) {
  cat("\nContig not found under that exact name. Here are the first 10 headers\n")
  cat("actually present in this file, for comparison:\n")
  print(head(header_ids, 10))
  stop("Stopping - contig ID not found. Check spelling/formatting against the list above.")
}

if (length(target_matches) > 1) {
  cat("\nWARNING: this contig ID appears MORE THAN ONCE in this file.\n")
  cat("Line numbers of each match:", header_positions[target_matches], "\n")
  stop("Stopping - duplicate contig ID found. This needs to be resolved before extracting.")
}

# --- Step 4: extract the sequence lines belonging to this ONE confirmed header ---
match_position_in_list <- target_matches[1]
start_line <- header_positions[match_position_in_list]

is_last_header <- match_position_in_list == length(header_positions)
end_line <- if (is_last_header) length(lines) else header_positions[match_position_in_list + 1] - 1


sequence_lines <- lines[(start_line + 1):end_line]
full_sequence <- paste(sequence_lines, collapse = "")
full_sequence <- toupper(full_sequence)


# --- Step 5: validate the sequence only contains real DNA characters ---
non_dna_chars <- gsub("[ACGTN]", "", full_sequence)
if (nchar(non_dna_chars) > 0) {
  cat("\nWARNING: found", nchar(non_dna_chars), "unexpected non-ACGTN characters:\n")
  print(table(strsplit(non_dna_chars, "")[[1]]))
  stop("Stopping - sequence contains invalid characters. Do not use this output for BLAST.")
} else {
}

# Step 6: trim to BLASTn's 1,000,000 bp limit if needed
max_length <- 500000
if (nchar(full_sequence) > max_length) {
  cat("\nContig exceeds", max_length, "bp - trimming to the first", max_length, "bp\n")
  full_sequence <- substr(full_sequence, 1, max_length)
  trimmed_note <- paste0("_TRIMMED_to_", max_length, "bp")
} else {
  trimmed_note <- ""
}

# --- Step 7: write out as properly-formatted FASTA (70 char line wrap) ---
wrap_width <- 70
n <- nchar(full_sequence)
starts <- seq(1, n, by = wrap_width)
wrapped <- substring(full_sequence, starts, pmin(starts + wrap_width - 1, n))

out_file <- paste0(target_contig, trimmed_note, ".fasta")
writeLines(c(paste0(">", target_contig, trimmed_note), wrapped), out_file)




## count number of ARGs per sample and number of unique ARGs per sample

# Using the RAW ARG file (all ARGs, not restricted to those matched to a
# genus-classified contig)
args_raw <- read.csv("ARGs_on_contigs.csv", stringsAsFactors = FALSE)

per_sample_summary <- args_raw %>%
  group_by(Sample) %>%
  summarise(
    n_ARG_hits = n(),
    n_unique_ARGs = n_distinct(GENE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_ARG_hits))

cat("=== ARG hits and unique ARGs per sample (all ARGs) ===\n\n")
for (i in seq_len(nrow(per_sample_summary))) {
  r <- per_sample_summary[i, ]
  cat(r$Sample, ": ", r$n_ARG_hits, " hits, ", r$n_unique_ARGs, " unique ARGs\n", sep = "")
}

cat("\n=== Totals across all samples ===\n")
cat("Total ARG hits:", nrow(args_raw), "\n")
cat("Total unique ARGs (dataset-wide):", n_distinct(args_raw$GENE), "\n")

write.csv(per_sample_summary, "fox_ARG_counts_per_sample_ALL.csv", row.names = FALSE)
cat("\nSaved to fox_ARG_counts_per_sample_ALL.csv\n")


# Summarise ARGs detected on Kraken2-classified contigs (fox dataset)

library(dplyr)

args_fox <- read.csv("args_fox_with_genus.csv", stringsAsFactors = FALSE)

cat("Total ARG hits (all contigs, classified or not):", nrow(args_fox), "\n\n")

# --- Break down by classification status ---
tryCatch(print(table(args_fox$classified, useNA = "ifany")), error = function(e) {
  cat("(console print failed due to session glitch - try options(na.print = NULL)\n")
  cat(" or restart R, then re-run this script)\n")
})


# --- Count of each ARG gene detected ---
gene_counts <- classified_args %>%
  count(GENE, sort = TRUE, name = "n_hits")

# --- Count of each resistance class detected ---
class_counts <- classified_args %>%
  count(RESISTANCE, sort = TRUE, name = "n_hits")

# Save first, so results are safely on disk even if printing fails
write.csv(gene_counts, "fox_ARG_gene_counts.csv", row.names = FALSE)
write.csv(class_counts, "fox_ARG_resistance_class_counts.csv", row.names = FALSE)

tryCatch(print(gene_counts, n = 50), error = function(e) {
  cat("(console print failed - open fox_ARG_gene_counts.csv to view results)\n")
})

tryCatch(print(class_counts, n = 50), error = function(e) {
  cat("(console print failed - open fox_ARG_resistance_class_counts.csv to view results)\n")
})

gene_counts <- read.csv("fox_ARG_gene_counts.csv")
for (i in 1:nrow(gene_counts)) {
  cat(gene_counts$GENE[i], ":", gene_counts$n_hits[i], "\n")
}


# Check whether excluding ARG hits on contigs without a genus assignment
# (classified, but only to a higher taxonomic rank) changes the count of
# UNIQUE ARG genes - i.e. are any genes found ONLY on those 18 rows,
# and nowhere else among genus-assigned contigs?


args_fox <- read.csv("args_fox_with_genus.csv", stringsAsFactors = FALSE)

classified <- args_fox %>% filter(classified == "C")
genus_assigned <- classified %>% filter(!is.na(genus_name))
no_genus <- classified %>% filter(is.na(genus_name))

genes_all_classified <- unique(classified$GENE)
genes_genus_assigned <- unique(genus_assigned$GENE)


# Genes that would be LOST if restricting to genus-assigned contigs only
lost_genes <- setdiff(genes_all_classified, genes_genus_assigned)

if (length(lost_genes) == 0) {
  cat("No genes are lost - every gene found on the 'no genus' contigs is\n")
  cat("ALSO found on at least one genus-assigned contig elsewhere.\n")
  cat("So the unique gene COUNT is the same either way, even though the\n")
  cat("total ROW count differs (908 vs 890).\n")
} else {
  cat("WARNING:", length(lost_genes), "gene(s) would be lost if restricting\n")
  cat("to genus-assigned contigs only - these genes ONLY appear on contigs\n")
  cat("without a genus assignment:\n")
  print(lost_genes)
}



## Cross-reference ARGs on plasmids

args_fox <- read.csv("args_fox_with_genus.csv", stringsAsFactors = FALSE)
plasmids_fox <- read.csv("plasmids_fox_with_genus.csv", stringsAsFactors = FALSE)

cat("ARG hits:", nrow(args_fox), "\n")
cat("Plasmid replicon hits:", nrow(plasmids_fox), "\n\n")

# --- Find contigs present in BOTH datasets (same Sample + same contig) ---
# NOTE: ARGs file uses "Contig" as its column name; Plasmids file uses
# "SEQUENCE" - matching this the same way as the earlier badger analysis.
plasmid_contigs <- plasmids_fox %>% distinct(Sample, SEQUENCE)

plasmid_mediated_args <- args_fox %>%
  inner_join(plasmid_contigs, by = c("Sample" = "Sample", "Contig" = "SEQUENCE"))

cat("=== Plasmid-mediated ARGs ===\n")
cat("ARG hits found on a contig that ALSO carries a plasmid replicon marker:",
    nrow(plasmid_mediated_args), "of", nrow(args_fox), "total ARG hits\n\n")

if (nrow(plasmid_mediated_args) > 0) {
  cat("Samples involved:", n_distinct(plasmid_mediated_args$Sample), "\n")
  cat("Unique ARG genes:", n_distinct(plasmid_mediated_args$GENE), "\n")
  cat("Unique genera involved:", n_distinct(plasmid_mediated_args$genus_name, na.rm = TRUE), "\n\n")
  
  # --- Pull in the actual plasmid replicon type for each matched contig ---
  plasmid_detail <- plasmids_fox %>%
    distinct(Sample, SEQUENCE, GENE) %>%
    rename(plasmid_replicon_gene = GENE)
  
  plasmid_mediated_args_detail <- plasmid_mediated_args %>%
    left_join(plasmid_detail, by = c("Sample" = "Sample", "Contig" = "SEQUENCE"))
  
  cat("=== Full detail ===\n")
  for (i in seq_len(nrow(plasmid_mediated_args_detail))) {
    r <- plasmid_mediated_args_detail[i, ]
    cat(r$Sample, "|", r$Contig, "|", r$genus_name, "| ARG:", r$GENE,
        "| Plasmid marker:", r$plasmid_replicon_gene, "\n")
  }
  
  write.csv(plasmid_mediated_args_detail, "fox_plasmid_mediated_ARGs.csv", row.names = FALSE)
  cat("\nSaved to fox_plasmid_mediated_ARGs.csv\n")
} else {
}



## -----------------------------------------------------------------------------

## Figure 5.9

library(dplyr)
library(ggplot2)
library(stringr)
library(dplyr)


args_fox <- read.csv("args_fox_with_genus.csv", stringsAsFactors = FALSE)

# --- 1. Select top genera by sample prevalence ---
n_top_genera <- 37
n_genes_per_genus <- 3  # top genes to pull FROM EACH genus, not overall

top_genera <- args_fox %>%
  filter(!is.na(genus_name)) %>%
  distinct(Sample, genus_name) %>%
  count(genus_name, sort = TRUE, name = "n_samples") %>%
  slice_head(n = n_top_genera) %>%
  pull(genus_name)

# For each top genus, find ITS OWN top genes by sample prevalence
top_genes <- args_fox %>%
  filter(genus_name %in% top_genera) %>%
  distinct(Sample, GENE, genus_name) %>%
  count(genus_name, GENE, name = "n_samples") %>%
  group_by(genus_name) %>%
  slice_max(n_samples, n = n_genes_per_genus, with_ties = FALSE) %>%
  ungroup() %>%
  pull(GENE) %>%
  unique()

cat("Selected", length(top_genes), "genes across", n_top_genera, "genera\n")

# --- 2. Build gene x genus matrix: number of SAMPLES with each combination ---
heatmap_data <- args_fox %>%
  filter(GENE %in% top_genes, genus_name %in% top_genera) %>%
  distinct(Sample, GENE, genus_name) %>%
  count(GENE, genus_name, name = "n_samples")

full_grid <- expand.grid(GENE = top_genes, genus_name = top_genera, stringsAsFactors = FALSE)
heatmap_data <- full_grid %>%
  left_join(heatmap_data, by = c("GENE", "genus_name")) %>%
  mutate(n_samples = ifelse(is.na(n_samples), 0, n_samples))

# Order genera by overall prevalence; order genes by which genus they were
# selected for (keeps each genus's genes grouped together visually)
gene_order <- args_fox %>%
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
argsgenera <- ggplot(heatmap_data, aes(x = genus_name, y = GENE, fill = n_samples)) +
  geom_tile(color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_samples > 0, n_samples, "")),
            size = 3, color = ifelse(heatmap_data$n_samples > 5, "white", "black")) +
  scale_fill_gradient(low = "white", high = "#8B3A3A", name = "No. of\nsamples",
                      trans = "sqrt") +
  labs(
    title = "Antimicrobial resistance gene carriage by bacterial genus (fox samples)",
    x = "Genus", y = "Resistance gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 22),
    axis.text.y = element_text(face = "italic", size = 15),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.title.position = "plot",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11)
  )


ggsave("fox_ARG_genus_heatmap_17.pdf", argsgenera, width = 20, height = 22)
print(argsgenera)
cat("\nSaved fox_VG_genus_heatmap.png and .pdf\n")

ggsave("fox_ARG_genus_heatmap.png", argsgenera, width = 16, height = 20, dpi = 400)




## linked contig of VGs, ARGs, plasmids


library(dplyr)

# --- Load the three genus-linked datasets ---
args_fox <- read.csv("args_fox_with_genus.csv", stringsAsFactors = FALSE)
vgs_fox <- read.csv("vgs_fox_with_genus.csv", stringsAsFactors = FALSE)
plasmids_fox <- read.csv("plasmids_fox_with_genus.csv", stringsAsFactors = FALSE)

# --- Standardise to one Sample + contig ID column per dataset ---
# (ARGs file uses "Contig"; VGs and Plasmids files use "SEQUENCE")
args_contigs <- args_fox %>% distinct(Sample, contig_id = Contig)
vgs_contigs <- vgs_fox %>% distinct(Sample, contig_id = SEQUENCE)
plasmid_contigs <- plasmids_fox %>% distinct(Sample, contig_id = SEQUENCE)

cat("Unique contigs with an ARG:", nrow(args_contigs), "\n")
cat("Unique contigs with a VG:", nrow(vgs_contigs), "\n")
cat("Unique contigs with a plasmid marker:", nrow(plasmid_contigs), "\n\n")

# --- Find contigs present in ALL THREE datasets ---
triple_hit_contigs <- args_contigs %>%
  inner_join(vgs_contigs, by = c("Sample", "contig_id")) %>%
  inner_join(plasmid_contigs, by = c("Sample", "contig_id"))

cat("=== Contigs carrying an ARG, a VG, AND a plasmid marker (all three) ===\n")
cat("Total:", nrow(triple_hit_contigs), "\n\n")

if (nrow(triple_hit_contigs) > 0) {
  print(triple_hit_contigs)
  
  # --- Pull the full detail from each dataset for these specific contigs ---
  arg_detail <- args_fox %>%
    inner_join(triple_hit_contigs, by = c("Sample", "Contig" = "contig_id")) %>%
    select(Sample, Contig, GENE, RESISTANCE, genus_name)
  
  vg_detail <- vgs_fox %>%
    inner_join(triple_hit_contigs, by = c("Sample", "SEQUENCE" = "contig_id")) %>%
    select(Sample, SEQUENCE, GENE, PRODUCT, genus_name)
  
  plasmid_detail <- plasmids_fox %>%
    inner_join(triple_hit_contigs, by = c("Sample", "SEQUENCE" = "contig_id")) %>%
    select(Sample, SEQUENCE, GENE, genus_name)
  
  cat("\n=== ARG detail ===\n")
  for (i in seq_len(nrow(arg_detail))) {
    r <- arg_detail[i, ]
    cat(r$Sample, "|", r$Contig, "|", r$genus_name, "| ARG:", r$GENE, "|", r$RESISTANCE, "\n")
  }
  
  cat("\n=== VG detail ===\n")
  for (i in seq_len(nrow(vg_detail))) {
    r <- vg_detail[i, ]
    cat(r$Sample, "|", r$SEQUENCE, "|", r$genus_name, "| VG:", r$GENE, "\n")
  }
  
  cat("\n=== Plasmid marker detail ===\n")
  for (i in seq_len(nrow(plasmid_detail))) {
    r <- plasmid_detail[i, ]
    cat(r$Sample, "|", r$SEQUENCE, "|", r$genus_name, "| Plasmid marker:", r$GENE, "\n")
  }
  
  write.csv(arg_detail, "fox_triple_hit_ARG_detail.csv", row.names = FALSE)
  write.csv(vg_detail, "fox_triple_hit_VG_detail.csv", row.names = FALSE)
  write.csv(plasmid_detail, "fox_triple_hit_plasmid_detail.csv", row.names = FALSE)
  cat("\nSaved fox_triple_hit_ARG_detail.csv, fox_triple_hit_VG_detail.csv,\n")
  cat("and fox_triple_hit_plasmid_detail.csv\n")
} else {
  cat("No contigs were found carrying an ARG, a VG, AND a plasmid marker\n")
  cat("simultaneously in this dataset.\n")
}



## Figure 5.10

library(dplyr)
library(ggplot2)

plasmid_args <- read.csv("fox_plasmid_mediated_ARGs.csv", stringsAsFactors = FALSE)

# --- 1. Count how many contigs link each ARG gene to each plasmid replicon type ---
heatmap_data <- plasmid_args %>%
  distinct(Sample, Contig, GENE, plasmid_replicon_gene) %>%
  count(GENE, plasmid_replicon_gene, name = "n_contigs")

# --- 2. Fill in missing gene-replicon combinations as 0, so the heatmap
#         grid is complete rather than showing gaps ---
all_genes <- sort(unique(heatmap_data$GENE))
all_replicons <- sort(unique(heatmap_data$plasmid_replicon_gene))

full_grid <- expand.grid(GENE = all_genes, plasmid_replicon_gene = all_replicons,
                         stringsAsFactors = FALSE)

heatmap_data <- full_grid %>%
  left_join(heatmap_data, by = c("GENE", "plasmid_replicon_gene")) %>%
  mutate(n_contigs = ifelse(is.na(n_contigs), 0, n_contigs))

# --- 3. Order genes and replicon types by overall prevalence, for a
#         cleaner-reading plot ---
gene_order <- heatmap_data %>%
  group_by(GENE) %>%
  summarise(total = sum(n_contigs), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(GENE)

replicon_order <- heatmap_data %>%
  group_by(plasmid_replicon_gene) %>%
  summarise(total = sum(n_contigs), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(plasmid_replicon_gene)

heatmap_data$GENE <- factor(heatmap_data$GENE, levels = rev(gene_order))
heatmap_data$plasmid_replicon_gene <- factor(heatmap_data$plasmid_replicon_gene, levels = replicon_order)

# --- 4. Plot ---
p <- ggplot(heatmap_data, aes(x = plasmid_replicon_gene, y = GENE, fill = n_contigs)) +
  geom_tile(color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_contigs > 0, n_contigs, "")),
            size = 3.2, color = ifelse(heatmap_data$n_contigs > 2, "white", "black")) +
  scale_fill_gradient(low = "white", high = "darkred", name = "No. of\ncontigs") +
  labs(
    title = "Resistance genes associated with plasmid replicon types",
    x = "Plasmid replicon type", y = "Resistance gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10, face = "italic"),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.title.position = "plot",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

ggsave("fox_ARG_plasmid_heatmap.png", p, width = 10, height = 8, dpi = 400)
ggsave("fox_ARG_plasmid_heatmap.pdf", p, width = 10, height = 8)
print(p)
cat("Saved fox_ARG_plasmid_heatmap.png and .pdf\n")


## -----------------------------------------------------------------------------


## Figure 5.11

library(ggalluvial)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Prepare data: one row per ARG hit, with genus, resistance class,
#         and plasmid replicon type
plasmid_args_fox <- read.csv("fox_plasmid_mediated_ARGs.csv", stringsAsFactors = FALSE)

alluvial_data <- plasmid_args_fox %>%
  filter(!is.na(genus_name)) %>%
  # RESISTANCE can list multiple classes separated by ";" - take the first
  # for simplicity in this lightweight visualisation
  mutate(resistance_class = sub(";.*", "", RESISTANCE)) %>%
  mutate(resistance_class = paste0(toupper(substr(resistance_class, 1, 1)),
                                   substr(resistance_class, 2, nchar(resistance_class)))) %>%
  rename(plasmid_type = plasmid_replicon_gene) %>%
  filter(!is.na(plasmid_type), plasmid_type != "NA", plasmid_type != "") %>%
  count(genus_name, resistance_class, plasmid_type, name = "freq")

cat("Alluvial data built:", nrow(alluvial_data), "unique genus-resistance-plasmid combinations\n")
print(alluvial_data)

# 2. Plot
n_genera <- n_distinct(alluvial_data$genus_name)

if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(RColorBrewer)
genus_palette <- colorRampPalette(brewer.pal(12, "Paired"))(n_genera)
names(genus_palette) <- sort(unique(alluvial_data$genus_name))

p <- ggplot(alluvial_data,
            aes(axis1 = genus_name, axis2 = resistance_class, axis3 = plasmid_type,
                y = freq)) +
  geom_alluvium(aes(fill = genus_name), width = 0.28, alpha = 0.85, color = "grey20", linewidth = 0.1) +
  geom_stratum(width = 0.28, fill = "grey92", color = "grey20", linewidth = 0.4) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5.5, fontface = "bold") +
  scale_fill_manual(values = genus_palette) +
  scale_x_discrete(limits = c("Genus", "Resistance class", "Plasmid replicon type"),
                   expand = c(0.12, 0.12)) +
  labs(
    title = "Plasmid-mediated ARGs (fox dataset): genus, resistance class, and plasmid type",
    x = NULL, y = "Number of ARG hits"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave("fox_ARG_plasmid_alluvial.png", p, width = 12, height = 10, dpi = 300)
ggsave("fox_ARG_plasmid_alluvial.pdf", p, width = 30, height = 15)
print(p)





## -----------------------------------------------------------------------------

## Urban versus Rural ARGs comparisons


urban <- read.csv("ARGs_on_reads.csv", check.names = FALSE) %>%
  mutate(Sample = as.character(Sample), Habitat = "Urban")

rural <- read.csv("Rural_foxes_ARGs.csv", check.names = FALSE) %>%
  mutate(Sample = as.character(Sample), Habitat = "Rural")

args_all <- bind_rows(urban, rural)

# Confirm sample counts match what's expected
args_all %>% distinct(Sample, Habitat) %>% count(Habitat)

# Build presence/absence matrix
gene_mat_all <- args_all %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Sample, Gene) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Gene, values_from = present, values_fill = 0) %>%
  column_to_rownames("Sample")

dim(gene_mat_all)  # samples x distinct genes, combined across both habitats

# Habitat lookup aligned to matrix row order
sample_habitat <- args_all %>%
  distinct(Sample, Habitat) %>%
  filter(Sample %in% rownames(gene_mat_all)) %>%
  arrange(match(Sample, rownames(gene_mat_all)))

stopifnot(all(sample_habitat$Sample == rownames(gene_mat_all)))


# Descriptive comparison

# Gene richness per sample
gene_richness_all <- args_all %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Sample, Gene) %>%
  count(Sample, name = "n_genes") %>%
  left_join(sample_habitat, by = "Sample")

gene_richness_all %>%
  group_by(Habitat) %>%
  summarise(
    n_samples  = n(),
    mean_genes = mean(n_genes),
    sd_genes   = sd(n_genes),
    min_genes  = min(n_genes),
    max_genes  = max(n_genes)
  )

# Which genes are shared vs unique to each habitat (core/accessory-style comparison)
genes_by_habitat <- args_all %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Gene, Habitat)

urban_genes <- genes_by_habitat %>% filter(Habitat == "Urban") %>% pull(Gene)
rural_genes <- genes_by_habitat %>% filter(Habitat == "Rural") %>% pull(Gene)

shared_genes     <- intersect(urban_genes, rural_genes)
urban_only_genes <- setdiff(urban_genes, rural_genes)
rural_only_genes <- setdiff(rural_genes, urban_genes)

length(shared_genes)      # genes found in both habitats
length(urban_only_genes)  # genes found only in urban foxes
length(rural_only_genes)  # genes found only in rural foxes


# Alpha diversity - is richness different between habitats?

wilcox.test(n_genes ~ Habitat, data = gene_richness_all)


# Beta diversity - does gene composition differ between habitats?

gene_jaccard_all <- vegdist(gene_mat_all, method = "jaccard", binary = TRUE)

bd_habitat <- betadisper(gene_jaccard_all, sample_habitat$Habitat)
anova(bd_habitat)

adonis_habitat <- adonis2(gene_mat_all ~ Habitat, data = sample_habitat,
                          method = "jaccard", binary = TRUE, permutations = 9999)
adonis_habitat


# 1. Simple lists: genes unique to each habitat
sort(urban_only_genes)
sort(rural_only_genes)

# 2. For SHARED genes, how common are they in each habitat?

gene_prevalence <- args_all %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Sample, Gene, Habitat) %>%
  count(Gene, Habitat, name = "n_samples_present") %>%
  left_join(
    sample_habitat %>% count(Habitat, name = "n_samples_total"),
    by = "Habitat"
  ) %>%
  mutate(prevalence = n_samples_present / n_samples_total) %>%
  select(Gene, Habitat, prevalence) %>%
  pivot_wider(names_from = Habitat, values_from = prevalence, values_fill = 0) %>%
  mutate(prevalence_diff = Urban - Rural) %>%
  arrange(desc(abs(prevalence_diff)))

gene_prevalence

# ---- 3. Formal per-gene test: Fisher's exact test for each gene ----

gene_presence_long <- args_all %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Sample, Gene, Habitat)

all_genes <- unique(gene_presence_long$Gene)

fisher_results <- map_dfr(all_genes, function(g) {
  present_samples <- gene_presence_long %>% filter(Gene == g) %>% pull(Sample)
  tab <- sample_habitat %>%
    mutate(Present = Sample %in% present_samples) %>%
    count(Habitat, Present) %>%
    pivot_wider(names_from = Present, values_from = n, values_fill = 0)
  # build 2x2 table safely even if a cell is entirely 0
  m <- matrix(c(
    tab$`TRUE`[tab$Habitat == "Urban"],  tab$`FALSE`[tab$Habitat == "Urban"],
    tab$`TRUE`[tab$Habitat == "Rural"],  tab$`FALSE`[tab$Habitat == "Rural"]
  ), nrow = 2, byrow = TRUE)
  p <- fisher.test(m)$p.value
  tibble(Gene = g, p_value = p)
})

fisher_results <- fisher_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

fisher_results  # top rows = genes most strongly associated with habitat
sum(fisher_results$p_adj < 0.05)  # how many survive multiple-testing correction




## Class-level summary


# ---- Load the combined file (already has Sample, Gene, Class, Cohort) ----
args_classes_all <- read.csv("urban_rural_classes.csv", check.names = FALSE)  # adjust filename

# Quick check before trusting it
args_classes_all %>% count(Cohort)
sum(is.na(args_classes_all$Class))  # should be 0 now

# ---- Rebuild gene lists per habitat from this file directly ----
urban_genes <- args_classes_all %>% filter(Cohort == "Urban") %>% distinct(Gene) %>% pull(Gene)
rural_genes <- args_classes_all %>% filter(Cohort == "Rural") %>% distinct(Gene) %>% pull(Gene)

shared_genes     <- intersect(urban_genes, rural_genes)
urban_only_genes <- setdiff(urban_genes, rural_genes)
rural_only_genes <- setdiff(rural_genes, urban_genes)

length(shared_genes); length(urban_only_genes); length(rural_only_genes)  # sanity check: should still be 55/105/48

# ---- Class-level summary table, now with full Class coverage ----
overlap_table <- tibble(
  Gene = c(shared_genes, urban_only_genes, rural_only_genes),
  Status = c(rep("Shared", length(shared_genes)),
             rep("Urban only", length(urban_only_genes)),
             rep("Rural only", length(rural_only_genes)))
) %>%
  left_join(args_classes_all %>% distinct(Gene, Class), by = "Gene")

sum(is.na(overlap_table$Class))  # should be 0 - flag if not

class_summary <- overlap_table %>%
  count(Class, Status) %>%
  pivot_wider(names_from = Status, values_from = n, values_fill = 0) %>%
  mutate(Total = Shared + `Urban only` + `Rural only`) %>%
  arrange(desc(Total))

class_summary

write.csv(class_summary, "ARG_class_summary_urban_rural.csv", row.names = FALSE)



args_classes_all <- read.csv("urban_rural_classes.csv", check.names = FALSE) %>%
  mutate(Gene = str_trim(Gene))  # trim whitespace defensively, same as before

args_classes_all %>% count(Cohort)  # should now show both Urban and Rural

top_genes <- gene_prevalence %>%
  mutate(Gene = str_trim(Gene)) %>%
  arrange(desc(abs(prevalence_diff))) %>%
  slice_head(n = 40) %>%
  left_join(args_classes_all %>% distinct(Gene, Class), by = "Gene") %>%
  mutate(Gene = fct_reorder(Gene, prevalence_diff))

sum(is.na(top_genes$Class))  # should now be 0


## Figure 5.12

class_colours <- c(
  "Tetracyclines"               = "#8DD3C7",
  "Aminoglycosides"              = "#E6B800",
  "Lincosamides"                  = "#BEBADA",
  "Phenicols"                     = "#FB8072",
  "Sulfonamides"                  = "#80B1D3",
  "Trimethoprim"                  = "#FDB462",
  "beta-lactams"                  = "#B3DE69",
  "MLSB"                          = "#FCCDE5",
  "Macrolides"                    = "#D9D9D9",
  "Oxazolidinones"                = "#BC80BD",
  "MDR"                           = "#CCEBC5",
  "Glycopeptides"                 = "#FFED6F",
  "Fosfomycin"                    = "#E41A1C",
  "Ionophore resistance"          = "#377EB8",
  "Streptogramin (A)"             = "#4DAF4A",
  "Fluoroquinolones/Quinolones"   = "gold"
)


p_prevalence <- ggplot(top_genes, aes(x = prevalence_diff, y = Gene, colour = Class)) +
  geom_segment(aes(x = 0, xend = prevalence_diff, y = Gene, yend = Gene), linewidth = 1.2) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  scale_colour_manual(values = class_colours, drop = TRUE) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0.5, face = "italic", colour = "grey40", size = 10)
  ) +
  labs(
    title = "Top 40 genes with the largest prevalence difference between habitats",
    x = "Prevalence difference",
    y = NULL,
    colour = "Class",
    caption = "Negative values: more common in Rural foxes   |   Positive values: more common in Urban foxes"
  )

p_prevalence
ggsave("top40_genes_prevalence_diff.pdf", plot = p_prevalence, width = 10, height = 10.5, units = "in")



## -----------------------------------------------------------------------------

## COMPARE URBAN AND RURAL MICROBIOMES

library(dplyr)
library(tidyr)

# --- 1. Load both merged phylum files ---
urban <- read.csv("bracken_phylum_merged.csv", stringsAsFactors = FALSE, check.names = FALSE)
rural <- read.csv("bracken_phylum_merged_rural.csv", stringsAsFactors = FALSE, check.names = FALSE)

# --- 2. Identify sample columns automatically (anything not name/taxonomy) ---
urban_sample_cols <- setdiff(colnames(urban), c("name", "taxonomy_id", "taxonomy_lvl"))
rural_sample_cols <- setdiff(colnames(rural), c("name", "taxonomy_id", "taxonomy_lvl"))

cat("Urban samples (n =", length(urban_sample_cols), "):\n")
print(urban_sample_cols)
cat("\nRural samples (n =", length(rural_sample_cols), "):\n")
print(rural_sample_cols)

# --- 3. Exclude Chordata (host DNA), consistent with earlier analyses ---
urban <- urban %>% filter(name != "Chordata")
rural <- rural %>% filter(name != "Chordata")

# --- 4. Reshape both to long format, add a group label ---
urban_long <- urban %>%
  select(name, all_of(urban_sample_cols)) %>%
  pivot_longer(all_of(urban_sample_cols), names_to = "sample", values_to = "reads") %>%
  mutate(group = "Urban")

rural_long <- rural %>%
  select(name, all_of(rural_sample_cols)) %>%
  pivot_longer(all_of(rural_sample_cols), names_to = "sample", values_to = "reads") %>%
  mutate(group = "Rural")

combined <- bind_rows(urban_long, rural_long) %>%
  mutate(present = reads > 0)

n_urban <- length(urban_sample_cols)
n_rural <- length(rural_sample_cols)

# --- 5. For each phylum, build a 2x2 table and run Fisher's exact test ---
phyla <- unique(combined$name)

results <- data.frame(
  phylum = character(), urban_present = integer(), urban_total = integer(),
  rural_present = integer(), rural_total = integer(),
  p_value = numeric(), stringsAsFactors = FALSE
)

for (ph in phyla) {
  sub <- combined %>% filter(name == ph)
  
  urban_present <- sum(sub$group == "Urban" & sub$present)
  rural_present <- sum(sub$group == "Rural" & sub$present)
  
  # 2x2 contingency table: rows = present/absent, columns = urban/rural
  tab <- matrix(
    c(urban_present, n_urban - urban_present,
      rural_present, n_rural - rural_present),
    nrow = 2,
    dimnames = list(c("Present", "Absent"), c("Urban", "Rural"))
  )
  
  test <- fisher.test(tab)
  
  results <- rbind(results, data.frame(
    phylum = ph,
    urban_present = urban_present, urban_total = n_urban,
    rural_present = rural_present, rural_total = n_rural,
    p_value = round(test$p.value, 4)
  ))
}

results <- results %>% arrange(p_value)

cat("\n=== Fisher's exact test results: phylum presence/absence, urban vs rural foxes ===\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(r$phylum, ": urban", r$urban_present, "/", r$urban_total,
      "| rural", r$rural_present, "/", r$rural_total,
      "| p =", r$p_value, "\n")
}

write.csv(results, "fox_phylum_presence_absence_fisher.csv", row.names = FALSE)
cat("\nSaved to fox_phylum_presence_absence_fisher.csv\n")

cat("\nPhyla with p < 0.05 (before any multiple-testing correction):\n")
print(results %>% filter(p_value < 0.05))

results$p_adjusted_BH <- p.adjust(results$p_value, method = "BH")

cat("Phyla with adjusted p < 0.05:\n")
print(results %>% filter(p_adjusted_BH < 0.05))










#### OTHER SUPPLEMENTARY FIGURES ####

## Figure S5.1

# strip commas and convert to numeric
metadata <- metadata %>%
  mutate(Reads = as.numeric(gsub(",", "", Reads)))

# One row per gene detected per sample here, so counting genes per sample
# is just counting rows per Sample
genes_per_sample <- args %>%
  count(Sample, name = "n_genes")

gene_depth <- metadata %>%
  left_join(genes_per_sample, by = "Sample")

# Quick look at the numbers themselves before plotting
gene_depth %>% select(Sample, Reads, n_genes) %>% arrange(Reads)

# Any sample in metadata with no matching rows in ARGs_on_reads.csv will
# show NA for n_genes here - worth checking rather than assuming 0 genes
sum(is.na(gene_depth$n_genes))

ggplot(gene_depth, aes(Reads, n_genes)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = Sample), size = 3) +
  geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
  theme_bw(base_size = 13) +
  labs(
    title = "Genes detected vs sequencing depth",
    x = "Total reads",
    y = "Number of genes detected"
  )

cor.test(gene_depth$Reads, gene_depth$n_genes, method = "spearman")



## Figure S5.2

library(tidyverse)

# ---- Gene richness per sample ----
# One row per gene hit per sample in args - count distinct genes per sample
gene_richness <- args %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct(Sample, Gene) %>%
  count(Sample, name = "n_genes") %>%
  arrange(desc(n_genes))

gene_richness

# ---- Add Location for context ----
gene_richness <- gene_richness %>%
  left_join(metadata %>% select(Sample, Location), by = "Sample")

gene_richness

# ---- Summary stats ----
summary(gene_richness$n_genes)

gene_richness %>%
  group_by(Location) %>%
  summarise(
    n_samples = n(),
    mean_genes = mean(n_genes),
    sd_genes = sd(n_genes),
    min_genes = min(n_genes),
    max_genes = max(n_genes)
  )

# ---- Quick bar plot ----
ggplot(gene_richness, aes(x = reorder(Sample, -n_genes), y = n_genes, fill = Location)) +
  geom_col() +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(
    title = "ARG richness per sample",
    x = "Sample",
    y = "Number of distinct ARGs detected"
  )



# ---- Wilcoxon rank-sum test: gene richness by Location ----
wilcox.test(n_genes ~ Location, data = gene_richness)

# Exact test explicitly, in case of ties triggering a normal approximation warning
wilcox.test(n_genes ~ Location, data = gene_richness, exact = TRUE)


################################################################################

## Figure S5.6

vgs_raw <- read.csv("VGs_on_contigs.csv", stringsAsFactors = FALSE)

# --- 1. Select top 40 genes by sample prevalence ---
n_top <- 40

top_genes <- vgs_raw %>%
  distinct(Sample, GENE) %>%
  count(GENE, sort = TRUE, name = "n_samples") %>%
  slice_head(n = n_top) %>%
  pull(GENE)

# --- 2. Build gene x sample matrix: number of contig hits per combination ---
heatmap_data <- vgs_raw %>%
  filter(GENE %in% top_genes) %>%
  count(Sample, GENE, name = "n_hits")

# Fill in missing combinations as 0 (gene not detected in that sample)
all_samples <- unique(vgs_raw$Sample)
full_grid <- expand.grid(Sample = all_samples, GENE = top_genes, stringsAsFactors = FALSE)
heatmap_data <- full_grid %>%
  left_join(heatmap_data, by = c("Sample", "GENE")) %>%
  mutate(n_hits = ifelse(is.na(n_hits), 0, n_hits))

# Order genes by prevalence (most common at top), samples alphabetically
heatmap_data$GENE <- factor(heatmap_data$GENE, levels = rev(top_genes))
heatmap_data$Sample <- factor(heatmap_data$Sample, levels = sort(all_samples))

# --- 3. Plot ---
p <- ggplot(heatmap_data, aes(x = Sample, y = GENE, fill = n_hits)) +
  geom_tile(color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_hits > 0, n_hits, "")),
            size = 3.2, color = ifelse(heatmap_data$n_hits > 3, "white", "black")) +
  scale_fill_gradient(low = "white", high = "#1B4D3E", name = "No. of\ncontig hits") +
  labs(
    title = "Most prevalent virulence genes across fox samples",
    x = "Sample", y = "Virulence gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(face = "italic", size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.title.position = "plot",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

ggsave("fox_top_VGs_heatmap.png", p, width = 8, height = 12, dpi = 400)
ggsave("fox_top_VGs_heatmap.pdf", p, width = 9, height = 10)
print(p)
cat("Saved fox_top_VGs_heatmap.png and .pdf\n")





# Top 40 most prevalent resistance genes across
# the 7 fox samples, showing number of contig hits per gene per sample.

## Figure S5.7

library(dplyr)
library(tidyr)
library(ggplot2)

args_raw <- read.csv("ARGs_on_contigs.csv", stringsAsFactors = FALSE)

# --- 1. Select top 40 genes by sample prevalence ---
n_top <- 40

top_genes <- args_raw %>%
  distinct(Sample, GENE) %>%
  count(GENE, sort = TRUE, name = "n_samples") %>%
  slice_head(n = n_top) %>%
  pull(GENE)

# --- 2. Build gene x sample matrix: number of contig hits per combination ---
heatmap_data <- args_raw %>%
  filter(GENE %in% top_genes) %>%
  count(Sample, GENE, name = "n_hits")

# Fill in missing combinations as 0 (gene not detected in that sample)
all_samples <- unique(args_raw$Sample)
full_grid <- expand.grid(Sample = all_samples, GENE = top_genes, stringsAsFactors = FALSE)
heatmap_data <- full_grid %>%
  left_join(heatmap_data, by = c("Sample", "GENE")) %>%
  mutate(n_hits = ifelse(is.na(n_hits), 0, n_hits))

# Order genes by prevalence (most common at top), samples alphabetically
heatmap_data$GENE <- factor(heatmap_data$GENE, levels = rev(top_genes))
heatmap_data$Sample <- factor(heatmap_data$Sample, levels = sort(all_samples))

# --- 3. Plot ---
commonargs <- ggplot(heatmap_data, aes(x = Sample, y = GENE, fill = n_hits)) +
  geom_tile(color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_hits > 0, n_hits, "")),
            size = 3.2, color = ifelse(heatmap_data$n_hits > 10, "white", "black")) +
  scale_fill_gradient(low = "white", high = "#8B3A3A", name = "No. of\ncontig hits") +
  labs(
    title = "Most prevalent genes across fox samples",
    x = "Sample", y = "Antimicrobial resistance gene"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(face = "italic", size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.title.position = "plot",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

ggsave("fox_top_ARGs_heatmap.png", commonargs, width = 8, height = 12, dpi = 400)
ggsave("fox_top_ARGs_heatmap.pdf", commonargs, width = 9, height = 10)
print(commonargs)
cat("Saved fox_top_VGs_heatmap.png and .pdf\n")









################################################################################
################################################################################



















































## Statistical analyses

library(vegan)

# ---- PERMANOVA on gene presence/absence (Jaccard distance) ----
# args_urban_mat: samples as rows, genes as columns (0/1) - already built earlier
# metadata: already loaded (fox_metadata_seq.csv)

sample_meta <- metadata %>%
  filter(Sample %in% rownames(args_urban_mat)) %>%
  arrange(match(Sample, rownames(args_urban_mat)))  # match row order to args_urban_mat

stopifnot(all(sample_meta$Sample == rownames(args_urban_mat)))  # confirm alignment before testing

# Check dispersion homogeneity first, same logic as the VITEK analysis -
# PERMANOVA assumes similar within-group spread between Location groups
gene_jaccard <- vegdist(args_urban_mat, method = "jaccard", binary = TRUE)

bd_gene_location <- betadisper(gene_jaccard, sample_meta$Location)
anova(bd_gene_location)

# PERMANOVA
adonis_gene_location <- adonis2(args_urban_mat ~ Location, data = sample_meta,
                                method = "jaccard", binary = TRUE, permutations = 9999)
adonis_gene_location



######################################################################









