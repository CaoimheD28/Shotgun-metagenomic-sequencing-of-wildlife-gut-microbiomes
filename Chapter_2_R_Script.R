setwd("C:/Users/cedoy/ScopingReview")
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


## -----------------------------------------------------------------------------
## Figure 2.2 Species counts, wild vs. captive


species_counts <- read.csv("species_counts.csv") %>%
  filter(Species != "" & !is.na(Species)) %>%
  replace_na(list(Wild = 0, Captive = 0)) %>%
  mutate(Total = Wild + Captive) %>%
  pivot_longer(
    cols = c(Wild, Captive),
    names_to = "Type",
    values_to = "Count"
  ) %>%
  mutate(Species = reorder(Species, -Total))



p1 <- ggplot(species_counts, aes(x = Species, y = Count, fill = Type)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = c("Wild" = "#2166AC", "Captive" = "#E8834D")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Species", y = "% of studies", fill = "Animal type") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.position = "top",
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  )

p1


ggsave(
  filename = "UPDATED_totalspeciescount.png",
  plot = p1,
  width = 24,
  height = 10,
  dpi = 300
)



## Figure 2.3 Regional heatmap by species and study context

region_species <- read.csv("region_spp.csv", check.names = FALSE)

region_counts <- region_species %>%
  group_by(Region, Species, Context) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(Region, Species, Context, fill = list(n = 0)) %>%
  mutate(Context = factor(Context, levels = c("Wild", "Captive")))

p2 <- ggplot(region_counts, aes(x = Region, y = Species, fill = n)) +
  geom_tile(color = "white", linewidth = 0.6) +
  facet_wrap(~ Context) +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkred",
    name = "Study count",
    guide = guide_colorbar(
      barwidth = 1,
      barheight = 8,
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  labs(x = "Region", y = "Species") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 15, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.spacing = unit(1, "lines"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )


p2

ggsave(
  filename = "UPDATED_geographic.png",
  plot = p2,
  width = 24,
  height = 14,
  dpi = 300
)

## ---------------------------------------------------------------------------
## Figure 2.4a Sample method x sample type, wild

df <- read.csv("Sample_method_type.csv")

df_clean <- df %>%
  separate_rows(Sample_method, sep = ",") %>%
  separate_rows(Sample_type, sep = ",") %>%
  mutate(
    Sample_method = trimws(Sample_method),
    Sample_type = trimws(Sample_type)
  )

p4 <- ggplot(df_clean, aes(x = Sample_method, fill = Sample_type)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2, width = 0.7) +
  scale_fill_viridis_d(option = "viridis") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sampling method", y = "Proportion of studies", fill = "Sample type") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 18, angle = 65, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

## -----------------------------------------------------------------------------
## Figure 2.4b Sample method x sample type, captive


df_captive <- read.csv("Captive_Sample_Method.csv")

df_captive_clean <- df_captive %>%
  separate_rows(Study_type, sep = ",") %>%
  separate_rows(Sample_type, sep = ",") %>%
  mutate(
    Study_type = trimws(Study_type),
    Sample_type = trimws(Sample_type)
  )

p5 <- ggplot(df_captive_clean, aes(x = Study_type, fill = Sample_type)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2, width = 0.7) +
  scale_fill_viridis_d(option = "viridis") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Study type", y = "Proportion of studies", fill = "Sample type") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 18, angle = 65, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

## ============================================================
## Combine with patchwork
## ============================================================

combined_45 <- p4 / p5 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

combined_45

ggsave(
  filename = "sample_type_method_combined.png",
  plot = combined_45,
  width = 20,
  height = 14,
  dpi = 300
)


## -----------------------------------------------------------------------------
## Figure 2.5: Microbial groups investigated, wild vs. captive


library(ggplot2)
library(dplyr)

df <- data.frame(
  Animal_Context = rep(c("Wild", "Captive"), each = 8),
  Category = rep(c("Enterobacterales", "Other Gram-negative", "Gram-positive",
                   "Anaerobes", "Mycobacteria", "Fungi", "Parasites", "N/A Molecular"), 2),
  Count = c(50.7, 8.5, 38, 0, 4.2, 1.4, 0, 11.3,
            48.2, 26.8, 39.3, 3.6, 1.8, 1.8, 5.4, 1.8)
)

# Reorder categories by total count (largest first, so they appear on the left)
category_order <- df %>%
  group_by(Category) %>%
  summarise(Total = sum(Count)) %>%
  arrange(desc(Total)) %>%
  pull(Category)

df$Category <- factor(df$Category, levels = category_order)

# Muted, professional colour palette
category_colors <- c(
  "Enterobacterales"    = "#A07C7B",
  "Other Gram-negative" = "#6B9FA8",
  "Gram-positive"       = "#8E9C87",
  "Anaerobes"           = "#D9A761",
  "Mycobacteria"        = "#A6D96A",
  "Fungi"               = "#D87B6B",
  "Parasites"           = "#7F6B9E",
  "N/A Molecular"       = "#66C2A5"
)

p6 <- ggplot(df %>% filter(Count > 0),
             aes(x = Category, y = Animal_Context, size = Count, fill = Category)) +
  geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.9) +
  scale_fill_manual(values = category_colors, guide = "none") +
  scale_size_continuous(range = c(4, 22), breaks = c(5, 10, 20, 30, 50)) +
  labs(
    x = "Microbial category",
    y = "Animal cohort",
    size = "% studies"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 13, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

p6

ggsave(
  filename = "microbial_groups_bubble.png",
  plot = p6,
  width = 11,
  height = 6,
  dpi = 300
)


## -----------------------------------------------------------------------------
## Figure 2.6a Detection/isolation methods, wild vs. captive

df <- read.csv("ID_ISO.csv", check.names = FALSE)

df_clean <- df %>%
  filter(`Detection_isolation method` != "" & !is.na(`Detection_isolation method`)) %>%
  replace_na(list(Wild = 0, Captive = 0)) %>%
  mutate(Total = Wild + Captive) %>%
  pivot_longer(
    cols = c(Wild, Captive),
    names_to = "Type",
    values_to = "Count"
  ) %>%
  mutate(`Detection_isolation method` = reorder(`Detection_isolation method`, -Total))

p7 <- ggplot(df_clean, aes(x = `Detection_isolation method`, y = Count, fill = Type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65,
           color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Wild" = "#2166AC", "Captive" = "#E8834D")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Detection / isolation method", y = "% of studies", fill = "Animal context") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )

## -----------------------------------------------------------------------------
## Figure 2.6b AMR detection methods, wild vs. captive


df2 <- read.csv("amr_test.csv", check.names = FALSE)

df2_clean <- df2 %>%
  filter(Detection_method != "" & !is.na(Detection_method)) %>%
  replace_na(list(Wild = 0, Captive = 0)) %>%
  mutate(Total = Wild + Captive) %>%
  pivot_longer(
    cols = c(Wild, Captive),
    names_to = "Type",
    values_to = "Count"
  ) %>%
  mutate(Detection_method = reorder(Detection_method, -Total))

p8 <- ggplot(df2_clean, aes(x = Detection_method, y = Count, fill = Type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65,
           color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Wild" = "#2166AC", "Captive" = "#E8834D")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Detection method", y = "% of studies", fill = "Animal context") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )


## Combine

combined_78 <- p7 / p8 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

combined_78

ggsave(
  filename = "detection_methods_combined.png",
  plot = combined_78,
  width = 22,
  height = 15,
  dpi = 300
)


## -----------------------------------------------------------------------------

## Figure 2.7a Sankey/alluvial plot: Microbial group -> Isolation method -> AMR method (wild)

library(ggalluvial)
library(RColorBrewer)


df <- read.csv("wild_microbe_amrtest.csv", check.names = FALSE)

df <- df[, !(is.na(names(df)) | names(df) == "")]

df_long <- df %>%
  filter(
    !is.na(Species), Species != "",
    !is.na(Isolation_Method), Isolation_Method != "",
    !is.na(AMR_Method), AMR_Method != ""
  ) %>%
  separate_rows(Species, sep = ";") %>%
  separate_rows(Isolation_Method, sep = ";") %>%
  separate_rows(AMR_Method, sep = ";") %>%
  mutate(
    Species = trimws(Species),
    Isolation_Method = trimws(Isolation_Method),
    AMR_Method = trimws(AMR_Method)
  ) %>%
  filter(
    Species != "", Species != "N/A",
    Isolation_Method != "", Isolation_Method != "N/A",
    AMR_Method != "", AMR_Method != "N/A"
  )


## Count flows and validate alluvial structure


alluvial_counts <- df_long %>%
  group_by(Species, Isolation_Method, AMR_Method) %>%
  summarise(Freq = n(), .groups = "drop")

stopifnot(is_alluvia_form(alluvial_counts, axes = 1:3, silent = TRUE))


## Order microbial groups by total frequency, largest first, so the
## biggest flows anchor the top of the diagram and are easiest to trace


group_order <- alluvial_counts %>%
  group_by(Species) %>%
  summarise(Total = sum(Freq)) %>%
  arrange(desc(Total)) %>%
  pull(Species)

alluvial_counts <- alluvial_counts %>%
  mutate(Species = factor(Species, levels = group_order))


## Qualitative, high-contrast palette designed for categorical distinctness


n_groups <- length(group_order)

if (n_groups <= 8) {
  base_colors <- brewer.pal(max(n_groups, 3), "Dark2")[1:n_groups]
} else {
  # extend beyond Dark2's 8 colours only when actually needed
  base_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_groups)
}

group_colors <- setNames(base_colors, group_order)


## Plot


library(ggrepel)

p_sankey <- ggplot(alluvial_counts,
                   aes(axis1 = Species,
                       axis2 = Isolation_Method,
                       axis3 = AMR_Method,
                       y = Freq)) +
  
  geom_alluvium(aes(fill = Species),
                width = 1/12,
                alpha = 0.8,
                color = "white",
                linewidth = 0.15) +
  
  geom_stratum(width = 1/8,
               fill = "grey92",
               color = "black",
               linewidth = 0.4) +
  
  
  geom_label_repel(stat = "stratum",
                   aes(label = after_stat(stratum)),
                   size = 4.6,
                   fontface = "bold",
                   label.size = 0,
                   label.padding = unit(0.2, "lines"),
                   fill = "white",
                   alpha = 0.95,
                   direction = "y",
                   nudge_x = 0.15,
                   box.padding = 0.3,
                   segment.size = 0.3,
                   segment.color = "grey40",
                   max.overlaps = Inf,
                   seed = 42) +
  
  scale_fill_manual(values = group_colors, name = "Microbial group") +
  
  scale_x_discrete(
    limits = c("Microbial group", "Isolation method", "AMR detection method"),
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  
  labs(y = "Number of studies") +
  
  
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

p_sankey

ggsave(
  "wild_three_layer_sankey.jpeg",
  plot = p_sankey,
  width = 28,
  height = 15,
  dpi = 300
)

## -----------------------------------------------------------------------------
## Figure 2.7b Sankey/alluvial plot: Microbial group -> Isolation method -> AMR method (captive)


df <- read.csv("captive_microbe_AMRtest.csv", check.names = FALSE)

df <- df[, !(is.na(names(df)) | names(df) == "")]

df_long <- df %>%
  filter(
    !is.na(Species), Species != "",
    !is.na(Isolation_Method), Isolation_Method != "",
    !is.na(AMR_Method), AMR_Method != ""
  ) %>%
  separate_rows(Species, sep = ";") %>%
  separate_rows(Isolation_Method, sep = ";") %>%
  separate_rows(AMR_Method, sep = ";") %>%
  mutate(
    Species = trimws(Species),
    Isolation_Method = trimws(Isolation_Method),
    AMR_Method = trimws(AMR_Method)
  ) %>%
  filter(
    Species != "", Species != "N/A",
    Isolation_Method != "", Isolation_Method != "N/A",
    AMR_Method != "", AMR_Method != "N/A"
  )


## Count flows and validate alluvial structure

alluvial_counts <- df_long %>%
  group_by(Species, Isolation_Method, AMR_Method) %>%
  summarise(Freq = n(), .groups = "drop")

stopifnot(is_alluvia_form(alluvial_counts, axes = 1:3, silent = TRUE))


## Order microbial groups by total frequency, largest first, so the
## biggest flows anchor the top of the diagram and are easiest to trace

group_order <- alluvial_counts %>%
  group_by(Species) %>%
  summarise(Total = sum(Freq)) %>%
  arrange(desc(Total)) %>%
  pull(Species)

alluvial_counts <- alluvial_counts %>%
  mutate(Species = factor(Species, levels = group_order))

## Qualitative, high-contrast palette designed for categorical distinctness

n_groups <- length(group_order)

if (n_groups <= 8) {
  base_colors <- brewer.pal(max(n_groups, 3), "Dark2")[1:n_groups]
} else {
  # extend beyond Dark2's 8 colours only when actually needed
  base_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_groups)
}

group_colors <- setNames(base_colors, group_order)

## Plot

p_sankey2 <- ggplot(alluvial_counts,
                    aes(axis1 = Species,
                        axis2 = Isolation_Method,
                        axis3 = AMR_Method,
                        y = Freq)) +
  
  geom_alluvium(aes(fill = Species),
                width = 1/12,
                alpha = 0.8,
                color = "white",
                linewidth = 0.15) +
  
  geom_stratum(width = 1/8,
               fill = "grey92",
               color = "black",
               linewidth = 0.4) +
  
  geom_label_repel(stat = "stratum",
                   aes(label = after_stat(stratum)),
                   size = 4.6,
                   fontface = "bold",
                   label.size = 0,
                   label.padding = unit(0.2, "lines"),
                   fill = "white",
                   alpha = 0.95,
                   direction = "y",
                   nudge_x = 0.15,
                   box.padding = 0.3,
                   segment.size = 0.3,
                   segment.color = "grey40",
                   max.overlaps = Inf,
                   seed = 42) +
  
  scale_fill_manual(values = group_colors, name = "Microbial group") +
  
  scale_x_discrete(
    limits = c("Microbial group", "Isolation method", "AMR detection method"),
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  
  labs(y = "Number of studies") +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

p_sankey2

ggsave(
  "captive_three_layer_sankey.jpeg",
  plot = p_sankey2,
  width = 28,
  height = 15,
  dpi = 300
)



## Patchwork the two together
library(patchwork)

combined_sankeys <- p_sankey / p_sankey2 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

combined_sankeys


ggsave(
  "wild_captive_sankey_combined.jpeg",
  plot = combined_sankeys,
  width = 26,
  height = 28,
  dpi = 300
)

ggsave(
  "wild_captive_sankey_combined.pdf",
  plot = combined_sankeys,
  width = 26,
  height = 28,
  dpi = 300
)



## -----------------------------------------------------------------------------
## Figure 2.8a: AMR outcome by species, wild


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

dataspecies_wild <- data.frame(
  Species = c("Red fox", "European badger", "Wolf spp.", "Eurasian otter", "Beech/Stone marten",
              "American mink", "Raccoon dog", "European polecat/Ferret", "European pine marten",
              "Golden jackal", "Least weasel", "European mink", "Stoat", "Steppe polecat", "'Fox'"),
  Total_studies = c(39, 21, 17, 14, 10, 3, 8, 7, 7, 2, 2, 1, 1, 1, 4),
  Resistant   = c(30, 11, 12, 10, 5, 2, 6, 2, 3, 2, 0, 1, 0, 0, 2),
  Susceptible = c(5, 7, 2, 2, 3, 1, 1, 3, 3, 0, 2, 0, 1, 1, 0),
  Unclear     = c(4, 3, 3, 2, 2, 0, 1, 2, 1, 0, 0, 0, 0, 0, 2)
)

data_wild_long <- dataspecies_wild %>%
  pivot_longer(
    cols = c("Resistant", "Susceptible", "Unclear"),
    names_to = "Result",
    values_to = "Count"
  ) %>%
  mutate(Result = factor(Result, levels = c("Resistant", "Susceptible", "Unclear")))

p10 <- ggplot(data_wild_long, aes(x = reorder(Species, Total_studies), y = Count, fill = Result)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Resistant"   = "#D6604D",
    "Susceptible" = "#4DAF4A",
    "Unclear"     = "#B3B3B3"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Host species", y = "Number of studies", fill = "Result") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

## -----------------------------------------------------------------------------
## Plot 2.8b: AMR outcome by species, captive


dataspecies_captive <- data.frame(
  Species = c("American mink", "'Mink'", "European polecat/Ferret", "Red fox",
              "European badger", "Eurasian otter", "Beech/Stone marten",
              "Raccoon dog", "'Fox'", "Wolf spp.", "Least weasel", "Arctic fox",
              "European mink"),
  Total_studies = c(23, 15, 9, 6, 5, 5, 5, 5, 3, 3, 3, 1, 1),
  Resistant   = c(15, 15, 4, 5, 3, 3, 4, 3, 3, 2, 0, 0, 1),
  Susceptible = c(2, 0, 0, 1, 2, 2, 1, 0, 0, 0, 3, 1, 0),
  Treatment   = c(6, 0, 4, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0)
)

data_captive_long <- dataspecies_captive %>%
  pivot_longer(
    cols = c("Resistant", "Susceptible", "Treatment"),
    names_to = "Result",
    values_to = "Count"
  ) %>%
  mutate(Result = factor(Result, levels = c("Resistant", "Susceptible", "Treatment")))

p11 <- ggplot(data_captive_long, aes(x = reorder(Species, Total_studies), y = Count, fill = Result)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Resistant"   = "#D6604D",
    "Susceptible" = "#4DAF4A",
    "Treatment"   = "#B3B3B3"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Host species", y = "Number of studies", fill = "Result") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )


## Combine


combined_1011 <- p10 / p11 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

combined_1011

ggsave(
  filename = "UPDATED_amr_by_species_combined.png",
  plot = combined_1011,
  width = 13,
  height = 12,
  dpi = 300
)

## -----------------------------------------------------------------------------

## Figure 2.9a — Wild: count antibiotic class occurrences per species


df_wild_raw <- read.csv("wild_host_resistance.csv", check.names = FALSE)

df_wild_raw <- df_wild_raw[, !(is.na(names(df_wild_raw)) | names(df_wild_raw) == "")]

wild_class_counts <- df_wild_raw %>%
  filter(!is.na(Species), Species != "",
         !is.na(Antibiotic_Class), Antibiotic_Class != "") %>%
  separate_rows(Antibiotic_Class, sep = ",") %>%
  mutate(
    Species = trimws(Species),
    Antibiotic_Class = trimws(Antibiotic_Class)
  ) %>%
  group_by(Species, Antibiotic_Class) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(Species, desc(Count))

write.csv(
  wild_class_counts,
  "antimicrobial_class_counts_wild.csv",
  row.names = FALSE
)

## Wild heatmap (built from the percentages file)


df_wild_pct <- read.csv("wild_host_resistance_percentages.csv", check.names = FALSE)

df_wild_pct <- df_wild_pct[, !(is.na(names(df_wild_pct)) | names(df_wild_pct) == "")]

df_wild_long <- df_wild_pct %>%
  filter(!is.na(Species_wild), Species_wild != "") %>%
  pivot_longer(
    cols = -Species_wild,
    names_to = "Antimicrobial_Class",
    values_to = "Percentage"
  ) %>%
  mutate(
    Species_wild = trimws(Species_wild),
    Antimicrobial_Class = trimws(Antimicrobial_Class),
    Antimicrobial_Class = factor(Antimicrobial_Class, levels = sort(unique(Antimicrobial_Class)))
  )

p12 <- ggplot(df_wild_long, aes(x = Antimicrobial_Class, y = Species_wild, fill = Percentage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    name = "% studies",
    limits = c(0, 100),
    guide = guide_colorbar(
      barwidth = 1,
      barheight = 8,
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  labs(x = "Antimicrobial class", y = "Species") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13, face = "bold", color = "black"),
    axis.text.y = element_text(size = 13, face = "bold", color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

## -----------------------------------------------------------------------------

## Figure 2.9b — Captive: count antibiotic class occurrences per species


df_captive_raw <- read.csv("captive_host_resistance.csv", check.names = FALSE)

df_captive_raw <- df_captive_raw[, !(is.na(names(df_captive_raw)) | names(df_captive_raw) == "")]

captive_class_counts <- df_captive_raw %>%
  filter(!is.na(Species), Species != "",
         !is.na(Antibiotic_Class), Antibiotic_Class != "") %>%
  separate_rows(Antibiotic_Class, sep = ",") %>%
  mutate(
    Species = trimws(Species),
    Antibiotic_Class = trimws(Antibiotic_Class)
  ) %>%
  group_by(Species, Antibiotic_Class) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(Species, desc(Count))

write.csv(
  captive_class_counts,
  "antimicrobial_class_counts_captive.csv",
  row.names = FALSE
)

## -----------------------------------------------------------------------------
## Figure 2.9b — Captive heatmap (built from the percentages file)


df_captive_pct <- read.csv("captive_host_resistance_percentages.csv", check.names = FALSE)

df_captive_pct <- df_captive_pct[, !(is.na(names(df_captive_pct)) | names(df_captive_pct) == "")]

df_captive_long <- df_captive_pct %>%
  filter(!is.na(Species_captive), Species_captive != "") %>%
  pivot_longer(
    cols = -Species_captive,
    names_to = "Antimicrobial_Class",
    values_to = "Percentage"
  ) %>%
  mutate(
    Species_captive = trimws(Species_captive),
    Antimicrobial_Class = trimws(Antimicrobial_Class),
    Antimicrobial_Class = factor(Antimicrobial_Class, levels = sort(unique(Antimicrobial_Class)))
  )

p13 <- ggplot(df_captive_long, aes(x = Antimicrobial_Class, y = Species_captive, fill = Percentage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    name = "% studies",
    limits = c(0, 100),
    guide = guide_colorbar(
      barwidth = 1,
      barheight = 8,
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) +
  labs(x = "Antimicrobial class", y = "Species") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13, face = "bold", color = "black"),
    axis.text.y = element_text(size = 13, face = "bold", color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

## Combine with ONE shared legend (both on the same 0-100% scale)


combined_1213 <- p12 / p13 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

combined_1213

ggsave(
  filename = "amr_profile_heatmaps_combined.png",
  plot = combined_1213,
  width = 12,
  height = 16,
  dpi = 300
)


################################################################################
################################################################################
################################################################################
################################################################################

## Stats 1. Chi-squared for region and species

chisq.test(table(regions$Region, regions$Species),
           simulate.p.value = TRUE)


unique(regions$Region)
unique(regions$Species)


## Stats 2. Association test between sampling method and sample type (wild)
chisq_wild <- chisq.test(table(df_clean$Sample_method, df_clean$Sample_type))
chisq_wild$expected
chisq.test(
  table(df_clean$Sample_method, df_clean$Sample_type),
  simulate.p.value = TRUE
)


## Stats 3: Do sample type OR AMR detection method influence resistance detection?

## Wild

# Read file
df <- read.csv("detection_resistance.csv", check.names = FALSE)


# Clean + split multiple entries
df_long <- df %>%
  separate_rows(Sample_type, sep = ",") %>%
  separate_rows(Test_method, sep = ";") %>%
  mutate(
    Sample_type = trimws(Sample_type),
    Test_method = trimws(Test_method),
    Resistance_detected = trimws(Resistance_detected)
  )

chisq.test(
  table(df_long$Sample_type, df_long$Resistance_detected),
  simulate.p.value = TRUE
)

chisq.test(
  table(df_long$Test_method, df_long$Resistance_detected),
  simulate.p.value = TRUE
)


## Captive

# Read file
df <- read.csv("captive_detection_resistance.csv", check.names = FALSE)


# Clean + split multiple entries
df_long <- df %>%
  separate_rows(Sample_type, sep = ",") %>%
  separate_rows(Test_method, sep = ";") %>%
  mutate(
    Sample_type = trimws(Sample_type),
    Test_method = trimws(Test_method),
    Resistance_detected = trimws(Resistance_detected)
  )



chisq.test(
  table(df_long$Sample_type, df_long$Resistance_detected),
  simulate.p.value = TRUE
)

chisq.test(
  table(df_long$Test_method, df_long$Resistance_detected),
  simulate.p.value = TRUE
)


chisq <- chisq.test(table(df_clean$Sample_method, df_clean$Sample_type))
chisq$expected

chisq.test(table(df_clean$Sample_method, df_clean$Sample_type),
           simulate.p.value = TRUE)


chisq <- chisq.test(table(df_captive_clean$Study_type, df_captive_clean$Sample_type))
chisq$expected

chisq.test(table(df_captive_clean$Study_type, df_captive_clean$Sample_type),
           simulate.p.value = TRUE)




################################################################################


# Convert all columns except Species to numeric

wild <- read.csv("wild_host_resistance_percentages.csv", check.names = FALSE)
captive <- read.csv("captive_host_resistance_percentages.csv", check.names = FALSE)

names(wild)
names(captive)

## calculate mean resistance for each class

wild_means <- c(
  mean(wild$`b-lactams`, na.rm = TRUE),
  mean(wild$Tetracyclines, na.rm = TRUE),
  mean(wild$Aminoglycosides, na.rm = TRUE),
  mean(wild$`Macrolides/Lincosamides`, na.rm = TRUE),
  mean(wild$`Sulfonamides/Trimethoprim`, na.rm = TRUE),
  mean(wild$`Quinolones/Fluoroquinolones`, na.rm = TRUE),
  mean(wild$Phenicols, na.rm = TRUE),
  mean(wild$Other, na.rm = TRUE)
)

captive_means <- c(
  mean(captive$`b-lactams`, na.rm = TRUE),
  mean(captive$Tetracyclines, na.rm = TRUE),
  mean(captive$Aminoglycosides, na.rm = TRUE),
  mean(captive$`Macrolides/Lincosamides`, na.rm = TRUE),
  mean(captive$`Sulfonamides/Trimethoprim`, na.rm = TRUE),
  mean(captive$`Quinolones/Fluoroquinolones`, na.rm = TRUE),
  mean(captive$Phenicols, na.rm = TRUE),
  mean(captive$Other, na.rm = TRUE)
)

## Stats 4: wilcox.test(wild_means, captive_means, paired = TRUE)
## Which antimicrobial classes differ most?

difference <- captive_means - wild_means
difference

median(wild_means)
median(captive_means)



# Read in files
wild <- read.csv("wild_resistance_percents.csv")
captive <- read.csv("captive_resistance_percents.csv")


# Make sure Percent_resistant is numeric
wild$Percent_resistant <- as.numeric(wild$Percent_resistant)
captive$Percent_resistant <- as.numeric(captive$Percent_resistant)

# Merge wild + captive by Species
merged <- merge(
  wild,
  captive,
  by = "Species"
)

## Stats 5: Species-level paired comparison
#
# Question:
# Is resistance generally higher in captive animals across species?


wilcox.test(
  merged$Percent_resistant.x,
  merged$Percent_resistant.y,
  paired = TRUE,
  exact = FALSE
)

## Stats 6: Spearman's correlation

# Question:
# Do species with high resistance in wild animals also show high resistance in captive animals?


cor.test(
  merged$Percent_resistant.x,
  merged$Percent_resistant.y,
  method = "spearman",
  exact = FALSE
)

## ---------------------------------------------------------------------------

## Supplementary materials

## Figure S2. put into supplementary: Sample types used, wild vs. captive (% of studies)

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

samples <- data.frame(
  Sample_type = c("Faecal sample", "Caecal contents", "Swab - Faecal",
                  "Swab - Perianal/Rectal/Perineal/Genital",
                  "Swab - Nasopharyngeal/Oral/Auditory canal",
                  "Swab - Skin/Wound/Internal organ",
                  "Tissue/Organ/Fluid sample", "Environmental sample",
                  "No details/Isolates from elsewhere", "Treatment paper"),
  Wild = c(53.5, 7, 4.2, 18.3, 8.5, 0, 19.7, 0, 4.2, 0),
  Captive = c(23.2, 0, 1.8, 14.3, 10.7, 5.4, 35.7, 12.5, 10.7, 17.9)
)

samples_long <- samples %>%
  pivot_longer(
    cols = c(Wild, Captive),
    names_to = "Source",
    values_to = "Count"
  )

p3 <- ggplot(samples_long, aes(x = reorder(Sample_type, Count), y = Count, fill = Source)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65,
           color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("Wild" = "#2166AC", "Captive" = "#E8834D")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Sample type", y = "% of studies", fill = "Study type") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.position = "top",
    legend.key.size = unit(0.5, "cm")
  )

p3





















