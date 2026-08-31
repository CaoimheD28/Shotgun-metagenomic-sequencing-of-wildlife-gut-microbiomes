## 1. Count how many unique species are present in the MAG dataset overall

library(tidyverse)

## 1. Read data 
df <- read_csv("MAGs_species.csv")

df <- df %>%
  rename(Species = `Species`)

## 2. Count unique species
n_unique_species <- n_distinct(df$Species)

cat("Total MAGs:", nrow(df), "\n")
cat("Unique species:", n_unique_species, "\n")

## 3. See the full list of unique species, with how many MAGs assigned to each
species_summary <- df %>%
  count(Species, name = "n_MAGs") %>%
  arrange(desc(n_MAGs))

species_summary


# Save

write_csv(species_summary, "unique_species_summary.csv")


## -----------------------------------------------------------------------------

## Figure 6.2. Species identity bar chart: top N species by prevalence

MAGs_species <- read_csv("MAGs_species.csv")

top_n_species <- 40


# How many samples each species was found in, overall
species_counts <- MAGs_species %>%
  count(Species, name = "n_samples") %>%
  arrange(desc(n_samples))

top_species <- species_counts %>%
  slice_max(n_samples, n = top_n_species, with_ties = FALSE) %>%
  pull(Species)

n_other_species <- length(setdiff(species_counts$Species, top_species))

# Collapse everything outside the top N into "Other"
df_plot <- MAGs_species %>%
  mutate(
    Species_grp = if_else(
      Species %in% top_species,
      as.character(Species),
      paste0("Other (", n_other_species, " species)")
    )
  )

plot_counts <- df_plot %>%
  count(Species_grp, Characteristic, name = "n")

# Order bars by total prevalence, but force "Other" to the bottom
other_label <- paste0("Other (", n_other_species, " species)")
species_order <- plot_counts %>%
  group_by(Species_grp) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(Species_grp == other_label, desc(total)) %>%
  pull(Species_grp)

plot_counts <- plot_counts %>%
  mutate(Species_grp = factor(Species_grp, levels = rev(species_order)))

pal <- c("#4F9D69", "#E58429", "#7d273a")
n_hosts <- n_distinct(plot_counts$Characteristic)

p5 <- ggplot(plot_counts, aes(x = Species_grp, y = n, fill = Characteristic)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = pal[seq_len(n_hosts)]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y   = element_text(size = 12, face = "italic"),
    axis.text.x   = element_text(size = 12),
    axis.title.x  = element_text(size = 11, margin = margin(t = 8)),
    axis.line     = element_line(linewidth = 0.4, color = "grey30"),
    axis.ticks    = element_line(linewidth = 0.3, color = "grey30"),
    legend.position = "top",
    legend.title  = element_text(size = 11, face = "bold"),
    legend.text   = element_text(size = 10),
    plot.title    = element_text(size = 12, face = "bold", margin = margin(b = 10)),
    plot.margin   = margin(10, 15, 10, 10)
  ) +
  labs(
    x = NULL, y = "Number of samples",
    title = paste0("Top ", top_n_species, " MAG species by sample prevalence"),
    fill = "Host type"
  )

p5
ggsave("species_prevalence_barchart.pdf", p5, width = 10, height = 10)


## -----------------------------------------------------------------------------

##  Figure 6.3. UpSet plot
##  Shows which species are shared vs unique between host types


install.packages("UpSetR")
library(UpSetR)
library(ggplot2)

upset_data <- MAGs_species %>%
  distinct(Characteristic, Species) %>%
  mutate(Present = 1) %>%
  pivot_wider(names_from = Characteristic, values_from = Present,
              values_fill = 0) %>%
  column_to_rownames("Species")

p3 <- upset(
  upset_data,
  nsets = ncol(upset_data),
  order.by = "freq",
  main.bar.color = "steelblue",
  sets.bar.color = "grey40",
  text.scale = 2
)


print(p3)

pdf("upset_plot.pdf", width = 8, height = 6)


## Find MAG species unique to each cohort (Rural_Badger, Rural_Fox,
## Urban_Fox) - i.e. the "single dot" columns of the UpSet plot



df <- read_csv("Overall_details.csv") %>%
  rename(Species = `Species identified by GTDB`)

## 1. Which cohort(s) is each species found in?
species_by_cohort <- df %>%
  distinct(Characteristic, Species)

species_n_cohorts <- species_by_cohort %>%
  group_by(Species) %>%
  summarise(n_cohorts = n_distinct(Characteristic), .groups = "drop")

## 2. Keep only species found in exactly ONE cohort
unique_species <- species_by_cohort %>%
  semi_join(
    species_n_cohorts %>% filter(n_cohorts == 1),
    by = "Species"
  )

## 3. Split into a named list, one vector per cohort
unique_by_cohort <- split(unique_species$Species, unique_species$Characteristic)

# View them individually, e.g.:
unique_by_cohort$Rural_Badger
unique_by_cohort$Rural_Fox
unique_by_cohort$Urban_Fox

# Quick counts per cohort
map_int(unique_by_cohort, length)

## 4. Save a clean table for supplementary material
unique_species %>%
  arrange(Characteristic, Species) %>%
  write_csv("unique_species_by_cohort.csv")

cat("\nSaved unique_species_by_cohort.csv -",
    nrow(unique_species), "species total, unique to a single cohort.\n")

## species SHARED across all three cohorts (the fully-connected


shared_all_three <- species_n_cohorts %>%
  filter(n_cohorts == n_distinct(df$Characteristic)) %>%
  pull(Species)

shared_all_three


## -----------------------------------------------------------------------------

## Figure 6.4. ARG presence/absence heatmap by species, split by Rural/Urban

## rural-specific/urban-specific/more-prevalent-in taxa.


library(tidyverse)


# ComplexHeatmap is on Bioconductor, not CRAN
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)   # for colour ramps


## STEP 1: Merge MAG species/metadata into the ARG file, once.

## 1. Read both files
args_df    <- read_csv("All_samples_ARGs_MAGs.csv")
species_df <- read_csv("MAGs_species.csv")

## 2. Check what the columns are actually called
cat("=== args_df columns ===\n")
print(names(args_df))
cat("\n=== species_df columns ===\n")
print(names(species_df))

## 3. Build the join key + environment column
## The ARG file's Bin_ID looks like "97005_bin10" (Sample_MAGID) -
## build the same key in species_df so they can be joined.
species_df <- species_df %>%
  mutate(
    Bin_ID = paste0(Sample, "_", MAGID),
    Environment = str_extract(Characteristic, "^[A-Za-z]+")  # "Rural_Badger" -> "Rural"
  )

## Sanity check: Bin_ID should be unique in species_df
if (any(duplicated(species_df$Bin_ID))) {
  warning("Duplicate Bin_IDs found in species_df - check Sample/MAGID combinations")
}

## 4. Join
## Note: args_df already has its own Species column, only pull
## the extra metadata (Characteristic/Environment/quality stats) from
## species_df - avoids a Species.x/Species.y clash from the join.
combined <- args_df %>%
  left_join(
    species_df %>%
      select(Bin_ID, Characteristic, Environment,
             Completeness, Contamination, MIMAG),
    by = "Bin_ID"
  )

## 5. Sanity check the join actually worked
n_unmatched <- sum(is.na(combined$Characteristic))
cat("\nRows in ARG file:", nrow(args_df), "\n")
cat("Rows successfully matched to MAG metadata:", nrow(combined) - n_unmatched, "\n")
cat("Unmatched rows:", n_unmatched, "\n")

if (n_unmatched > 0) {
  cat("\nExample unmatched Bin_IDs (check these against species_df$Bin_ID):\n")
  combined %>% filter(is.na(Characteristic)) %>% distinct(Bin_ID) %>% head(10) %>% print()
}

combined <- combined %>% filter(!is.na(Characteristic))
cat("Dropped", n_unmatched, "unmatched rows. Rows remaining:", nrow(combined), "\n")

## 6. Recode intrinsic efflux/regulatory genes
## These genes are typically intrinsic multidrug efflux pumps or
## regulators (not classic acquired resistance determinants), so their
## many-drug-class listing in the raw data isn't biologically
## meaningful for this figure. Recode to a single category.

intrinsic_efflux_genes <- c(
  "acrB", "acrD", "acrE", "acrF", "acrS", "bacA", "baeR", "baeS", "BcII",
  "cpxA", "CRP", "efmA", "efrA", "efrB", "Enterobacter_cloacae_acrA",
  "Escherichia_coli_acrA", "Escherichia_coli_mdfA", "evgA", "evgS", "gadW",
  "gadX", "H-NS", "Klebsiella_pneumoniae_acrA", "Klebsiella_pneumoniae_KpnG",
  "Klebsiella_pneumoniae_KpnH", "lmrB", "lmrD", "lmrP", "marA", "mdtA",
  "mdtB", "mdtC", "mdtE", "mdtF", "mdtG", "mdtH", "mdtM", "mdtN", "mdtO",
  "mdtP", "msbA", "msrC", "optrA", "oqxA", "oqxB", "ramA", "tolC", "vmlR"
)

recode_label <- "Efflux pump / intrinsic MDR"   # edit this label to taste

# Sanity check: flag any genes in the list that don't actually appear
# in the data - catches typos/capitalization mismatches before they
# silently do nothing
not_found <- setdiff(intrinsic_efflux_genes, unique(combined$ARG))
if (length(not_found) > 0) {
  cat("\nWARNING: these genes in intrinsic_efflux_genes were not found",
      "in combined$ARG (check spelling/capitalization):\n")
  print(not_found)
}

n_recoded <- sum(combined$ARG %in% intrinsic_efflux_genes)
cat("\nRecoding", n_recoded, "rows (",
    length(intersect(intrinsic_efflux_genes, unique(combined$ARG))),
    "distinct genes) to '", recode_label, "'\n")

combined <- combined %>%
  mutate(
    Resistance = if_else(ARG %in% intrinsic_efflux_genes, recode_label, Resistance)
  )

## 7. Save the clean merged file
write_csv(combined, "combined_ARGs_MAGs.csv")
cat("\nSaved combined_ARGs_MAGs.csv -", nrow(combined), "rows.\n")

## 8. Save the full species/MAG list separately
## Important: combined_ARGs_MAGs.csv only contains MAGs that had at
## least one ARG hit (since it's built from the ARG file). For other
## counts - need EVERY MAG, including ones with zero ARGs
species_clean <- species_df %>%
  select(Species, Characteristic, Environment, Bin_ID, Sample, MAGID,
         Completeness, Contamination, MIMAG)

write_csv(species_clean, "species_clean.csv")
cat("Saved species_clean.csv -", nrow(species_clean), "rows (full MAG list).\n")


combined <- read_csv("combined_ARGs_MAGs.csv")




## STEP 2. Build the species x ARG matrix (presence/absence)

hit_summary <- combined %>%
  distinct(Species, ARG) %>%
  mutate(present = 1)

mat <- hit_summary %>%
  pivot_wider(names_from = ARG, values_from = present, values_fill = 0) %>%
  column_to_rownames("Species") %>%
  as.matrix()


## STEP 3. Row annotations: MAG counts + rural/urban specificity

species_env_counts <- combined %>%
  distinct(Species, Bin_ID, Environment) %>%
  count(Species, Environment) %>%
  pivot_wider(names_from = Environment, values_from = n, values_fill = 0)

# Make sure both columns exist even if one environment is missing entirely
if (!"Rural" %in% names(species_env_counts)) species_env_counts$Rural <- 0
if (!"Urban" %in% names(species_env_counts)) species_env_counts$Urban <- 0

species_env_counts <- species_env_counts %>%
  mutate(
    n_MAGs = Rural + Urban,
    Pattern = case_when(
      Rural > 0 & Urban == 0 ~ "Rural-specific",
      Urban > 0 & Rural == 0 ~ "Urban-specific",
      Rural > Urban          ~ "More prevalent Rural",
      Urban > Rural          ~ "More prevalent Urban",
      TRUE                   ~ "Equal"
    )
  ) %>%
  column_to_rownames("Species")

# Align to the matrix row order
species_env_counts <- species_env_counts[rownames(mat), ]

Pattern_colors <- c(
  "Rural-specific"        = "#0F6E56",
  "More prevalent Rural"  = "#7FC7B4",
  "Urban-specific"        = "#B24C63",
  "More prevalent Urban"  = "#E4A6B4",
  "Equal"                 = "grey85"
)

row_anno <- rowAnnotation(
  `n MAGs` = anno_text(species_env_counts$n_MAGs, gp = gpar(fontsize = 8)),
  Pattern = species_env_counts$Pattern,
  col = list(Pattern = Pattern_colors),
  annotation_legend_param = list(
    Pattern = list(title = "Rural/Urban pattern")
  ),
  simple_anno_size = unit(3, "mm")
)


## STEP 4. Column annotation: colour ARGs by resistance drug class

## Each ARG can have several classes listed in `Resistance` (semicolon separated)
## take the FIRST listed class as the "primary" one for colouring the top strip

arg_class <- combined %>%
  distinct(ARG, Resistance) %>%
  mutate(primary_class = str_split(Resistance, ";") %>% map_chr(1)) %>%
  distinct(ARG, .keep_all = TRUE) %>%
  filter(ARG %in% colnames(mat))

# Align to matrix column order
arg_class <- arg_class[match(colnames(mat), arg_class$ARG), ]

class_colors <- c(
  "aminoglycoside" = "#E8A33D",
  "cephamycin"     = "#4C72B0",
  "lincosamide"    = "#D65F5F",
  "macrolide"      = "#8172B2",
  "tetracycline"   = "#55A868"
)

# The full (unfiltered) dataset will have more classes than the 5 above
# (e.g. glycopeptide, phenicol, rifamycin, fluoroquinolone...) - fill
# in colours for any classes not already named, so nothing ends up
# grey/NA in the plot
unique_classes <- unique(arg_class$primary_class)
missing_classes <- setdiff(unique_classes, names(class_colors))

if (length(missing_classes) > 0) {
  extra_colors <- scales::hue_pal()(length(missing_classes))
  names(extra_colors) <- missing_classes
  class_colors <- c(class_colors, extra_colors)
}

cat("Drug classes found:", length(unique_classes), "\n")
cat(paste(unique_classes, collapse = ", "), "\n")

col_anno <- HeatmapAnnotation(
  `Drug class` = arg_class$primary_class,
  col = list(`Drug class` = class_colors),
  annotation_name_side = "left",
  simple_anno_size = unit(3, "mm")
)


## STEP 5. Plot heatmap

ht <- Heatmap(
  mat,
  name = "ARG presence",
  col = c("0" = "grey95", "1" = "#0F4C81"),
  rect_gp = gpar(col = "white", lwd = 0.5),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,   # we'll build a custom legend below
  top_annotation = col_anno,
  left_annotation = row_anno,
  row_names_gp = gpar(fontface = "italic", fontsize = 12),
  column_names_gp = gpar(fontsize = 11),
  column_names_rot = 45,
  border = FALSE
)

# Custom legend explaining presence/absence
presence_legend <- Legend(
  labels = c("Absent", "Present"),
  legend_gp = gpar(fill = c("grey95", "#0F4C81")),
  title = "ARG presence"
)

pdf_width  <- max(11, ncol(mat) * 0.25)
pdf_height <- max(16, nrow(mat) * 0.2)

pdf("arg_heatmap.pdf", width = pdf_width, height = pdf_height)
draw(ht, annotation_legend_list = list(presence_legend),
     merge_legend = TRUE)
dev.off()


## -----------------------------------------------------------------------------


## Figure 6.5 Plasmid replicon vs ARG/VFG overlap plot


## STEP 1. SET PATH
input_file <- "C:/PATH/TO/plasmid_arg_vfg_table.csv"
output_dir <- "C:/PATH/TO/OUTPUT_FOLDER"


## STEP 2. CLASSIFY EACH ROW INTO ARG/VFG OVERLAP CATEGORY
has_arg <- trimws(tolower(df[[arg_col]])) == "yes"
has_vfg <- trimws(tolower(df[[vfg_col]])) == "yes"

overlap_category <- ifelse(has_arg & has_vfg, "Both ARG & VFG",
                           ifelse(has_arg & !has_vfg, "ARG only",
                                  ifelse(!has_arg & has_vfg, "VFG only",
                                         "Neither")))
df$overlap_category <- factor(overlap_category,
                              levels = c("Neither", "VFG only", "ARG only", "Both ARG & VFG"))

## STEP 3. TALLY COUNTS PER REPLICON TYPE x OVERLAP CATEGORY
tally <- table(df[[plasmid_col]], df$overlap_category)

# Order replicon types by total count (most common at top of the plot)
replicon_totals <- rowSums(tally)
tally <- tally[order(replicon_totals), , drop = FALSE]  # ascending

cat("--- MAG counts per replicon type and ARG/VFG overlap category ---\n\n")
print(tally)

## STEP 4. SAVE TABLE 
tally_out <- as.data.frame.matrix(tally)
tally_out <- cbind(Replicon = rownames(tally_out), tally_out)
write.csv(tally_out, file.path(output_dir, "replicon_arg_vfg_overlap.csv"), row.names = FALSE)

## STEP 5. PLOT
plot_path <- file.path(output_dir, "replicon_arg_vfg_overlap.png")
n_replicons <- nrow(tally)
fig_height <- max(5, 0.5 * n_replicons + 2)

png(plot_path, width = 10, height = fig_height, units = "in", res = 200)
layout(matrix(1:2, nrow = 2), heights = c(fig_height - 1.3, 1.3))

par(mar = c(4, 16, 3, 2))  # wide left margin for full replicon names
colors <- c("Neither" = "#d3d3d3", "VFG only" = "#f1c40f",
            "ARG only" = "#e74c3c", "Both ARG & VFG" = "#8e44ad")

barplot(
  t(tally),
  horiz = TRUE,
  las = 1,               # horizontal axis labels (readable replicon names)
  col = colors[colnames(tally)],
  border = NA,
  xlab = "Number of MAGs",
  main = "Plasmid replicon carriage vs. ARG/VFG presence",
  cex.names = 0.85
)

par(mar = c(0, 1, 0, 1))
plot.new()
legend(x = "center", ncol = 4, bty = "n", cex = 0.9,
       legend = names(colors), fill = colors)


# ------------------------------------------------------------------------------


## Figure 6.6 MAG Functional Profile Heatmap (R)

# STEP 1. SET PATHS 
input_dir  <- "C:/PATH/TO/FEATURES.TSV"
output_dir <- "C:/PATH/TO/OUTPUT_FOLDER"

library(pheatmap)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## STEP 2 FUNCTIONAL CATEGORY DEFINITIONS
# Keyword-based bins applied to each gene's DFAST "product" description.

CATEGORY_KEYWORDS <- list(
  # Checked first: MGE-identifying terms are highly specific (plasmid,
  # transposase, phage, etc.), so they should take priority over more generic
  # keywords in later categories that could otherwise intercept them first --
  # e.g. "plasmid segregation protein" would match Replication/Repair's
  # "segregation" keyword before ever reaching this category, if this were
  # checked later. Putting MGE first avoids that.
  "Mobile Genetic Elements" = c("transposase", "conjugal transfer", "integrase",
                                "insertion sequence", "IS3 family", "IS30 family",
                                "prophage", "phage", "plasmid"),
  "Translation/Ribosome" = c("ribosomal protein", "ribosome", "translation",
                             "aminoacyl-tRNA", "tRNA ligase", "tRNA synthetase",
                             "tRNA", "elongation factor", "initiation factor", "release factor",
                             "GTPase", "trigger factor", "SsrA-binding"),
  "Transcription" = c("RNA polymerase", "sigma factor", "transcriptional regulator",
                      "transcription regulator", "transcription factor",
                      "transcriptional repressor", "transcriptional activator",
                      "transcription termination", "transcription antitermination",
                      "helix-turn-helix", "catabolite control protein"),
  "Replication/Repair/Recombination" = c("DNA polymerase", "DNA gyrase", "DNA topoisomerase",
                                         "DNA repair", "helicase", "primase", "DNA ligase",
                                         "recombinase", "DNA methyltransferase",
                                         "restriction endonuclease", "nuclease", "resolvase",
                                         "Holliday junction", "DNA-binding protein",
                                         "DNA-processing", "segregation", "condensation protein",
                                         "SMC-Scp", "glycosylase", "mismatch repair",
                                         "recombination regulator"),
  "Transport" = c("ABC transporter", "transporter", "permease", "efflux", "exporter",
                  "symporter", "antiporter", "channel", "porin", "ATP-binding cassette",
                  "ATPase", "solute carrier", "solute-binding", "ATP-binding protein",
                  "transport protein", "iron transport", "preprotein translocase",
                  "signal recognition particle"),
  "Signal Transduction" = c("sensor histidine kinase", "response regulator",
                            "two-component system", "chemotaxis", "signal transduction",
                            "cyclic di-GMP", "GAF domain", "pheromone",
                            "nitrogen regulator", "P-II family"),
  "Cell Wall/Membrane/Envelope" = c("peptidoglycan", "penicillin-binding", "outer membrane",
                                    "lipopolysaccharide", "glycosyl transferase", "glycosyltransferase",
                                    "glycoside hydrolase", "cell wall", "membrane protein",
                                    "lipoprotein", "teichoic acid", "cell division",
                                    "transglycosylase", "septation", "cell cycle protein",
                                    "FtsK", "envelope stress"),
  "Motility" = c("flagell", "pilus", "pilin", "fimbria"),
  "Energy Metabolism" = c("NADH dehydrogenase", "cytochrome", "ATP synthase",
                          "electron transport", "oxidoreductase", "dehydrogenase", "oxidase",
                          "reductase", "iron-sulfur", "Fe-S", "thioredoxin", "ferredoxin",
                          "dioxygenase"),
  "Central/Carbon Metabolism" = c("pyruvate", "acetyl-CoA", "citrate", "glycolysis",
                                  "kinase", "synthase", "carboxylase", "transferase",
                                  "isomerase", "hydrolase", "epimerase", "racemase",
                                  "phosphatase", "phosphoesterase", "esterase", "deaminase",
                                  "ligase", "lyase", "aldolase", "deacetylase", "peptidase",
                                  "osidase", "dehydratase", "transketolase", "transacylase",
                                  "ethanolamine utilization", "ethanolamine utilisation", "microcompartment",
                                  "mutase", "transaminase", "phosphorylase", "protease",
                                  "acyl carrier protein"),
  "Defence/Stress Response" = c("CRISPR", "toxin", "antitoxin", "heat shock", "chaperon",
                                "superoxide dismutase", "catalase", "peroxidase",
                                "resistance protein", "multidrug resistance", "stress protein",
                                "cold shock", "cold-shock", "nucleotide exchange factor"),
  "Secondary Metabolism" = c("polyketide synthase", "non-ribosomal peptide synthetase",
                             "terpene synthase", "bacteriocin", "secondary metabolite")
)
UNCATEGORIZED_LABEL <- "Other/Unclassified"
HYPOTHETICAL_LABEL  <- "Hypothetical/Unknown"
DOMAIN_ONLY_LABEL   <- "Uncharacterised (domain-only, e.g. DUF)"
ALL_CATEGORIES <- c(names(CATEGORY_KEYWORDS), DOMAIN_ONLY_LABEL, UNCATEGORIZED_LABEL, HYPOTHETICAL_LABEL)

categorize_product <- function(product) {
  if (is.na(product) || trimws(product) == "") return(HYPOTHETICAL_LABEL)
  p <- tolower(product)
  # NOTE: "uncharacterized protein" (American spelling) is matched here
  # deliberately, alongside the British spelling, because it needs to match
  # DFAST's actual output text, which follows standard NCBI/GenBank
  # nomenclature (American spelling) -- this isn't a comment/label, so it's
  # left as real annotation tools would write it, with both forms covered.
  if (grepl("hypothetical", p) || trimws(p) %in% c("unknown", "uncharacterized protein", "uncharacterised protein")) {
    return(HYPOTHETICAL_LABEL)
  }
  for (cat_name in names(CATEGORY_KEYWORDS)) {
    kws <- tolower(CATEGORY_KEYWORDS[[cat_name]])
    if (any(vapply(kws, function(k) grepl(k, p, fixed = TRUE), logical(1)))) {
      return(cat_name)
    }
  }
  # Conserved-domain-but-unknown-function fallback (Pfam DUF families, generic
  # "XXX family protein" / "XXX domain-containing protein" / repeat-scaffold
  # proteins like TPR or pentapeptide-repeat with no other functional cue).
  # Unanchored since a gene symbol is often appended after the generic
  # description, e.g. "zinc ribbon domain-containing protein YjdM". These are
  # functionally uncharacterised, but distinct from a plain "hypothetical
  # protein" hit since a specific conserved domain WAS detected -- just not
  # one with a known role yet.
  if (grepl("\\bduf[0-9]+\\b", p) ||
      grepl("family protein", p) ||
      grepl("domain[- ]containing protein", p) ||
      grepl("domain protein", p) ||
      grepl("repeat[- ]?containing protein", p) ||
      grepl("repeat protein", p)) {
    return(DOMAIN_ONLY_LABEL)
  }
  UNCATEGORIZED_LABEL
}

# ---- STEP 3. FIND & PARSE FILES
# First, list EVERYTHING in the folder (including subfolders) so we can see
# exactly what's there and compare against what matches our expected naming
# pattern -- this reveals any files that are being silently skipped because
# their name doesn't end in "_features.tsv" (different extension, typo,
# sitting in a subfolder, etc.)
all_files_in_dir <- list.files(input_dir, full.names = TRUE, recursive = TRUE)
# Exclude anything inside output_dir (e.g. results from a previous run of this
# script) so they don't get flagged as unexpected/non-matching files below.
output_dir_norm <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
all_files_norm <- normalizePath(all_files_in_dir, winslash = "/", mustWork = FALSE)
all_files_in_dir <- all_files_in_dir[!startsWith(all_files_norm, output_dir_norm)]
cat(sprintf("Total files found anywhere under %s (excluding output folder): %d\n", input_dir, length(all_files_in_dir)))

files <- list.files(input_dir, pattern = "_features\\.tsv$", full.names = TRUE,
                    ignore.case = TRUE, recursive = TRUE)
if (length(files) == 0) {
  stop(paste("No '*_features.tsv' files found in", input_dir))
}
cat(sprintf("Files matching '*_features.tsv' pattern: %d\n", length(files)))

if (length(all_files_in_dir) != length(files)) {
  non_matching <- setdiff(all_files_in_dir, files)
  cat(sprintf("\n*** %d file(s) in the folder did NOT match the '*_features.tsv' pattern: ***\n",
              length(non_matching)))
  print(non_matching)
  cat("Check these filenames above -- likely a different extension, spelling,\n")
  cat("or they're feature.tsv files with different naming than expected.\n\n")
}

mag_name_from_path <- function(path) {
  sub("_features\\.tsv$", "", basename(path), ignore.case = TRUE)
}

parse_features_tsv <- function(path) {
  df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE,
                   colClasses = "character", check.names = FALSE)
  colnames(df) <- tolower(trimws(colnames(df)))
  if ("feature" %in% colnames(df)) {
    df <- df[toupper(df$feature) == "CDS", , drop = FALSE]
  }
  df
}

# ---- 4b. CHECK FOR DUPLICATE MAG NAMES BEFORE PROCESSING ----------------------

mag_names_all <- vapply(files, mag_name_from_path, character(1))
dup_names <- unique(mag_names_all[duplicated(mag_names_all)])
if (length(dup_names) > 0) {
  cat("\n*** WARNING: duplicate MAG names detected -- these files would overwrite each other! ***\n")
  for (d in dup_names) {
    cat(sprintf("  MAG name '%s' is shared by:\n", d))
    for (f in files[mag_names_all == d]) cat(sprintf("      %s\n", f))
  }
  cat("Rename these files so each produces a unique MAG identifier, then re-run.\n\n")
}

# ---- 5. BUILD MAG x CATEGORY MATRIX ------------------------------------------
rows <- list()
failed_files <- character(0)
for (f in files) {
  mag <- mag_name_from_path(f)
  result <- tryCatch({
    df <- parse_features_tsv(f)
    if (!"product" %in% colnames(df)) {
      stop("no 'product' column found (check file has a header row with a 'product' column)")
    }
    cats <- vapply(df$product, categorize_product, character(1))
    counts <- table(factor(cats, levels = ALL_CATEGORIES))
    rows[[mag]] <- as.integer(counts)
    cat(sprintf("  %s: %d CDS processed\n", mag, nrow(df)))
    TRUE
  }, error = function(e) {
    cat(sprintf("  FAILED: %s -- %s\n", f, conditionMessage(e)))
    FALSE
  })
  if (!isTRUE(result)) failed_files <- c(failed_files, f)
}

mat <- do.call(rbind, rows)
colnames(mat) <- ALL_CATEGORIES
mat <- as.data.frame(mat, check.names = FALSE)

# ---- 5b. EXPLICIT ACCOUNTING: did every input file make it into the matrix? ---
n_unique_mags_expected <- length(unique(mag_names_all)) -
  (if (length(failed_files) > 0) length(unique(mag_name_from_path(failed_files))) else 0)
cat(sprintf("\n%d files found (%d unique MAG names) -> %d MAGs in final matrix.\n",
            length(files), length(unique(mag_names_all)), nrow(mat)))
if (nrow(mat) < n_unique_mags_expected) {
  cat("*** Some MAGs did not make it into the matrix and it wasn't from duplicate names or failed files above ***\n")
} else if (length(dup_names) > 0) {
  cat(sprintf("(The gap between files found and MAGs in the matrix is explained by the %d duplicate name(s) flagged above.)\n", length(dup_names)))
}
if (length(failed_files) > 0) {
  cat("\nThe following files FAILED to process and were skipped:\n")
  print(failed_files)
}

# Save the underlying matrices (raw counts and % of annotated CDS per MAG)
write.csv(mat, file.path(output_dir, "genome_x_category_counts.csv"))
rel_mat <- sweep(as.matrix(mat), 1, rowSums(mat), FUN = "/") * 100
write.csv(rel_mat, file.path(output_dir, "genome_x_category_percent.csv"))
cat(sprintf("\nSaved matrix (%d MAGs x %d categories) to %s\n", nrow(mat), ncol(mat), output_dir))

# ---- STEP 6. VISUALISATION
# Panel A: annotation completeness per MAG (% Hypothetical + % Uncharacterised
#          domain-only vs. % actually assigned to a named functional category).
#          These "categories of ignorance" tend to have the largest values and
#          visually swamp a shared colour scale if plotted alongside real
#          functional categories, so they get their own plot instead.
# Panel B: heatmap of ONLY the named functional categories, scaled per-column
#          (z-score) so that real variation in smaller categories (e.g.
#          Motility, Defence/Stress) is visible rather than being compressed
#          near zero by categories with much larger raw percentages.

named_categories <- names(CATEGORY_KEYWORDS)  # excludes Hypothetical/Uncharacterised/Other

# --- Panel A: annotation completeness - Figure S6.1
completeness_df <- data.frame(
  MAG = rownames(mat),
  Hypothetical = rel_mat[, HYPOTHETICAL_LABEL],
  Uncharacterised = rel_mat[, DOMAIN_ONLY_LABEL],
  Other = rel_mat[, UNCATEGORIZED_LABEL]
)
completeness_df$Annotated <- 100 - completeness_df$Hypothetical -
  completeness_df$Uncharacterised - completeness_df$Other
completeness_df <- completeness_df[order(-completeness_df$Annotated), ]

completeness_path <- file.path(output_dir, "annotation_completeness.png")
png(completeness_path, width = 9, height = max(6, 0.08 * nrow(mat)) + 1.4, units = "in", res = 200)
layout(matrix(1:2, nrow = 2), heights = c(6, 1.6))

par(mar = c(5, 4, 4, 2))
bar_mat <- t(as.matrix(completeness_df[, c("Annotated", "Hypothetical", "Uncharacterised", "Other")]))
barplot(
  bar_mat,
  horiz = TRUE,
  col = c("#31688e", "#c0392b", "#f1c40f", "#95a5a6"),
  border = NA,
  names.arg = rep("", nrow(completeness_df)),  # MAG labels hidden, consistent with the heatmap
  xlab = "% of total CDS",
  main = sprintf("Annotation completeness across %d MAGs", nrow(mat))
)

par(mar = c(0, 1, 0, 1))
plot.new()
legend(x = "center", ncol = 2, bty = "n", cex = 0.95, y.intersp = 1.6,
       legend = c("Annotated (named category)", "Hypothetical/Unknown",
                  "Uncharacterised (domain-only)", "Other/Unclassified"),
       fill = c("#31688e", "#c0392b", "#f1c40f", "#95a5a6"))
dev.off()
cat(sprintf("Saved annotation completeness plot: %s\n", completeness_path))

# --- Panel B: functional-category-only heatmap, per-column (z-score) scaled,
#              with a genus annotation strip alongside the rows ---
functional_only_mat <- rel_mat[, named_categories, drop = FALSE]

# z-scoring requires non-zero variance in each column; a category that's 0%
# (or identical) across every MAG can't be scaled and would break clustering.
col_sd <- apply(functional_only_mat, 2, sd)
zero_var_cols <- names(col_sd)[col_sd == 0 | is.na(col_sd)]
if (length(zero_var_cols) > 0) {
  cat(sprintf("\nNote: dropping %d categories with no variation across MAGs (all identical/zero), so they can't be z-scored: %s\n",
              length(zero_var_cols), paste(zero_var_cols, collapse = ", ")))
  functional_only_mat <- functional_only_mat[, !(colnames(functional_only_mat) %in% zero_var_cols), drop = FALSE]
}

# Genus extraction from MAG names
extract_genus <- function(mag_name) {
  tokens <- strsplit(mag_name, "_")[[1]]
  if (length(tokens) < 2) return("Unmatched")
  tokens[length(tokens) - 1]  # second-to-last "_"-chunk = genus
}

genus_vec <- vapply(rownames(functional_only_mat), extract_genus, character(1))
unmatched <- rownames(functional_only_mat)[genus_vec == "Unmatched"]
if (length(unmatched) > 0) {
  cat(sprintf("\nNote: %d MAG name(s) don't have a genus_species separator and are labelled 'Unmatched':\n",
              length(unmatched)))
  print(unmatched)
}

annotation_row <- data.frame(Genus = genus_vec, row.names = rownames(functional_only_mat))

# Build a distinct colour for each genus present

genus_palette <- c("#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
                   "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990")
genus_levels <- sort(unique(genus_vec))
genus_colors <- setNames(genus_palette[seq_along(genus_levels)], genus_levels)
annotation_colors <- list(Genus = genus_colors)

heatmap_path <- file.path(output_dir, "functional_profile_heatmap.png")
pheatmap(
  functional_only_mat,
  scale = "column",   # <- key fix: z-score each category independently so
  #    smaller categories aren't visually flattened by
  #    categories with larger raw percentages
  color = colorRampPalette(c("#440154", "#31688e", "#35b779", "#fde725"))(100),
  clustering_method = "average",
  show_rownames = FALSE,      # <- individual MAGs not labelled, as requested
  show_colnames = TRUE,
  fontsize_col = 10,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  main = sprintf("Functional profile across %d MAGs\n(relative variation per category, named categories only)", nrow(mat)),
  filename = heatmap_path,
  width = 10,
  height = 7
)


###############################################################################

## Genus clustering order check


# 1. SET PATH
results_dir <- "C:/PATH/TO/OUTPUT_FOLDER"


percent_path <- file.path(results_dir, "genome_x_category_percent.csv")
if (!file.exists(percent_path)) {
  stop(paste("Could not find", percent_path, "-- run analyze_mags_LAB.R first."))
}

percent_df <- read.csv(percent_path, row.names = 1, check.names = FALSE)

# 2. RESTRICT TO THE 12 NAMED FUNCTIONAL CATEGORIES 
# (excludes Hypothetical/Unknown, Uncharacterised, Other/Unclassified, exactly
# as the heatmap does)
excluded_cols <- c("Hypothetical/Unknown", "Uncharacterised (domain-only, e.g. DUF)", "Other/Unclassified")
functional_only_mat <- percent_df[, !(colnames(percent_df) %in% excluded_cols), drop = FALSE]

# Drop zero-variance columns exactly as the main script does (can't z-score these)
col_sd <- apply(functional_only_mat, 2, sd)
zero_var_cols <- names(col_sd)[col_sd == 0 | is.na(col_sd)]
if (length(zero_var_cols) > 0) {
  cat(sprintf("Dropping %d zero-variance categories (same as main script): %s\n\n",
              length(zero_var_cols), paste(zero_var_cols, collapse = ", ")))
  functional_only_mat <- functional_only_mat[, !(colnames(functional_only_mat) %in% zero_var_cols), drop = FALSE]
}

# 3. REPRODUCE PHEATMAP'S EXACT CLUSTERING 
# pheatmap scales (z-scores) columns FIRST when scale="column", then clusters
# rows on Euclidean distance of the SCALED matrix, using the clustering_method
# given (we used "average"). Replicating that exactly here.
scaled_mat <- scale(as.matrix(functional_only_mat))  # column-wise z-score, same as pheatmap scale="column"
row_dist <- dist(scaled_mat, method = "euclidean")
row_hclust <- hclust(row_dist, method = "average")

# This is the exact row order the heatmap dendrogram displays, top to bottom
ordered_mags <- rownames(functional_only_mat)[row_hclust$order]

# 4. EXTRACT GENUS FOR EACH MAG (same logic as main script)
extract_genus <- function(mag_name) {
  tokens <- strsplit(mag_name, "_")[[1]]
  if (length(tokens) < 2) return("Unmatched")
  tokens[length(tokens) - 1]
}
ordered_genus <- vapply(ordered_mags, extract_genus, character(1))

# 5. PRINT THE FULL ORDERED LIST
cat("--- MAGs in dendrogram order (top to bottom of the heatmap) ---\n\n")
for (i in seq_along(ordered_mags)) {
  cat(sprintf("%3d. [%s] %s\n", i, ordered_genus[i], ordered_mags[i]))
}

# 6. DETECT CONTIGUOUS SAME-GENUS BLOCKS
run_lengths <- rle(ordered_genus)
block_df <- data.frame(
  genus = run_lengths$values,
  block_size = run_lengths$lengths,
  stringsAsFactors = FALSE
)

# Total count of each genus overall, to compute "X of Y clustered together"
genus_totals <- table(ordered_genus)

cat("\n\n--- Contiguous same-genus blocks detected in the dendrogram order ---\n\n")
for (i in seq_len(nrow(block_df))) {
  g <- block_df$genus[i]
  block_size <- block_df$block_size[i]
  total_for_genus <- genus_totals[[g]]
  if (block_size > 1) {
    cat(sprintf("  %s: block of %d contiguous MAGs (out of %d total %s MAGs in the dataset)\n",
                g, block_size, total_for_genus, g))
  }
}

cat("\n(Genera with no line above appear only as single, isolated MAGs in the dendrogram order --\n")
cat(" i.e. never adjacent to another MAG of the same genus.)\n")

# ---- 7. SUMMARY PER GENUS: LARGEST BLOCK VS TOTAL COUNT ---------------------
cat("\n--- Summary: largest contiguous block per genus vs. total MAGs of that genus ---\n\n")
all_genera <- sort(unique(ordered_genus))
for (g in all_genera) {
  blocks_for_g <- block_df$block_size[block_df$genus == g]
  largest_block <- max(blocks_for_g)
  total <- genus_totals[[g]]
  cat(sprintf("  %-20s largest contiguous block = %d / %d total MAGs\n", g, largest_block, total))
}


################################################################################
# Functional consistency across all LAB MAGS
# do these MAGs occupy a similar functional niche overall?

# STEP 1. SET PATH
results_dir <- "C:/PATH/TO/OUTPUT_FOLDER"


percent_path <- file.path(results_dir, "genome_x_category_percent.csv")
if (!file.exists(percent_path)) {
  stop(paste("Could not find", percent_path, "-- run analyze_mags_LAB.R first."))
}

percent_df <- read.csv(percent_path, row.names = 1, check.names = FALSE)

# Focus on the 12 named functional categories (the "real" biology), excluding
# the three annotation-quality categories which aren't functions themselves.
excluded_cols <- c("Hypothetical/Unknown", "Uncharacterised (domain-only, e.g. DUF)", "Other/Unclassified")
functional_only <- percent_df[, !(colnames(percent_df) %in% excluded_cols), drop = FALSE]

# STEP 2. COMPUTE MEAN, SD, AND COEFFICIENT OF VARIATION PER CATEGORY
summary_df <- data.frame(
  Category = colnames(functional_only),
  Mean_pct = round(colMeans(functional_only), 2),
  SD = round(apply(functional_only, 2, sd), 2),
  Min = round(apply(functional_only, 2, min), 2),
  Max = round(apply(functional_only, 2, max), 2)
)
summary_df$CV_percent <- round(100 * summary_df$SD / summary_df$Mean_pct, 1)
summary_df <- summary_df[order(summary_df$CV_percent), ]
rownames(summary_df) <- NULL

cat(sprintf("--- Functional category consistency across all %d MAGs ---\n", nrow(functional_only)))
cat("(sorted from most consistent/similar to most variable)\n\n")
print(summary_df, row.names = FALSE)

# STEP 3. OVERALL SUMMARY
mean_cv <- round(mean(summary_df$CV_percent), 1)
n_low_cv <- sum(summary_df$CV_percent < 30)
n_high_cv <- sum(summary_df$CV_percent >= 50)


