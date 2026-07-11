# Load necessary libraries
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(tidyr)
library(forcats)
library(openxlsx)
library(tidytext)       # for reorder_within
library(scales)         # for ggplot2 scales

# Set base directory relative to the current script location
base_dir <- ".."  # Parent folder of 'scripts' directory

# Load data files
ecm_enriched1 <- read_csv(file.path(base_dir, "data", "ecm_mapped_combined.csv"))

# Define strain order: decreasing genome size
strain_order <- c(
  "B. ovatus", "B. thetaiotaomicron", "B. intestinalis", "B. xylanisolvens",
  "B. salyersiae", "B. fragilis", "B. acidifaciens", "B. uniformis",
  "B. caccae", "B. stercoris", "B. cellulolyticus"
)

# ---- Dotplot preparation for Protease families ----
merops_families <- ecm_enriched1 %>%
  filter(!is.na(Peptidase_Family)) %>%
  mutate(
    family_ID = str_extract(Peptidase_Family, "^[^.]+"),
    family_group = substr(family_ID, 1, 1),
    family_group = ifelse(family_group %in% c("U", "X"), "other", family_group),
    Strain = factor(Strain, levels = strain_order)
  ) %>%
  distinct(Strain, Gene_ID, family_ID, family_group) %>%
  group_by(Strain, family_ID, family_group) %>%
  summarise(n_genes = n_distinct(Gene_ID), .groups = "drop") %>%
  group_by(family_ID, family_group) %>%
  mutate(prevalence = n_distinct(Strain),
         total_genes = sum(n_genes)) %>%
  ungroup() %>%
  group_by(family_group) %>%
  arrange(desc(prevalence), desc(total_genes), .by_group = TRUE) %>%
  mutate(rank_in_group = row_number()) %>%
  ungroup() %>%
  mutate(
    level_label = paste0(family_group,"_", rank_in_group),
    family_order = factor(level_label, levels = unique(level_label))
  ) %>%
  select(Strain, family_ID, family_group, n_genes, prevalence, total_genes, family_order)

# Plot Protease dotplot
palette_prot <- c(
  "S" = "#543f4c", "M" = "#8a5c7d", "C" = "#cd9ab7", "A" = "#c9ae8f",
  "N" = "#f37d1a", "T" = "#b8525c", "other" = "#ADADAD"
)

ggplot(merops_families, aes(x = Strain, y = family_order, size = n_genes, color = family_group)) +
  geom_count(shape = 16, show.legend = TRUE) +
  scale_y_reordered() +
  scale_size_continuous(name = "Gene count", range = c(2, 6), breaks = c(1, 5, 9)) +
  scale_x_discrete(position = "top", expand = c(0.05, 0)) +
  scale_color_manual(values = palette_prot, name = "Protease class") +
  facet_grid(family_group ~ ., scales = "free_y", space = "free") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4.5))) +
  theme_minimal() +
  theme(
    strip.background.y = element_blank(),
    strip.text.y.left  = element_blank(),
    strip.text.y.right = element_text(angle = 270, hjust = 0.5, size = 12),
    panel.grid.major = element_line(color="gray70", linetype = "dotted", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.text.x.top = element_text(angle = 50, vjust = 0, hjust = 0, face = "oblique", size = 13),
    panel.spacing.y = unit(0.6, "lines"),
    axis.text.y = element_text(color = "black", size = 9),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

# Save protease dotplot
ggsave(file.path(base_dir, "/bubbleplot_protease_families.pdf"),
       width = 6, height = 12, dpi = 300)


# ---- Dotplot preparation for CAZymes ----

# Load matched annotation with CAZymes
ecm_matched <- read_csv(file.path(base_dir, "KarenMA/Papers/ECM_database/ECM_all_enzymes_matched_annotation.csv"))

cazymes <- ecm_enriched1 %>%
  filter(Is_CAZyme == TRUE) %>%
  select(Strain, Gene_ID, EC_combined, uniprot_cazy, brenda_recommended_name,
         brenda_Substrate_Category, `dbCAN_Recommend Results`)

cazy_dotplot_df <- cazymes %>%
  separate_rows(`dbCAN_Recommend Results`, sep = "\\|") %>%
  mutate(
    cazy_col = sub("_.*", "", `dbCAN_Recommend Results`),
    cazy_col = trimws(cazy_col),
    cazy_ID = cazy_col,
    cazy_group = str_extract(as.character(cazy_ID), "^[A-Za-z]{1,3}"),
    Strain = factor(Strain, levels = strain_order)
  ) %>%
  distinct(Strain, Gene_ID, cazy_ID, cazy_group) %>%
  group_by(Strain, cazy_ID, cazy_group) %>%
  summarise(n_genes = n_distinct(Gene_ID), .groups = "drop") %>%
  filter(!cazy_group %in% c("GT", "CBM"))

# Plot CAZyme dotplot
palette_cazy <- c(
  "GH" = "#2a5c64",
  "PL" = "#86acae",
  "CE" = "#808080"
)

cazy_dotplot_df <- cazy_dotplot_df %>%
  group_by(cazy_group) %>%
  mutate(
    prevalence = n_distinct(Strain),
    total_genes = sum(n_genes)
  ) %>%
  arrange(desc(prevalence), desc(total_genes)) %>%
  ungroup() %>%
  mutate(
    rank_in_group = row_number(),
    level_label = paste0(cazy_group, "_", rank_in_group),
    family_order = factor(level_label, levels = unique(level_label)),
    family_ordered = tidytext::reorder_within(cazy_ID, rank_in_group, cazy_group)
  ) %>%
  mutate(cazy_group = factor(cazy_group, levels = names(palette_cazy)))

ggplot(cazy_dotplot_df, aes(x = Strain, y = family_ordered, size = n_genes, color = cazy_group)) +
  geom_count(shape = 16) +
  scale_y_reordered() +
  scale_size_continuous(name = "Gene count", range = c(2, 6), breaks = c(1, 13, 27)) +
  scale_color_manual(values = palette_cazy, name = "CAZyme class") +
  facet_grid(cazy_group ~ ., scales = "free_y", space = "free") +
  scale_x_discrete(position = "top", expand = c(0.05, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    strip.background.y = element_blank(),
    strip.text.y.left = element_blank(),
    strip.text.y.right = element_text(angle = 270, hjust = 0.5, size = 12),
    panel.grid.major = element_line(color = "gray70", linetype = "dotted", linewidth = 0.3),
    panel.grid.minor.x = element_blank(),
    axis.text.x.top = element_text(angle = 50, hjust = 0, vjust = 0, face = "oblique", size = 13),
    panel.spacing.y = unit(0.6, "lines"),
    axis.text.y = element_text(color = "black", size = 9),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  )

# Save CAZyme dotplot
ggsave(file.path("/bubbleplot_cazy_filtered_color.pdf"),
       width = 6, height = 12, dpi = 300)

