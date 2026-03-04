library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(openxlsx)

# Set project base directory relative to your script location
# Adjust this as needed to your folder structure
base_dir <- ".."

# ---- Export GAG CAZymes ----
gag_substrates <- c("dermatan", "heparan", "hyaluronic acid", "keratan", "chondroitin")
cazy_gags <- ecm_enriched1 %>%
  filter(brenda_Substrate_Category %in% gag_substrates) %>%
  select(Gene_ID, Strain, brenda_Substrate_Category, brenda_recommended_name,
         dbCAN_Unique, Is_CAZyme, EC_number, Peptidase_Class, Peptidase_Family)

write.xlsx(cazy_gags,
           file.path(base_dir, "KarenMA/Papers/CMBE_Bacteroides/Figures/Figure4_substrate_breakdown/cazyme_gags.xlsx"))

# ---- Export ECM proteins and proteases ----
protein_substrates <- c("collagen", "elastin", "fibrinogen", "osteonectin")
proteins_proteases <- ecm_enriched1 %>%
  filter(brenda_Substrate_Category %in% protein_substrates) %>%
  select(Gene_ID, Strain, brenda_Substrate_Category, brenda_recommended_name,
         EC_number, Peptidase_Class, Peptidase_Family, dbCAN_Unique)

write.xlsx(proteins_proteases,
           file.path(base_dir, "KarenMA/Papers/CMBE_Bacteroides/Figures/Figure4_substrate_breakdown/proteins_proteases.xlsx"))

# ---- Export proteoglycans ----
proteoglycan_substrates <- c("agrin", "biglycan", "decorin", "fibromodulin",
                             "perlecan", "serglycin", "versican")
proteoglycans_all_enzymes <- ecm_enriched1 %>%
  filter(brenda_Substrate_Category %in% proteoglycan_substrates) %>%
  select(Gene_ID, Strain, brenda_Substrate_Category, brenda_recommended_name,
         EC_number, Peptidase_Class, Peptidase_Family, dbCAN_Unique, Is_CAZyme)

write.xlsx(proteoglycans_all_enzymes,
           file.path(base_dir, "KarenMA/Papers/CMBE_Bacteroides/Figures/Figure4_substrate_breakdown/proteoglycans_all_enzymes.xlsx"))

# ---- Export glycoproteins ----
glycoprotein_substrates <- c("fibrillin", "fibronectin", "laminin", "nidogen",
                             "tenascin", "vitronectin", "thrombospo")
glycoprot_all_enzymes <- ecm_enriched1 %>%
  filter(brenda_Substrate_Category %in% glycoprotein_substrates) %>%
  select(Gene_ID, Strain, brenda_Substrate_Category, brenda_recommended_name,
         EC_number, Peptidase_Class, Peptidase_Family, dbCAN_Unique, Is_CAZyme)

write.xlsx(glycoprot_all_enzymes,
           file.path(base_dir, "KarenMA/Papers/CMBE_Bacteroides/Figures/Figure4_substrate_breakdown/glycoprot_all_enzymes.xlsx"))

# ---- Process CAZyme families from dbCAN and UniProt combined columns ----
split_families <- function(db, up) {
  if (!is.na(db) && db != "") {
    return(str_split(db, ",")[[1]])
  } else if (!is.na(up) && up != "") {
    return(str_split(up, ";")[[1]])
  }
  character(0)
}

ecm_long_cazfamily_2_sources <- ecm_enriched1 %>%
  mutate(
    dbCAN_Unique = as.character(dbCAN_Unique),
    uniprot_cazy = as.character(uniprot_cazy)
  ) %>%
  rowwise() %>%
  mutate(families = list(split_families(dbCAN_Unique, uniprot_cazy))) %>%
  ungroup() %>%
  unnest(families, keep_empty = TRUE) %>%
  rename(family = families) %>%
  mutate(
    family = str_replace_all(family, " ", ""),
    broad_cat = str_replace_all(gsub("[0-9]", "", sub("^([A-Za-z]+).*", "\\1", family)), " ", "")
  ) %>%
  distinct()

# ---- Filter to remove CBM, GT, and Peptidase Classes X, U ----
new_caz_counts1 <- ecm_long_cazfamily_2_sources %>%
  filter(!(brenda_Substrate_Category %in% c("collagen", "elastin", "fibrinogen",
                                            "osteonectin", "laminin", "fibronectin"))) %>%
  filter(
    !is.na(broad_cat),
    broad_cat != "",
    broad_cat != "CBM",
    broad_cat != "GT"
  )

# Remove X and U protease classes and filter relevant substrates
minus_x_u <- ecm_long_cazfamily_2_sources %>%
  filter(!(brenda_Substrate_Category %in% c("chondroitin", "dermatan", "heparan",
                                            "keratan", "hyaluronic acid"))) %>%
  filter(!Peptidase_Class %in% c("X", "U"),
         !is.na(Peptidase_Class),
         Peptidase_Class != "")

# Join filtered CAZymes and proteases
joint <- bind_rows(new_caz_counts1, minus_x_u) %>%
  distinct()


# ---- Prepare heatmap count matrix including certain substrates ----
intestine_subs <- c(
  "agrin", "collagen", "elastin", "laminin", "fibrillin", "periostin",
  "tenascin", "thrombospondin", "biglycan", "decorin", "perlecan",
  "serglycin", "versican", "fibronectin", "fibrinogen", "fibromodulin",
  "nidogen", "osteonectin", "vitronectin", "hyaluronic acid",
  "heparan", "keratan", "chondroitin", "dermatan"
)

count_mat_raw_minus <- joint %>%
  filter(brenda_Substrate_Category %in% intestine_subs) %>%
  select(Strain, brenda_Substrate_Category, total_enzyme) %>%
  pivot_wider(names_from = brenda_Substrate_Category, values_from = total_enzyme,
              values_fill = list(total_enzyme = 0)) %>%
  column_to_rownames("Strain")

long_counts_minus <- count_mat_raw_minus %>%
  rownames_to_column("Strain") %>%
  pivot_longer(-Strain, names_to = "substrate", values_to = "Gene_Count") %>%
  mutate(Gene_Count = ifelse(Gene_Count == 0, NA_real_, Gene_Count))

# Categories for substrates
category_map <- c(
  agrin = "Proteoglycan", biglycan = "Proteoglycan", decorin = "Proteoglycan",
  versican = "Proteoglycan", perlecan = "Proteoglycan", fibromodulin = "Proteoglycan",
  serglycin = "Proteoglycan", collagen = "Protein", elastin = "Protein",
  fibrinogen = "Protein", osteonectin = "Protein", fibronectin = "Glycoprotein",
  laminin = "Glycoprotein", fibrillin = "Glycoprotein", periostin = "Glycoprotein",
  vitronectin = "Glycoprotein", nidogen = "Glycoprotein", tenascin = "Glycoprotein",
  thrombospondin = "Glycoprotein", chondroitin = "Glycosaminoglycan",
  `hyaluronic acid` = "Glycosaminoglycan", heparan = "Glycosaminoglycan",
  keratan = "Glycosaminoglycan", dermatan = "Glycosaminoglycan"
)

col_annot_minus <- tibble(
  substrate = names(count_mat_raw_minus),
  category = category_map[names(count_mat_raw_minus)]
)

long_anno_minus <- long_counts_minus %>%
  left_join(col_annot_minus, by = c("substrate")) %>%
  mutate(Strain = factor(Strain, levels = c(
    "B. ovatus", "B. xylanisolvens", "B. salyersiae", "B. thetaiotaomicron",
    "B. caccae", "B. intestinalis", "B. uniformis", "B. acidifaciens",
    "B. fragilis", "B. stercoris", "B. cellulolyticus"
  )))

# ---- Plot heatmap ----
library(ggplot2)
library(RColorBrewer)

ggplot(long_anno_minus, aes(x = substrate, y = Strain, fill = Gene_Count)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = brewer.pal(7, "YlGnBu"),
    limits = c(0, max(long_anno_minus$Gene_Count, na.rm = TRUE))
  ) +
  geom_text(aes(label = Gene_Count), color = "black", size = 3) +
  facet_wrap(~ category, scales = "free_x", ncol = 4, strip.position = "bottom") +
  scale_x_discrete(position = "top") +
  theme_classic() +
  theme(
    axis.text.y = element_text(face = "italic", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 0, size = 10),
    axis.title = element_blank(),
    strip.text = element_text(size = 10),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white", color = "white"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  )

# ---- Save heatmap plot ----
ggsave(file.path(base_dir, "KarenMA/Papers/CMBE_Bacteroides/Preliminary Figs/heatmap_all_substrates_minus_nondegrading_UsXs.pdf"),
       width = 10, height = 4.7, dpi = 300)
