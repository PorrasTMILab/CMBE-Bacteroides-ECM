# Load necessary libraries
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tidyr)
library(readxl)
library(writexl)

# ------------------------------
# Step 1: Load combined annotation data and harmonize strain names
# ------------------------------

# Define base path to your data folder (adjust as needed)
base_path <- "../data"

# Load combined CAZyme + EggNOG annotation file
annot_df <- read_csv(file.path(base_path, "CAZyme_EggNOG_Combined.csv"))

# Harmonize strain names
harmonize_strain <- function(df, col) {
  mutate(df,
         !!sym(col) := case_when(
           !!sym(col) == "GCA_000613385.1"                ~ "B. acidifaciens",
           !!sym(col) == "GCA_000169015.1"                ~ "B. caccae",
           !!sym(col) == "GCA_000154525.1"                ~ "B. stercoris",
           !!sym(col) == "GCA_000381365.1"                ~ "B. salyersiae",
           !!sym(col) == "GCA_000210075.1"                ~ "B. xylanisolvens",
           !!sym(col) == "GCA_000172175.1"                ~ "B. intestinalis",
           !!sym(col) == "GCA_900107315.1"                ~ "B. uniformis",
           !!sym(col) == "GCA_028539425.1"                ~ "B. thetaiotaomicron",
           !!sym(col) == "GCA_900107475.1"                ~ "B. ovatus",
           !!sym(col) == "GCA_001997325.1"                ~ "B. fragilis",
           !!sym(col) == "GCA_045689915.1"                ~ "B. fragilis 43858",
           !!sym(col) == "GCA_000598905.1"                ~ "B. fragilis 43859",
           !!sym(col) == "GCA_001699875.1"                ~ "B. fragilis 43860",
           !!sym(col) == "GCA_025567135.1"                ~ "B. cellulolyticus",
           TRUE ~ !!sym(col)
         ))
}



# ------------------------------
# Step 2: Process EC annotations - combine dbCAN and eggNOG EC numbers
# ------------------------------

# Flatten or clean list-like EC columns to character strings
annot_df_clean <- annot_df %>%
  mutate(
    `dbCAN_EC#` = as.character(`dbCAN_EC#`),
    eggNOG_EC = as.character(eggNOG_EC),
    `dbCAN_EC#` = map_chr(`dbCAN_EC#`, ~ toString(.[[1]])),
    eggNOG_EC = map_chr(eggNOG_EC, ~ toString(.[[1]]))
  )

# Function to combine unique, valid EC numbers from both sources
combine_unique_ecs <- function(db, eg) {
  db <- ifelse(is.na(db) | db == "", NA_character_, db)
  eg <- ifelse(is.na(eg) | eg == "", NA_character_, eg)
  db_clean <- db %>%
    str_remove_all(":\\d+") %>%
    str_replace_all("\\|", ",") %>%
    str_remove_all("-")
  eg_clean <- str_remove_all(eg, "-")
  db_vec <- str_split(db_clean, ",")[[1]]
  eg_vec <- str_split(eg_clean, ",")[[1]]
  valid_ec <- "^\\d+\\.\\d+\\.\\d+\\.\\d+$"
  ecs <- unique(c(db_vec, eg_vec))
  ecs <- ecs[!is.na(ecs) & ecs != ""]
  ecs <- ecs[str_detect(ecs, valid_ec)]
  if (length(ecs) == 0) return(NA_character_)
  paste(sort(ecs), collapse = ",")
}

# Apply EC combination on each row
cazy_proc <- annot_df_clean %>%
  rowwise() %>%
  mutate(
    EC_combined = combine_unique_ecs(`dbCAN_EC#`, eggNOG_EC),
    dbCAN_Unique = map_chr(`dbCAN_Recommend Results`, ~ {
      if (is.na(.x) || .x == "") return(NA_character_)
      unique_vals <- str_split(.x, "\\|")[[1]] %>%
        str_remove_all("_e?\\d+") %>% 
        unique() %/% sort()
      paste(unique_vals, collapse = ",")
    })
  ) %>%
  ungroup() %>%
  filter(is.na(EC_combined) | !str_detect(EC_combined, "-")) %>%
  mutate(source = "CAZyme")

# ------------------------------
# Step 3: Load MEROPS data, filter and merge
# ------------------------------

merops_tbl <- read_csv(file.path(base_path, "all_csv_peptidases.csv")) %>%
  rename(E_value = `E-value`,
         Strain = Genome,
         Gene_ID = Query_ID) %>%
  filter(E_value <= 0.01,
         Peptidase_Class != "I",
         !str_detect(Peptidase_Family, "\\.[^.]*UN", ignore_case = TRUE)) %>%
  mutate(Unassigned_peptidase = if_else(
    str_detect(Peptidase_Family, "\\.[^.]*UP", ignore_case = TRUE),
    "unassigned", "known"))

merops_merged <- cazy_proc %>%
  full_join(merops_tbl, by = "Gene_ID") %>%
  mutate(
    merops_identified = if_else(!is.na(MEROPS_ID), "peptidase", NA_character_),
    source = "Protease" # overwrite source flag for merged MEROPS rows
  )

final_df <- merops_merged %>%
  separate_rows(EC_combined, sep = ",") %>%
  filter(EC_combined != "-") %>%
  select(c(
    "Source_File", "Gene_ID", "eggNOG_CAZy", "dbCAN_EC#", "dbCAN_dbCAN_sub",
    "dbCAN_Recommend Results", "Is_CAZyme", "Strain.x", "EC_combined", "dbCAN_Unique",
    "source", "MEROPS_ID", "E_value", "Description", "Peptidase_Family",
    "Peptidase_Class", "Strain.y", "Unassigned_peptidase", "merops_identified"
  ))

# Save cleaned mapping for downstream use
write_csv(final_df, file.path(base_path, "clean_EC_meropsID.csv"))

# ------------------------------
# Step 4: Load UniProt + BRENDA ECM database and clean
# ------------------------------

short_db_ECM <- read_csv(file.path(base_path, "UniProt_PLUS_BRENDA_CollapsedByEC.csv")) %>%
  select(EC_number, brenda_recommended_name, uniprot_cazy, uniprot_merops,
         uniprot_Protease_Type, `uniprot_protein names`, `uniprot_ec number`,
         brenda_Substrate_Category, brenda_substrate_product, brenda_natural_substrate_product) %>%
  mutate(brenda_recommended_name = str_remove(brenda_recommended_name, "^RN\\s+")) %>%
  separate_rows(brenda_Substrate_Category, sep = ",\\s*") %>%
  filter(!is.na(brenda_Substrate_Category) & brenda_Substrate_Category != "")

# Remove unwanted substrates
drop_subs <- c("adiponectin", "calpain", "mitogen", "cadherin")
short_db_ECM <- short_db_ECM %>%
  filter(!brenda_Substrate_Category %in% drop_subs)

# ------------------------------
# Step 5: Split and clean EC columns in ECM database
# ------------------------------

db_EC <- short_db_ECM %>%
  mutate(db_EC = coalesce(`uniprot_ec number`, EC_number)) %>%
  mutate(db_EC = na_if(db_EC, "-")) %>%
  separate_rows(db_EC, sep = ";") %>%
  mutate(db_EC = str_trim(db_EC))

# ------------------------------
# Step 6: Match final annotation dataframe to ECM database by EC
# ------------------------------

matched_by_ec <- final_df %>%
  inner_join(db_EC, by = c("EC_combined" = "db_EC")) %>%
  distinct()

# ------------------------------
# Step 7: Match by MEROPS family where applicable
# ------------------------------

short_long_mer <- short_db_ECM %>%
  mutate(uniprot_merops = str_trim(uniprot_merops)) %>%
  separate_rows(uniprot_merops, sep = ";") %>%
  mutate(uniprot_merops = na_if(uniprot_merops, ""))

final_df_proteases <- final_df %>%
  filter(!is.na(MEROPS_ID)) %>%
  mutate(Peptidase_Family = str_trim(Peptidase_Family))

match_by_merops <- final_df_proteases %>%
  inner_join(short_long_mer, by = c("Peptidase_Family" = "uniprot_merops")) %>%
  distinct()

# ------------------------------
# Step 8: Combine EC and MEROPS matches into one dataframe
# ------------------------------

ecm_enriched <- matched_by_ec %>%
  full_join(match_by_merops, by = "Gene_ID") %>%
  distinct()

# Resolve duplicate columns by preference for non-NA values
dup_cols <- gsub("\\.x$", "", grep("\\.x$", names(ecm_enriched), value = TRUE))
for (col in dup_cols) {
  ecm_enriched[[col]] <- coalesce(ecm_enriched[[paste0(col, ".x")]], ecm_enriched[[paste0(col, ".y")]])
}
ecm_enriched_clean <- ecm_enriched %>%
  select(-matches("\\.x$"), -matches("\\.y$"))

# Load gene-to-strain mapping for join
cazy_map <- annot_df %>%
  select(Gene_ID, Strain) %>%
  distinct()

ecm_enriched1 <- left_join(ecm_enriched_clean, cazy_map, by = "Gene_ID")

# ------------------------------
# Step 9: Collapse data per gene for final output
# ------------------------------

ecm_genes <- ecm_enriched1 %>%
  group_by(Gene_ID) %>%
  summarize(
    Strain = paste(unique(Strain), collapse = ", "),
    EC_combined = paste(unique(EC_combined), collapse = ", "),
    dbCAN_Unique = paste(unique(dbCAN_Unique), collapse = ", "),
    uniprot_cazy = paste(unique(uniprot_cazy), collapse = ", "),
    uniprot_merops = paste(unique(uniprot_merops), collapse = ", "),
    mr_family_matched = paste(unique(Peptidase_Family), collapse = ", "),
    brenda_Substrate_Category = paste(unique(brenda_Substrate_Category), collapse = ", "),
    brenda_recommended_name = paste(unique(brenda_recommended_name), collapse = ", "),
    .groups = "drop"
  ) %>%
  select(-Strain.y) %>%  # Remove redundant column
  rename(Strain = Strain.x)

# ------------------------------
# Step 10: Save combined mapping files
# ------------------------------

write_csv(ecm_enriched1, file.path(base_path, "ecm_mapped_combined.csv"))
write_csv(ecm_genes, file.path(base_path, "ecm_genes_collapsed.csv"))

