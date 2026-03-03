# Load required packages
library(dplyr)
library(purrr)
library(stringr)
library(httr)
library(readr)
library(tidyr)
library(jsonlite)

# ------------------------------
# Step 1: Set project directories (relative paths)
# ------------------------------
data_dir <- "../data"           # Folder for raw & intermediate files
output_dir <- "../results"      # Folder to save output CSVs
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------
# Step 2: Download and process BRENDA database (if needed)
# ------------------------------
# Note: If you already have 'full_BRENDA.csv', skip download step

# library(brendaDb)
# brenda_filepath <- DownloadBrenda()  # This opens license agreement
# brenda_df <- ReadBrenda(brenda_filepath)
# brenda_wide <- brenda_df %>%
#   pivot_wider(id_cols = ID,
#               names_from = "field",
#               values_from = "description",
#               values_fn = list(description = ~ paste(., collapse = "; "))) %>%
#   janitor::clean_names()
# write_csv(brenda_wide, file.path(data_dir, "full_BRENDA.csv"))

# ------------------------------
# Step 3: Load existing full BRENDA dataset CSV
# ------------------------------
full_brenda_df <- read_csv(file.path(data_dir, "full_BRENDA.csv"), show_col_types = FALSE)

# ------------------------------
# Step 4: Define extracellular matrix (ECM) core terms and regex pattern
# ------------------------------
ecm_terms <- c(
  "collagen", "elastin", "fibrillin", "matrix protein", "chondroadherin",
  "discoidin", "otogelin", "otolin", "periostin", "tectorin",
  "tenascin", "thrombospondin", "aggrecan", "biglycan", "brevican",
  "decorin", "keratocan", "perlecan", "serglycin", "versican",
  "agrin", "ameloblastin", "amelogenin", "asporin", "dermatopontin",
  "extracellular phosphoglycoprotein", "fibronectin", "fibrinogen",
  "fibromodulin", "fibulin", "gliomedin", "hemicentin", "hevin",
  "kielin", "laminin", "matrilin", "multimerin", "nephrocan",
  "nephronectin", "neurocan", "nidogen", "podocan", "opticin",
  "osteoglycin", "osteomodulin", "osteonectin", "sialoprotein",
  "vitrin", "vitronectin", "chondroitin", "hyaluron", "hyaluronic acid",
  "dermatan", "heparan", "keratan", "focal adhesion", "integrin",
  "netrin", "paxillin", "talin", "vinculin", "zyxin",
  "dentin", "extracellular matrix", "matrix phosphoglycoprotein"
)

ecm_regex <- regex(paste0("\\b(", paste(ecm_terms, collapse = "|"), ")\\b"), ignore_case = TRUE)

# ------------------------------
# Step 5: Filter BRENDA for ECM substrates/products & extract matched terms
# ------------------------------
substrate_hits <- full_brenda_df %>%
  filter(str_detect(natural_substrate_product, ecm_regex) |
         str_detect(substrate_product, ecm_regex)) %>%
  mutate(
    matched_natural = str_extract_all(natural_substrate_product, ecm_regex),
    matched_substrate = str_extract_all(substrate_product, ecm_regex)
  )

# Save filtered substrate hits to CSV for record
write_csv(substrate_hits, file.path(output_dir, "substrate_hits_raw_BRENDA.csv"))

# ------------------------------
# Step 6: UniProt API fields (URL-encoded)
# ------------------------------
fields_enc <- paste0(
  c("accession","reviewed","id","protein_name","gene_names","organism_name",
    "length","xref_merops","xref_cazy","xref_pfam",
    "ec","cc_function","cc_activity_regulation","cc_catalytic_activity",
    "ph_dependence","protein_existence","cc_interaction","cc_disease",
    "xref_embl","xref_glygen","xref_glyconnect","xref_glycosmos",
    "xref_carbonyldb","xref_unicarbkb","xref_eggnog","xref_panther",
    "kinetics","go_p"),
  collapse = "%2C"
)

# ------------------------------
# Step 7: Function to fetch UniProt data for a given EC number
# ------------------------------
fetch_uniprot_ec <- function(ec_number) {
  query <- sprintf("(ec:%s)", ec_number)
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?fields=",
    fields_enc,
    "&format=tsv&query=",
    URLencode(query, reserved = TRUE)
  )
  
  res <- httr::GET(url)
  if (httr::status_code(res) != 200) return(tibble())
  
  txt <- httr::content(res, "text", encoding = "UTF-8")
  if (nchar(txt) == 0) return(tibble())
  
  df <- read_tsv(txt, show_col_types = FALSE,
                 col_types = cols(.default = col_character())) %>%
    rename_with(tolower, everything()) %>%
    mutate(original_ec = ec_number)
  return(df)
}

# ------------------------------
# Step 8: Get unique EC numbers from substrate hits and query UniProt
# ------------------------------
ec_vec <- substrate_hits %>%
  pull(id) %>%
  unique()

all_results_tb <- purrr::map_dfr(ec_vec, fetch_uniprot_ec, .id = "call_nr")

# Collapse multivalued columns (optional, but recommended for clarity)
all_results_tb_collapsed <- all_results_tb %>%
  group_by(original_ec) %>%
  summarise(across(everything(), ~paste(sort(unique(.x)), collapse = "; ")), .groups = "drop")

# Save UniProt API results collapsed by EC number
write_csv(all_results_tb_collapsed, file.path(output_dir, "UniProt_API_results_full_CollapsedByEC.csv"))

# ------------------------------
# Step 9: Helper functions to classify protease and CAZyme types
# ------------------------------
protease_detect <- function(merops_string) {
  tokens <- str_split(merops_string, "\\s*;\\s*")[[1]]
  tokens <- tokens[tokens != ""]
  first_letters <- toupper(substr(tokens, 1, 1))
  letter_to_type <- c(
    S = "Serine", C = "Cysteine", M = "Metallopeptidase", T = "Threonine",
    A = "Aspartic", P = "Mixed/Unknown", G = "Glutamic",
    N = "Asparagine Lyases", U = "Unknown", I = "Inhibitor"
  )
  types <- letter_to_type[first_letters]
  types <- na.omit(types)
  if(length(types) == 0) NA_character_ else paste(unique(types), collapse = ", ")
}

cazy_detect <- function(cazy_id) {
  pats <- list(
    "Glycoside Hydrolase" = "GH\\d+",
    "Glycosyl Transferase" = "GT\\d+",
    "Polysaccharide Lyase" = "PL\\d+",
    "Carbohydrate Esterase" = "CE\\d+",
    "Carbohydrate-Binding Module" = "CBM\\d+",
    "Auxiliary Activity" = "AA\\d+"
  )
  
  tokens <- str_split(cazy_id, "\\s*;\\s*")[[1]] %>%
    stringr::str_trim() %>%
    discard(~ .x == "")
  
  matches <- names(pats)[purrr::map_lgl(pats, ~ any(str_detect(tokens, .x)))]
  if (length(matches) == 0) NA_character_ else paste(matches, collapse = ", ")
}

# ------------------------------
# Step 10: Annotate UniProt results with protease and CAZyme classes
# ------------------------------
all_results_tb_collapsed_flat <- all_results_tb_collapsed %>%
  mutate(
    Protease_Type = purrr::map_chr(xref_merops, protease_detect),
    CAZyme_Class = purrr::map_chr(xref_cazy, cazy_detect)
  )

write_csv(all_results_tb_collapsed_flat, file.path(output_dir, "UniProt_API_results_full_CollapsedByEC_annotated.csv"))

# ------------------------------
# Step 11: Prepare BRENDA data & merge with UniProt annotations
# ------------------------------
substrate_brenda <- substrate_hits %>%
  select(
    id, recommended_name, systematic_name, synonyms, matched_natural,
    matched_substrate, inhibitors, natural_substrate_product,
    substrate_product, reaction, reaction_type, ph_optimum, ph_range,
    temperature_optimum, temperature_range, application
  ) %>%
  rename(EC_number = id)

# Add prefixes to distinguish source columns
substrate_brenda_prefixed <- substrate_brenda %>%
  rename_with(~ paste0("brenda_", .x), -EC_number)
uniprot_prefixed <- all_results_tb_collapsed_flat %>%
  rename(EC_number = original_ec) %>%
  rename_with(~ paste0("uniprot_", .x), -EC_number)

# Merge datasets on EC number, keeping all BRENDA entries
merged_tbl <- left_join(substrate_brenda_prefixed, uniprot_prefixed, by = "EC_number")

write_csv(merged_tbl, file.path(output_dir, "UniProt_PLUS_BRENDA_CollapsedByEC.csv"))

# ------------------------------
# Step 12: Extract key columns for easier downstream use and save
# ------------------------------
merged_summary <- merged_tbl %>%
  select(
    EC_number,
    brenda_recommended_name,
    brenda_matched_natural,
    brenda_matched_substrate,
    uniprot_Protease_Type,
    uniprot_xref_merops,
    uniprot_CAZyme_Class,
    uniprot_xref_cazy,
    uniprot_xref_pfam,
    uniprot_cc_activity_regulation
  )

write_csv(merged_summary, file.path(output_dir, "shortened_UniProt_PLUS_BRENDA_CollapsedByEC.csv"))
