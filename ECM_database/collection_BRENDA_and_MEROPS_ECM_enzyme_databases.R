

#packages loading
```{r}
library(XML)
library(readr)
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(plotly)
library(gapminder)
library(paletteer)
library(httr)
library(cluster)
library(pheatmap)
library(ggdendro)
library(cowplot)
library(stringr)
library(ggrepel)
library(ggrounded)
library(forcats)
library(grid) 
library(ggplot2)  
library(ggtext)
library(fillpattern)
library(ggpattern)
library(httr)
library(jsonlite)
library(DBI)
library(RSQLite)
library(tidyverse)
library(rvest)
library(RSQLite)
library(dplyr)
library(httr)
library(readr)
library(dplyr)
library(stringr)
library(ggalluvial)
library(viridis)
library(igraph)
library(ggraph)
library(scales)
library(circlize) 
 library(magick)
 library(schrute)
 library(treemapify)
library(dplyr)
library(stringr)
library(purrr)

```



#directory and load files
```{r}


```

#Download full BRENDA dataset
```{r}
#generating full BRENDA database
library(brendaDb)


# Download the file (this will open the license agreement link in your browser)
brenda.filepath <- DownloadBrenda()


df <- ReadBrenda(brenda.filepath)


df_substrates <- df %>%
  pivot_wider(id_cols = ID,
              names_from = "field",
              values_from = "description",
              values_fn = list(description = function(x) paste(x, collapse = "; "))
) %>%
  clean_names()


write_csv(df_substrates, file.path("folder/full_BRENDA.csv"))

```


##select columns within full BRENDA dataset
```{r}

#selecting only the columns that will be used to expand their content into their own separate columns or to collapse all of the information (usually short content like names)
df_with_subcols <- full_brenda_df %>%
  select( c(
    "recommended_name",
    "systematic_name",
    "synonyms",
    "reaction",
    "natural_substrate_product",
    "substrate_product",
    "reference"
    
  ))
```

#Parsing NATURAL_SUBSTRATE_PRODUCT
```{r Parsing NATURAL_SUBSTRATE_PRODUCT}

# -----------------------------
# 1. Start with the full dataset
# -----------------------------
parsed_nsp_table <- full_brenda_df %>%
  
  # -----------------------------------------
  # 2. Select only the columns you want to keep
  #    (id + the NSP column)
  # -----------------------------------------
  select(id,recommended_name, natural_substrate_product) %>%
  
  # ---------------------------------------------------------
  # 3. Split the NSP column into individual NSP blocks
  #    Each block begins with "NSP #"
  # ---------------------------------------------------------
  mutate(blocks = str_split(
    natural_substrate_product,
    regex("(?=NSP\\s*#)")
  )) %>%
  unnest(blocks) %>%
  filter(blocks != "") %>%
  mutate(blocks = str_trim(blocks)) %>%
  
  # -----------------------------------------
  # 4. Parse each NSP block using your function
  # -----------------------------------------
  mutate(parsed = map(blocks, parse_NSP_block)) %>%
  
  # ---------------------------------------------------------
  # 5. Unnest parsed tibble with a prefix for clarity
  #    This avoids name collisions and keeps provenance clear
  # ---------------------------------------------------------
  unnest(parsed, names_sep = "_nsp_") %>%
  # Rename parsed columns to your preferred names
  rename(
    nsp_ids        = parsed_nsp_ids,
    nsp_reaction   = parsed_nsp_reaction,
    nsp_commentary = parsed_nsp_commentary,
    nsp_flag       = parsed_nsp_flag,
    nsp_literature = parsed_nsp_literature
  ) %>%

  
  # -----------------------------------------
  # 6. Keep only id + parsed NSP columns
  # -----------------------------------------
  # 6. Keep only id + the new NSP columns
  select(id,recommended_name, nsp_ids, nsp_reaction, nsp_commentary, nsp_flag, nsp_literature)






```

#SUBSTRATE_PRODUCT
```{r SUBSTRATE_PRODUCT}

parse_SP_block <- function(block) {

  # 1. Extract IDs between the first two hashes
  ids <- str_match(block, "^SP\\s*#([^#]+)#")[,2] %>% str_trim()
  ids <- if_else(is.na(ids) | ids == "", NA_character_, ids)

  # 2. Remove the SP header + IDs
  body <- str_replace(block, "^SP\\s*#[^#]+#\\s*", "")

  # 3. Extract reaction: everything before first (, {, <, or end
  react <- str_match(body, "^(.*?)\\s*(\\(|\\{|<|$)")[,2] %>% str_trim()

  # fallback if parentheses missing but '#' appears
  if (is.na(react) || react == "") {
    react <- str_match(body, "^(.*?)\\s*(#|\\{|<)")[,2] %>% str_trim()
  }

  # 4. Remove reaction from body
  body_no_react <- str_remove(body, fixed(react))

  # 5. Commentary inside parentheses
  commentary <- str_match(body_no_react, "\\(([^)]*)\\)")[,2] %>% str_trim()
  commentary <- if_else(is.na(commentary) | commentary == "", NA_character_, commentary)

  # 6. Flag inside {}
  flag <- str_match(body, "\\{([^}]*)\\}")[,2] %>% str_trim()
  flag <- if_else(is.na(flag) | flag == "", NA_character_, flag)

  # 7. Literature IDs inside last <>
  literature <- str_match(body, "<([^>]+)>$")[,2] %>% str_trim()
  literature <- if_else(is.na(literature) | literature == "", NA_character_, literature)

  tibble(
    ids        = ids,
    reaction   = react,
    commentary = commentary,
    flag       = flag,
    literature = literature,
    raw_block  = block
  )
}


parsed_sp_table <- full_brenda_df %>%
  
  # 1. Select only columns of interest
  select(id,recommended_name, substrate_product) %>%
  
  # 2. Split into SP blocks
  mutate(blocks = str_split(
    substrate_product,
    regex("(?=SP\\s*#)")
  )) %>%
  unnest(blocks) %>%
  filter(blocks != "") %>%
  mutate(blocks = str_trim(blocks)) %>%
  
  # 3. Parse each block
  mutate(parsed = map(blocks, parse_SP_block)) %>%
  
  # 4. Unnest with prefix
  unnest(parsed, names_sep = "_sp_") %>%
  
  # 5. Rename parsed columns
  rename(
    sp_ids        = parsed_sp_ids,
    sp_reaction        = parsed_sp_reaction,
    sp_commentary = parsed_sp_commentary,
    sp_flag       = parsed_sp_flag,
    sp_literature = parsed_sp_literature
  ) %>%
  
  # 6. Keep only id + parsed columns
  select(id,recommended_name, sp_ids, sp_reaction, sp_commentary, sp_flag, sp_literature)



```

#Unified substrate table
```{r}
unified_substrate_table <- bind_rows(
  parsed_nsp_table %>%
    mutate(type = "natural_substrate") %>%
    rename(
      ids        = nsp_ids,
      reaction   = nsp_reaction,
      commentary = nsp_commentary,
      flag       = nsp_flag,
      literature = nsp_literature
    ),
  
  parsed_sp_table %>%
    mutate(type = "substrate_product") %>%
    rename(
      ids        = sp_ids,
      reaction   = sp_reaction,
      commentary = sp_commentary,
      flag       = sp_flag,
      literature = sp_literature
    )
) %>%
  select(id,recommended_name, type, ids, reaction, commentary, flag, literature)


write_csv(unified_substrate_table, file.path( ONE,"folder/parsed_substrates_BRENDA.csv"))



```

#ECM terms
```{r}
ecm_substrates <- c(
  "collagen", "elastin", "fibrillin", "matrix protein", "chondroadherin",
  "discoidin", "otogelin", "otolin", "periostin", "tectorin", "tenascin",
  "thrombospondin",
  "aggrecan", "biglycan", "brevican", "decorin", "keratocan", "perlecan",
  "serglycin", "versican",
  "agrin", "ameloblastin", "amelogenin", "asporin", "dermatopontin",
  "extracellular phosphoglycoprotein", "fibronectin", "fibrinogen",
  "fibromodulin", "fibulin", "gliomedin", "hemicentin", "hevin", "kielin",
  "laminin", "matrilin", "multimerin", "nephrocan", "nephronectin",
  "neurocan", "nidogen", "podocan", "opticin", "osteoglycin",
  "osteomodulin", "osteonectin", "sialoprotein", "vitrin", "vitronectin",
  "chondroitin", "hyaluron", "hyaluronic acid", "dermatan", "heparan",
  "keratan",
  "focal adhesion", "integrin", "netrin", "paxillin", "talin", "vinculin",
  "zyxin",
  "dentin", "extracellular matrix", "matrix phosphoglycoprotein"
)

filtered_substrates <- unified_substrate_table %>%
  filter(
    str_detect(
      str_to_lower(reaction),
      str_c(str_to_lower(ecm_substrates), collapse = "|")
    )
  )

unified_with_match <- filtered_substrates %>%
  mutate(
    matched_term = map_chr(
      reaction,
      ~ {
        rxn <- str_to_lower(.x)
        hit <- ecm_substrates[str_detect(rxn, str_to_lower(ecm_substrates))]
        if (length(hit) == 0) NA_character_ else hit[1]
      }
    )
  ) %>%
  select(c(
    "EC_Number"="id",
    "recommended_name",
    "matched_term",
    "reaction",     
    "commentary",
    "type"
    
  )
  )

write_csv(unified_with_match, "folder/ECMsubstrate_matches_reaction_BRENDA.csv")

ecm_summary <- unified_with_match %>%
  filter(!is.na(matched_term)) %>%
  group_by(matched_term) %>%
  summarize(
    n_unique_EC = n_distinct(id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_unique_EC))

write_csv(ecm_summary, "folder/ECMsubstrate_matches_summary_BRENDA.csv")


```



#Download and filter MEROPS dataset
```{r}
# -------------------------------------------------------------------
#  The list of substrings you want to find
# -------------------------------------------------------------------
ecm_substrates <- c(
  "collagen", "elastin", "fibrillin", "matrix protein", "chondroadherin",
  "discoidin", "otogelin", "otolin", "periostin", "tectorin", "tenascin",
  "thrombospondin",
  "aggrecan", "biglycan", "brevican", "decorin", "keratocan", "perlecan",
  "serglycin", "versican",
  "agrin", "ameloblastin", "amelogenin", "asporin", "dermatopontin",
  "extracellular phosphoglycoprotein", "fibronectin", "fibrinogen",
  "fibromodulin", "fibulin", "gliomedin", "hemicentin", "hevin", "kielin",
  "laminin", "matrilin", "multimerin", "nephrocan", "nephronectin",
  "neurocan", "nidogen", "podocan", "opticin", "osteoglycin",
  "osteomodulin", "osteonectin", "sialoprotein", "vitrin", "vitronectin",
  "chondroitin", "hyaluron", "hyaluronic acid", "dermatan", "heparan",
  "keratan",
  "focal adhesion", "integrin", "netrin", "paxillin", "talin", "vinculin",
  "zyxin",
  "dentin", "extracellular matrix", "matrix phosphoglycoprotein"
)

# -------------------------------------------------------------------
#  Grab the file from the FTP server -----------------------------
# -------------------------------------------------------------------
url <- "https://ftp.ebi.ac.uk/pub/databases/merops/current_release/database_files/Substrate_search.txt"

# The file is tab‑delimited, has no header, and every field is wrapped
# in single quotes.  read_tsv() does the heavy lifting for you.


raw <- read_delim(
  file = url,
  delim = "\t",
  col_names = FALSE,
  escape_double = FALSE,
  trim_ws = FALSE
)

# ------------------------------------------------------------------
#   Strip the surrounding single quotes from every column
# ------------------------------------------------------------------
# The whole file is wrapped in single‑quote characters; we drop them.
remove_quotes <- function(x) {
  str_remove_all(x, "^'|'$")
}
# ------------------------------------------------------------------
#   Keep only the columns of interest, and rename them
# ------------------------------------------------------------------
subset_df <- raw %>% 
  select(Cleavage_event_ID = X1, 
    MEROPS_identifier  = X2, 
    Substrate_name  = X3, 
    Substrate_identifier = X4, 
    UniProt_accession = X14, 
    Position_in_substrate = X15, 
    organism = X16, 
    Peptidase_name = X17, 
    Method = X21, 
    Inhibitor_status = X22, 
    Physiological_or_artificial = X23) %>%
  mutate(across(everything(), remove_quotes))  %>%
  
  #Replace the literal string "NULL" with TRUE NA everywhere
  mutate(across(everything(), ~ na_if(.x, "NULL"))) 
    
    



# ------------------------------------------------------------------
#   Search for your substrates in the *Substrate_name* column
# ------------------------------------------------------------------
# Build a single regex that matches any of the listed substrings,
# case‑insensitively and with word boundaries.

substr_pat <- paste0("(?i)\\b(", paste(ecm_substrates, collapse = "|"), ")\\b")

merops_db_clean1 <- subset_df %>%
  mutate(
    contains_term = str_detect(Substrate_name, substr_pat)
  ) %>%
  filter(contains_term == TRUE) %>%

  mutate(
    matched_term = map_chr(
      Substrate_name,
      ~ {
        txt <- str_to_lower(.x)
        hit <- ecm_substrates[str_detect(txt, str_to_lower(ecm_substrates))]
        if (length(hit) == 0) NA_character_ else hit[1]
      }
    )
  ) %>%
  select(c(
    MEROPS_identifier,
    matched_term,
    Substrate_name,
    Peptidase_name,
    
  )
    
  ) %>%
  distinct()



# ------------------------------------------------------------------
#   Collapse filtered MEROPS dataset by family, grouping all substrates in the family
# ------------------------------------------------------------------

collapsed_merops <- merops_db_clean1 %>%
  mutate( 
    family_merops   = ifelse(is.na(MEROPS_identifier) | MEROPS_identifier == "",
                             NA_character_,
                             sub("\\..*", "", MEROPS_identifier)),   # everything before the dot
    subfamily_merops = ifelse(is.na(MEROPS_identifier) | !grepl("\\.", MEROPS_identifier),
                              NA_character_,
                              sub(".*\\.", "", MEROPS_identifier))    # everything after the dot
  ) %>%
  group_by(family_merops) %>%
  summarise(across(
    .cols = everything(),
    .fns = ~ paste(unique(na.omit(.x)), collapse = "; "),
    .names = "{.col}"
  )) %>%
  ungroup()
   


# write_xlsx(collapsed_merops, file.path(ONE, "folder/MEROPS_families_substrates_db.xlsx"))
```


