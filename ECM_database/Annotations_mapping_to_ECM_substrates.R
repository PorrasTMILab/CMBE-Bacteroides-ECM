


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
library(openxlsx)
library(writexl)

```





```{r}
rescued <- read_xlsx(file.path( ONE,"KarenMA/Papers/ECM_database/removed_rows_filtered.xlsx"))


new_db_minus_rescue <- anti_join(
  new_db, rescued,
  by = "EC_number"
  
)



new_db <- read_csv(file.path( ONE,"KarenMA/Papers/ECM_database/parsed_ECM_db.csv"))

new_db_full <- read_csv(file.path( ONE,"KarenMA/Papers/ECM_database/parsed_UniProt_API_results_full_CollapsedByEC.csv"))

#cazyme_eggnog <- read_csv(file.path( ONE,"KarenMA/Papers/CMBE_Bacteroides/Database_files/CAZyme_EggNOG_Combined_renamed.csv"))

#eggnog in the full 
new_db_cog <- new_db_full %>%
  select(c(
    "original_ec",
    "protein names"   , 
    `ec number`,
    merops,
    Protease_Type,
    cazy,
    CAZyme_class,
    "interacts with",
    
    "eggnog"
  ))

```


```{r}
#Load database files

db_ECM_BRENDA_revised

MEROPS_families_collapsed <- db_ECM_MEROPS_protease_families

```



#Load and bind annotation files for each species
```{r}
#=================================================================================
#PROKKA FILES
#=================================================================================

# List all CSV files in the folder
prokka_csv_files <- list.files(
  path = file.path("folder/prokka_results"),
  pattern = ".tsv",
  full.names = TRUE
)

# Read, bind, and add filename column
prokka_combined_df <- prokka_csv_files %>%
  map_df(~ read_tsv(.x, show_col_types = FALSE) %>%
           mutate(source_prokka = basename(.x))) %>%
  rename("Gene_ID" = "locus_tag")

#=================================================================================
#MEROPS FILES
#=================================================================================


# List all CSV files in the folder
merops_csv_files <- list.files(
  path = file.path("folder/merops_results"),
  pattern = ".csv",
  full.names = TRUE
)

# Read, bind, and add filename column
merops_combined_df <- merops_csv_files %>%
  map_df(~ read_csv(.x, show_col_types = FALSE) %>%
           mutate(source_merops = basename(.x))) %>%
  rename("Gene_ID" = "Query_ID")


#=================================================================================
#DBCAN FILES
#=================================================================================

# List all CSV files in the folder
dbcan_csv_files <- list.files(
  path = file.path("folder/dbcan_results"),
  pattern = ".tsv",
  full.names = TRUE
)

# Read, bind, and add filename column
dbcan_combined_df <- dbcan_csv_files %>%
  map_df(~ read_tsv(.x, show_col_types = FALSE) %>%
           mutate(source_dbcan = basename(.x))) %>%
  rename("Gene_ID" = "Gene ID")

##=================================================================================


# 1. Original number of unique Gene_IDs
total_gene_ids_original <- dbcan_combined_df %>%
  distinct(Gene_ID) %>%            # keep only unique IDs
  nrow()

# 2. Filter rows with #ofTools >= 2
filtered_df <- dbcan_combined_df %>%
  filter(`#ofTools` >= 2)

# 3. Number of unique Gene_IDs after filtering
total_gene_ids_filtered <- filtered_df %>%
  distinct(Gene_ID) %>% nrow()

# 4. How many unique Gene_IDs are removed
gene_ids_dropped <- total_gene_ids_original - total_gene_ids_filtered


cat("Total unique Gene_IDs before filtering:", total_gene_ids_original, "\n")
cat("Total unique Gene_IDs after filtering (#ofTools >= 2):", total_gene_ids_filtered, "\n")
cat("Number of unique Gene_IDs dropped:", gene_ids_dropped, "\n")

#=================================================================================
#EGGNOG FILES
#=================================================================================

# List all CSV files in the folder
eggnog_csv_files <- list.files(
  path = file.path("folder/eggnog_results"),
  pattern = ".annotations",
  full.names = TRUE
)

# Read, bind, and add filename column
eggnog_combined_df <- eggnog_csv_files %>%
  map_df(~ read_tsv(.x,
                    skip=4,
                    show_col_types = FALSE) %>%
           mutate(source_eggnog = basename(.x))) %>%
  rename("Gene_ID" = "#query")

#=================================================================================
#SIGNALP FILES
#=================================================================================

# List all CSV files in the folder
signalp_csv_files <- list.files(
  path = file.path("folder/signalp"),
  pattern = ".txt",
  full.names = TRUE
)

# Read, bind, and add filename column
signalp_combined_df <- signalp_csv_files %>%
  map_df(~ read_tsv(.x,skip = 1, show_col_types = FALSE) %>%
           mutate(source_signalp = basename(.x))) %>%
  rename("Gene_ID" = `# ID`) %>%
  separate(
    Gene_ID,
    into = c("Gene_ID", "description"),
    sep = " ",
    extra = "merge",
    fill = "right"
  )


#=================================================================================
#Bind FILES
#=================================================================================

bind_annotations <- prokka_combined_df %>%
  rename_with(~ paste0("prokka_", .x), -Gene_ID) %>%

  full_join(eggnog_combined_df %>%
  rename_with(~ paste0("eggnog_", .x), -Gene_ID)
            
            , by = "Gene_ID")  %>%
  
  full_join(merops_combined_df %>%
  rename_with(~ paste0("merops_", .x), -Gene_ID)
            
            , by = "Gene_ID") %>%
  
  full_join(dbcan_combined_df %>%
  rename_with(~ paste0("dbcan_", .x), -Gene_ID)
            
            , by = "Gene_ID")  %>%
  
  full_join(signalp_combined_df %>%
  rename_with(~ paste0("signalp_", .x), -Gene_ID)
            
            , by = "Gene_ID") 




bind_annotations_essentials <- bind_annotations %>%
  select(
    Gene_ID,
    prokka_ftype,
    prokka_EC_number,
    # prokka_COG,
    prokka_product,

    eggnog_evalue,
    eggnog_score,
    # eggnog_eggNOG_OGs,
    eggnog_COG_category,
    eggnog_Description,
    eggnog_Preferred_name,
    # eggnog_GOs,
    eggnog_EC,
    # eggnog_KEGG_ko,
    # eggnog_KEGG_Pathway,
    eggnog_CAZy,
    # eggnog_PFAMs,

    merops_MEROPS_ID,
    `merops_E-value`,
    merops_Score,
    merops_Description,
    merops_Peptidase_Family,
    merops_Peptidase_Class,

    `dbcan_EC#`,
    dbcan_dbCAN_hmm,
    dbcan_dbCAN_sub,
    dbcan_DIAMOND,
    `dbcan_#ofTools`,
    `dbcan_Recommend Results`,

    signalp_description,
    signalp_Prediction,
    `signalp_CS Position`
  )


```



editing bind data to format ec number columns
```{r}

combine_unique_ecs <- function(db, eg, pr) {
  ## 1️⃣  Guard against missing/empty values
  db  <- if (is.na(db)  || db  == "") NA_character_ else db
  eg  <- if (is.na(eg)  || eg  == "") NA_character_ else eg
  pr  <- if (is.na(pr)  || pr  == "") NA_character_ else pr


  ## 2️⃣  Clean the dbCAN side
  db_clean <- str_remove_all(db, ":\\d+")  %>%      # remove “:number”
              str_replace_all("\\|", ",") %>%      # unify delimiter
              str_remove_all("-")                 # drop stray hyphen(s)

  ## 3️⃣  Clean the eggNOG side
  eg_clean <- str_remove_all(eg, "-")              # also drop hyphens

  ## 4️⃣  Split into vectors
  db_vec <- str_split(db_clean, ",")[[1]]
  eg_vec <- str_split(eg_clean, ",")[[1]]
  pr_vec <- str_split(pr, ",")[[1]]

    ## 5️⃣  Combine and drop empty values
  ecs <- unique(c(db_vec, eg_vec, pr_vec))
  ecs <- ecs[!is.na(ecs) & ecs != ""]

  ## 6️⃣  Drop anything that STILL contains a wildcard (hyphen)
  # ecs <- ecs[!str_detect(ecs, "-")]

  if (length(ecs) == 0) return(NA_character_)
  paste(sort(ecs), collapse = ",")
}



#==============================================================================
# utility: clean any vector of raw strings, return a comma‑separated string
clean_families <- function(raw) {
  raw %>%
    # split on +, |, or comma
    str_split("[+|,]") %>% unlist() %>%
    str_trim() %>%
    # strip everything after "_" and parenthesis, keep only alphanum
    str_remove_all("_[^\\(]*") %>%
    str_remove_all("\\(.*\\)") %>%
    str_remove_all("[^A-Za-z0-9]") %>%
    unique() %>% sort() %>%
    # keep only actual CAZyme family names
    keep(~ grepl("^(GH|GT|CE|PL|CBM|AA)[0-9]+$", .)) %>%
    # collapse back to a single string (empty string if nothing left)
    paste(collapse = ",")
}

# 2. The main mutate() that builds `dbCAN_Unique`
# ------------------------------------------------------------------
cazy_proc_better <- bind_annotations_essentials %>%
  rowwise() %>%                                 # per‑row processing
  mutate(
    EC_combined = combine_unique_ecs(`dbcan_EC#`, eggnog_EC, prokka_EC_number),

    dbCAN_Unique = {
      # 2.1  Prefer the priority column – clean it first
      rec_raw <- `dbcan_Recommend Results`

      if (!is.na(rec_raw) && rec_raw != "") {
        rec_clean <- clean_families(rec_raw)
        if (rec_clean != "") rec_clean else {
          # no valid family was found – fall back
          src_cols <- c(`dbcan_dbCAN_hmm`, `dbcan_dbCAN_sub`, `dbcan_DIAMOND`)
          src_cols <- src_cols[!is.na(src_cols) & src_cols != ""]
          if (length(src_cols) == 0) NA_character_
          else clean_families(src_cols)
        }
      } else {
        # 2.2  No value in the priority column –  all in any of the tool columns
        src_cols <- c(`dbcan_dbCAN_hmm`, `dbcan_dbCAN_sub`, `dbcan_DIAMOND`)
        src_cols <- src_cols[!is.na(src_cols) & src_cols != ""]
        if (length(src_cols) == 0) NA_character_
        else clean_families(src_cols)
      }
    }
  ) %>% ungroup()



#keep every row, but blank out MEROPS fields when ANY of the following conditions would have removed the row: E‑value > 0.01; Non‑peptidase homologues (MEROPS family contains .UN…); Inhibitors (MEROPS family begins with I, e.g., I02.A01)
#want to set both: merops_MEROPS_ID = NA and merops_Peptidase_Family = NA

temp_merops_clean_cazy_done <- cazy_proc_better %>%
  mutate(
    # condition 1: bad E-value
    bad_evalue = `merops_E-value` > 0.01,

    # condition 2: non-peptidase homologues (UN)
    is_homologue_UN = grepl("\\.[^.]*UN", merops_Peptidase_Family, ignore.case = TRUE),

    # condition 3: inhibitors (MEROPS families starting with I)
    is_inhibitor = grepl("^I", merops_Peptidase_Family, ignore.case = TRUE),
    
    # condition 4: dbCAN says 3 tools hit
    is_dbcan_three     = replace_na(`dbcan_#ofTools` == 3, FALSE),
    # apply NA to MEROPS fields if ANY condition is true
    merops_MEROPS_ID = if_else(
      bad_evalue | is_homologue_UN | is_inhibitor | is_dbcan_three,
      NA_character_,
      merops_MEROPS_ID
    ),

    merops_Peptidase_Family = if_else(
      bad_evalue | is_homologue_UN | is_inhibitor | is_dbcan_three,
      NA_character_,
      merops_Peptidase_Family
    ),

    merops_Peptidase_Class = if_else(
      bad_evalue | is_homologue_UN | is_inhibitor | is_dbcan_three,
      NA_character_,
      merops_Peptidase_Class
  )) 


full_raw_annotations_all_species <-  temp_merops_clean_cazy_done %>%
  separate_rows(EC_combined,
                sep = ",") %>%
  filter(!str_starts(trimws(as.character(Gene_ID)), "#")) %>%
  mutate(
    family_only_merops   = ifelse(is.na(merops_Peptidase_Family) | merops_Peptidase_Family == "",
                             NA_character_,
                             sub("\\..*", "", merops_Peptidase_Family)),   # everything before the dot
    subfamily_only_merops = ifelse(is.na(merops_Peptidase_Family) | !grepl("\\.", merops_Peptidase_Family),
                              NA_character_,
                              sub(".*\\.", "", merops_Peptidase_Family))    # everything after the dot
  )

# write_csv(full_raw_annotations_all_species, file.path( ONE,"folder/full_raw_annotations_all_species.csv"))  

```



#Mapping by EC number

```{r}

collapse_substrates_db_new <- db_ECM_BRENDA_revised %>%
  group_by(EC_number) %>%
  reframe(
    across(-matched_term, ~ first(.x)),
    matched_term = paste(matched_term, collapse = "; ")
  )


# ----  Match by EC (any  EC to any ECM EC) ----------------------------

new_matched_by_ec1 <- full_raw_annotations_all_species %>%
  # rename_with(~ paste0(.x, "_annotation"), -any_of(c("Strain", "Gene_ID", "EC_combined"))) %>%
  left_join(collapse_substrates_db_new %>%
             select(-synonyms) , 
            by = c("EC_combined" = "EC_number")) %>%         # join on the EC number

  rename("EC_number" = "EC_combined" ) %>%
  distinct()




```




#Mapping by MEROPS db
```{r}

merops_family_match <- inner_join(collapsed_merops,
                                  full_raw_annotations_all_species %>%
                                    select(
                                        Gene_ID,
                                        family_only_merops,
                                        subfamily_only_merops,
                                        merops_MEROPS_ID,
                                        `merops_E-value`,
                                        EC_combined
                                        
                                        ),
                                  by = c("family_merops" = "family_only_merops"))

ec_merops_family <- left_join(new_matched_by_ec1,
                              collapsed_merops
                                  ,
                                  by = c("family_only_merops" = "family_merops" )) 

               

ec_merops_family1 <- ec_merops_family %>%
  mutate(
  matched_term_joined = map2_chr(`matched_term.x`, 
                                 `matched_term.y`, ~ {
      # helper: split a cell, trim spaces, drop empties
      split_clean <- function(x) {
        if (is.na(x) | x == "") character(0)
        else  str_split(x, "[,;]")[[1]] %>% 
              str_trim() %>%
              .[ . != "" ]               # remove empty strings
      }

      parts_a <- split_clean(.x)
      parts_b <- split_clean(.y)

      # union of the two lists, sorted for consistency
      union <- unique(c(parts_a, parts_b))
      if (length(union) == 0) NA_character_ else paste(sort(union), collapse = ";")
    })
  ) %>%
  mutate(
    # First decide which rows have an “empty” Strain.
    # Adjust the test if you use NA, an empty string, or both.
    Strain = case_when(
      # 
      startsWith(Gene_ID, "CJALCLJH") ~ "B. acidifaciens",
      startsWith(Gene_ID, "GIPNGPCC") ~ "B. caccae",
       startsWith(Gene_ID, "OEKDHLHL") ~ "B. cellulolyticus",
       startsWith(Gene_ID, "LIKKJHMI") ~ "B. fragilis",
       startsWith(Gene_ID, "MMBIAKFC") ~ "B. intestinalis",
       startsWith(Gene_ID, "IHKDOBEK") ~ "B. ovatus",
       startsWith(Gene_ID, "OCFCDACC") ~ "B. thetaiotaomicron",
       startsWith(Gene_ID, "OIGBMDAE") ~ "B. uniformis",
       startsWith(Gene_ID, "JAPKKLMI") ~ "B. stercoris",
       startsWith(Gene_ID, "IABPPNAN") ~ "B. salyersiae",
       startsWith(Gene_ID, "EHKGMGKO") ~ "B. xylanisolvens"
      
     
    )
  ) 


```


  
```{r}  
  

# 1.  Synonym map for substrates that must be merged
# ------------------------------------------------------------------
substrate_map <- list(
  "hyaluron"          = "hyaluronan",
  "hyaluronic acid"   = "hyaluronan"
  # add more pairs here in the same pattern
)

## helper that normalises a single term:
## trim blanks, lower‑case, then replace by the canonical value if present
norm_sub <- function(term) {
  t <- str_trim(term)          # remove stray spaces
  l <- tolower(t)              # lower‑case for matching
  if (l %in% names(substrate_map)) substrate_map[[l]] else t
}

# ------------------------------------------------------------------
# 2.  Gene‑level summarisation
# ------------------------------------------------------------------
ecm_gene_summary <- ec_merops_family2 %>%
  group_by(Gene_ID) %>%                 # one row per gene

  summarise(
    # # 2a.  Collapsed EC numbers
    # EC_number = paste(
    #   sort(unique(EC_number[!is.na(EC_number)])),
    #   collapse = ","
    # ),

    # 2b.  Unique, normalised substrates
    matched_term_joined = {
      # break into individual substrings, flatten, drop NAs
      subs_vec <- 
        str_split(matched_term_joined[!is.na(matched_term_joined)], ";" ) %>%


        unlist()

      # normalise every substring and keep only one copy of each
      subs_clean <- unique(
        vapply(subs_vec, norm_sub, character(1), USE.NAMES = FALSE)
      )

      # collapse back into a semi‑colon string
      paste(subs_clean, collapse = ";")
    },

    # 2c.  All remaining columns – first non‑missing value from the group
    across(
      .cols = -c(matched_term_joined),
      .fns  = ~ first(.x[!is.na(.x)]),
      .names = "{.col}"
    ),

    .groups = "drop"
  ) %>%
  mutate(
    cazy_class = str_remove_all(dbCAN_Unique, "\\d")) %>%
# 5.  Append "_sulfatase" if either text column mentions it

  mutate(
    sulfatases = 
      (str_detect(prokka_product,     regex("sulfatase", ignore_case = TRUE)) |
       str_detect(eggnog_Description, regex("sulfatase", ignore_case = TRUE)))  
    ) 





ecm_enriched <- ecm_gene_summary %>%
  filter(!(is.na(matched_term_joined) | matched_term_joined == ""))
  
# write_csv(ecm_enriched, file.path( ONE,"folder/new_ecm_mapped_combined.csv")) 


```


```{r}

# helper – split a comma‑separated string into a *clean* vector
split_clean <- function(x) {
  if (is.na(x) || x == "") character(0) else {
    x %>%                     # split, trim, drop empties
      str_split(",") %>% unlist() %>%
      str_trim() %>% .[ . != "" ]
  }
}

# ---------------------------------------------------------------
#   Summaries per strain (wide format)
# ---------------------------------------------------------------
strain_summary_wide <- ecm_gene_summary %>%
  
  # 1. Mark rows that actually have each annotation
  mutate(
    has_ecm      = matched_term_joined != c("", NA),
    has_prot     = !is.na(family_only_merops) &
                   family_only_merops != "" & 
      matched_term_joined != "",
    # ----- 3. CAZymes -------------------------------------------------
    #  - the usual rule: non‑empty dbCAN_Unique AND a substrate is listed
    #  - plus: the substrList has one of the four special words
    has_cazy = (
  
        
        !is.na(dbCAN_Unique) &
       dbCAN_Unique != "" &
       matched_term_joined != "") 
        
    
  ) %>%
  
  # 2. Group by strain and aggregate
  group_by(Strain) %>%
  summarise(

    # # ECMs ----------------------------------------------------------
    # ecm_genes = n_distinct(Gene_ID[has_ecm]),

    # Proteases ------------------------------------------------------
    protease_genes = n_distinct(Gene_ID[has_prot]),


    # CAZymes -------------------------------------------------------
    cazyme_genes = n_distinct(Gene_ID[has_cazy]),

    .groups = "drop"                            # final flat data frame
  )

# ---------------------------------------------------------------
#   classes
# ---------------------------------------------------------------
strain_summary <- ecm_gene_summary %>%
  
  # 1. Mark rows that actually have each annotation
  mutate(
    has_ecm      = matched_term_joined != c("", NA),
    has_prot     = !is.na(merops_Peptidase_Family) &
                   merops_Peptidase_Family != "" & 
      matched_term_joined != "",
    # ----- 3. CAZymes -------------------------------------------------
    #  - the usual rule: non‑empty dbCAN_Unique AND a substrate is listed

    has_cazy = (
        
        !is.na(dbCAN_Unique) &
       dbCAN_Unique != "" &
       matched_term_joined != "")   )) %>%
  separate_rows(dbCAN_Unique, sep = ";") %>%   # explode the list
  separate_rows(dbCAN_Unique, sep = ",") %>%
  separate_rows(cazy_class, sep = ";")%>%
  separate_rows(cazy_class, sep = ",")
  
  
ecm_classes_caz <- strain_summary %>%
  filter(has_cazy == TRUE ) %>%
  
    
  group_by(Strain, cazy_class) %>%
  
  summarise(
    gene_count = n_distinct(Gene_ID)
  ) %>%
  pivot_wider(
    names_from = cazy_class,
    values_from = gene_count
  )


ecm_classes_prot <- strain_summary %>%
  mutate(merops_Peptidase_Class = substr(merops_Peptidase_Class, 1, 1)) %>%
  filter(has_prot == TRUE ) %>%

  group_by(Strain, merops_Peptidase_Class) %>%
  
  summarise(
    gene_count = n_distinct(Gene_ID)
  ) %>%
  pivot_wider(
    names_from = merops_Peptidase_Class,
    values_from = gene_count
  )




# ---------------------------------------------------------------
#   Summaries per strain (wide format)
# ---------------------------------------------------------------
substrate_summary <- ecm_gene_summary %>%
  
  # Mark rows that actually have each annotation
  mutate(
    has_ecm      = matched_term_joined != c("", NA),
    
    
    
    has_prot     = !is.na(family_only_merops) &
                   family_only_merops != "" & 
      matched_term_joined != "" & 
      !is.na(matched_term_joined),
    #  CAZymes -------------------------------------------------
    #  - the usual rule: non‑empty dbCAN_Unique AND a substrate is listed
    has_cazy = (

        
        !is.na(dbCAN_Unique) &
       dbCAN_Unique != "" &
       matched_term_joined != "" & 
      !is.na(matched_term_joined)) 
        
    
  ) %>%
  filter( has_prot == TRUE | has_cazy == TRUE
    
    )



# ---------------------------------------------------------------
#   Substrate‑target table (per strain & per substrate)
# ---------------------------------------------------------------
substrate_targets <- substrate_summary %>%
  separate_rows(matched_term_joined, sep = ";") %>%   # explode the list
  separate_rows(matched_term_joined, sep = ",") %>%
  mutate(substrate = str_trim(matched_term_joined)) %>%
  filter(substrate != "" & !is.na(substrate)) %>%
  group_by(Strain, substrate) %>%
  summarise(
    ## list of Gene_IDs that contain this substrate
    unique_genes    = n_distinct(Gene_ID),            # count
    gene_ids        = list(unique(Gene_ID)),           # vector
    gene_ids_str    = paste(unique(Gene_ID), collapse = ", "),  # comma‑separated text
    
    .groups = "drop"
  ) %>%
  rename(Substrate = substrate)



substrates_gene_list <- substrate_summary %>%
  separate_rows(matched_term_joined, sep = ";") %>%   # explode the list
    separate_rows(matched_term_joined, sep = ",") %>%   # explode the list

  distinct() %>%
  filter(!(matched_term_joined == "" | is.na(matched_term_joined)))

substrates_gene_list_families <- substrates_gene_list %>%
  select(c(
  Gene_ID,               
matched_term_joined,
  EC_number,
  dbCAN_Unique,
cazy_class,
  family_only_merops,
  subfamily_only_merops,
merops_Peptidase_Class,
  Peptidase_name,
  Substrate_name,
Strain
)) %>%
  group_by(Strain, matched_term_joined, merops_Peptidase_Class) %>%
  summarise(
    number_genes = n_distinct(Gene_ID)
  ) 





```


#substrate specific list of enzymes, elastin, collagen 
```{r}

#full list
elastin_full_list <- susbstrates_gene_list %>%
  filter(matched_term_joined == "elastin")

#list of classes by strain
elastin_classes <- substrate_targets_without_cbm_gt %>%
  filter(matched_term_joined == "elastin") %>%
  group_by(Strain, merops_Peptidase_Class) %>%
  summarise(
    ## list of Gene_IDs that contain this substrate
    unique_genes    = n_distinct(Gene_ID),  # comma‑separated text
    
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = merops_Peptidase_Class,
              values_from = unique_genes) 

#full list
collagen_full_list <- susbstrates_gene_list %>%
  filter(matched_term_joined == "collagen") 

#list of classes by strain
collagen_classes <- substrate_targets_without_cbm_gt %>%
  filter(matched_term_joined == "collagen") %>%
  group_by(Strain, merops_Peptidase_Class) %>%
  summarise(
    ## list of Gene_IDs that contain this substrate
    unique_genes    = n_distinct(Gene_ID),  # comma‑separated text
    
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = merops_Peptidase_Class,
              values_from = unique_genes) 



```




#table for family to substrates relationship

```{r}

unique_terms_in_cell <- function(cell, sep = ";") {
  # 1.  Split on the column‑separator (e.g. ";")
  terms <- strsplit(cell, sep, fixed = TRUE)[[1]]

  # 2.  Trim leading / trailing whitespace
  terms <- trimws(terms)

  # 3.  Remove empty strings
  terms <- terms[terms != ""]

  # 4.  Keep only the distinct terms (in the order they first appear)
  unique_terms <- terms[!duplicated(terms)]

  # 5.  Re‑join them
  paste(unique_terms, collapse = sep)
}




table_families_to_EC_substrates <-   ecm_gene_summary %>%
separate_rows(dbCAN_Unique, sep = ",") %>%
  
  select(c(
    dbCAN_Unique,
    family_only_merops,
      EC_number,
      matched_term_joined,
      recommended_name,
      Peptidase_name,
      prokka_product,
      eggnog_Description,
      merops_Description,
      matched_term.x,
      matched_term.y,
      Substrate_name
    
  )) %>%
  
  distinct() %>%
  group_by(
    dbCAN_Unique,
    family_only_merops
  ) %>%
  summarise(
    EC_numbers =  paste(unique(EC_number), collapse = ", "),
    substrates_joined = paste(unique(matched_term_joined), collapse = ";"),
    BRENDA_recommended_name =  paste(unique(recommended_name), collapse = ", "),
    MEROPSdb_peptidase_name =  paste(unique(Peptidase_name), collapse = ", "),
      
      annotation_prokka_product =  paste(unique(prokka_product), collapse = ", "),
      annotation_eggnog_Description =  paste(unique(eggnog_Description), collapse = ", "),
      annotation_merops_Description =  paste(unique(merops_Description), collapse = ", "),
      BRENDA_matched_term.x =  paste(unique(matched_term.x), collapse = ", "),
      MEROPSdb_matched_term.y =  paste(unique(matched_term.y), collapse = ", "),
      MEROPSdb_exactSubstrate_name =  paste(unique(Substrate_name), collapse = ", ")
    ) %>%
   mutate(
    unique_substrates  = vapply(substrates_joined, unique_terms_in_cell, character(1))) %>%
  
  filter(!(substrates_joined == "")) %>%
  select(c(
    "CAZy family" = dbCAN_Unique,
    "MEROPS family" = family_only_merops,
    "Associated EC number(s)" = EC_numbers,
    "Mapped substrate(s)" =unique_substrates,
    "Enzyme name" = MEROPSdb_peptidase_name,
    BRENDA_recommended_name, 
    "Prokka annotation" =annotation_prokka_product,
    "eggNOG annotation" =annotation_eggnog_Description,
    "MEROPS annotation" =annotation_merops_Description))
    


```



#removing substrates from non-degradation classes: CBM and GT CAZy classes annotations 

```{r}


# ---------------------------------------------------------------
# Substrate‑target table (per strain & per substrate)
# ---------------------------------------------------------------
substrate_targets_without_cbm_gt <- substrate_summary %>%
  separate_rows(cazy_class, sep = ",") %>%
  
  mutate(
    matched_term_joined = case_when(
      # When the class is “CBM” or “GT” *and* we have no MEROPS family
      cazy_class %in% c("CBM", "GT") & is.na(family_only_merops) ~ NA_character_,  # clear it

      # otherwise keep the original value
      TRUE ~ matched_term_joined
    )
  ) %>%
  
  
  separate_rows(matched_term_joined, sep = ";") %>%   # explode the list
  separate_rows(matched_term_joined, sep = ",") %>%
  mutate(substrate = str_trim(matched_term_joined)) %>%
  filter(substrate != "" & !is.na(substrate)) %>%
  group_by(Strain, substrate) %>%
  summarise(
    ## list of Gene_IDs that contain this substrate
    unique_genes    = n_distinct(Gene_ID),            # count
    gene_ids        = list(unique(Gene_ID)),           # vector
    gene_ids_str    = paste(unique(Gene_ID), collapse = ", "),  # comma‑separated text
    
    .groups = "drop"
  ) %>%
  rename(Substrate = substrate)

```



#Summary annotation families (Supplementary table 2)
```{r}

supplement_table2 <-  table_families_to_EC_substrates %>%
  distinct() %>%
  separate_rows(`Associated EC number(s)`,
                `Mapped substrate(s)`,
                `Enzyme name`,
                sep = "[,;]",
                trim_ws = TRUE)


                    
       # Define which columns contain comma or semicolon separated values
cols_to_clean <- c(
  "Associated EC number(s)",
  "Mapped substrate(s)",
  "Enzyme name",
  "BRENDA_recommended_name",
  "Prokka annotation",
  "eggNOG annotation",
  "MEROPS annotation"
)

# Function to replace separators with newlines
clean_cell <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    TRUE ~ x %>%
      str_replace_all(",\\s*", "\n") %>%
      str_replace_all(";\\s*", "\n") %>%
      str_trim()
  )
}

# Apply cleaning to target columns
supplement_table2_clean <- supplement_table2 %>%
  ungroup() %>%
  mutate(across(all_of(cols_to_clean), clean_cell))

# Create workbook with wrap text style
wb <- createWorkbook()
addWorksheet(wb, "Supplementary Table 2")

# Write data
writeData(wb, sheet = 1, x = supplement_table2_clean)

# Apply wrap text style to entire table
wrap_style <- createStyle(wrapText = TRUE, valign = "top")
addStyle(wb, sheet = 1,
         style = wrap_style,
         rows = 1:(nrow(supplement_table2_clean) + 1),
         cols = 1:ncol(supplement_table2_clean),
         gridExpand = TRUE)             
                    


clean_brenda <- function(x) {
  if (is.na(x)) return(NA_character_)
  x %>%
    str_split("\n") %>%                        # split into individual lines
    unlist() %>%
    str_remove("^RN\\t") %>%                   # remove "RN\t" at start of line
    str_trim() %>%                             # trim whitespace
    .[. != "NA" & . != ""] %>%                 # remove lines that are just "NA" or empty
    paste(collapse = "\n")                     # rejoin with newline
}

supplement_table2_clean2 <- supplement_table2_clean %>%
  mutate(BRENDA_recommended_name = sapply(BRENDA_recommended_name, clean_brenda)) 



supplement_table2_clean22 <- supplement_table2_clean2 %>%
   # Step 1: Merge CAZy family and MEROPS family into one column
  mutate(
    Family = case_when(
      !is.na(`CAZy family`) & !is.na(`MEROPS family`) ~ paste(`CAZy family`, `MEROPS family`, sep = "\n"),
      !is.na(`CAZy family`) ~ `CAZy family`,
      !is.na(`MEROPS family`) ~ `MEROPS family`,
      TRUE ~ NA_character_
    )
  ) %>%
  # Step 2: Drop the original two columns
  select(-`CAZy family`, -`MEROPS family`) %>%
  # Step 3: Rename remaining columns
  
  # Step 1: Rename columns that already exist
  rename(
    "Family"                          = Family,
    "EC number(s)"                    = `Associated EC number(s)`,
        "MEROPS recommended enzyme name(s)" =    "Enzyme name",
    "BRENDA recommended enzyme name(s)" = BRENDA_recommended_name,
    "ECM substrate(s)"                = `Mapped substrate(s)`,
    "Prokka annotation(s)"            = `Prokka annotation`,
    "eggNOG functional annotation(s)" = `eggNOG annotation`,
    "MEROPS annotation(s)"            = `MEROPS annotation`
  ) %>%
  # Step 2: Add missing columns as placeholders
  mutate(
    "Number of species it appears in" = NA_real_,
    "Enzyme class"                    = NA_character_,
    "MEROPS recommended enzyme name(s)" = NA_real_,   # or NA if not derivable
  ) %>%
  # Step 3: Reorder columns to match desired order
  select(
    "Family",
    "Number of species it appears in",
    "Enzyme class",
    "EC number(s)",
    "BRENDA recommended enzyme name(s)",
    "MEROPS recommended enzyme name(s)",
    "ECM substrate(s)",
    "Prokka annotation(s)",
    "eggNOG functional annotation(s)",
    "MEROPS annotation(s)"
  ) %>%
  mutate(`Enzyme class` = sapply(Family, clean_enzyme_class))



clean_enzyme_class <- function(x) {
  if (is.na(x)) return(NA_character_)
  
  # Split by newline in case two families are merged
  families <- str_split(x, "\n")[[1]] %>% str_trim()
  
  classes <- case_when(
    # CAZyme classes first (more specific, must come before single letters)
    str_detect(families, "(?i)^GH")  ~ "Glycoside hydrolase",
    str_detect(families, "(?i)^GT")  ~ "Glycosyltransferase",
    str_detect(families, "(?i)^PL")  ~ "Polysaccharide lyase",
    str_detect(families, "(?i)^CE")  ~ "Carbohydrate esterase",
    str_detect(families, "(?i)^AA")  ~ "Auxiliary activity",
    str_detect(families, "(?i)^CBM") ~ "Carbohydrate-binding module",
    # MEROPS protease classes (single letter prefix)
    str_detect(families, "(?i)^M")   ~ "Metalloprotease",
    str_detect(families, "(?i)^C")   ~ "Cysteine protease",
    str_detect(families, "(?i)^S")   ~ "Serine protease",
    str_detect(families, "(?i)^A")   ~ "Aspartic protease",
    str_detect(families, "(?i)^T")   ~ "Threonine protease",
    str_detect(families, "(?i)^G")   ~ "Glutamic protease",
    str_detect(families, "(?i)^P")   ~ "Mixed protease",
    TRUE ~ NA_character_
  )
  
  # Remove NAs and duplicates, then join with newline
  classes %>%
    .[!is.na(.)] %>%
    unique() %>%
    paste(collapse = "\n")
}






write.csv(supplement_table2_clean22 %>% ungroup() %>% as.data.frame() %>% `rownames<-`(NULL),
          file = file.path(ONE, "folder/Supplementary_Table2_final.csv"),
          row.names = FALSE)



```

