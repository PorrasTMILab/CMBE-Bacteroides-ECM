
#libraries
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
library(pheatmap)
library(ggdendro)
library(cowplot)
library(stringr)
library(ggrepel)
library(ggrounded)
library(forcats)
library(grid) 
library(xfun)
library(ggtext)
library(fillpattern)
library(ggpattern)
library(tidyverse)
library(ggalluvial)
library(viridis)
library(igraph)
library(ggraph)
library(scales)
library(magick)
library(schrute)
library(treemapify)
library(RColorBrewer)
library(tidytext)
library(vegan)     # Bray‑Curtis, hclust…
library(Hmisc)
library(vegan)
library(openxlsx)   # install.packages("tidytext")
#devtools::install_github("vqf/nVennR2")

```

#Dendrogram building
```{r}
##  Build the raw count matrix
count_mat_raw <- substrate_targets %>%
  filter(Substrate %in% intestine_subs) %>%   # keep only the substrates we care about
  select(Strain, Substrate, unique_genes) %>%      # pick raw gene count (unique genes)
  pivot_wider(
    names_from   = Substrate,
    values_from  = unique_genes,
    values_fill  = list(unique_genes = 0)          # 0 when a strain has no genes for that substrate
  ) %>%
  column_to_rownames("Strain")                  # rows = species, cols = substrates

##  log10(x + 1) – keep the transform in a separate matrix
count_mat_log <- log10(count_mat_raw + 1)



##   Bray–Curtis distance on the log‑transformed data
# vegan::vegdist expects a data frame or matrix with rows = samples
bray_dist <- vegdist(count_mat_log, method = "bray")

# convert to a matrix for later use or inspection
bray_dist_mat <- as.matrix(bray_dist)
rownames(bray_dist_mat) <- rownames(count_mat_log)
colnames(bray_dist_mat) <- rownames(count_mat_log)

hc_upgma   <- hclust(bray_dist, method = "average")

```



#Run Spearman Correlations (Per Substrate)
```{r}

# Combine  data frames column-wise 
# We want per-substrate correlations between dfA and dfB columns

substrates <- colnames(dfA)
rho <- numeric(length(substrates)); pval <- numeric(length(substrates))
names(rho) <- substrates; names(pval) <- substrates

for (sub in substrates) {
  x <- as.numeric(as.character(dfA[[sub]]))
  y <- as.numeric(as.character(dfB[[sub]]))
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) {
    rho[sub] <- NA_real_; pval[sub] <- NA_real_
    next
  }
  tmp <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  rho[sub] <- as.numeric(tmp$estimate)
  pval[sub] <- tmp$p.value
}

result_table <- data.frame(
  Substrate = substrates,
  Spearman_rho = rho,
  P_value = pval
)
print(result_table)

# FDR correction
result_table$P_adj <- p.adjust(result_table$P_value, method = "fdr")

```




#Mantel test
```{r}
# For matrixA
# Extract species and numeric substrate data separately
species <- matrixA$Strain
Xa <- matrixA %>% select(-Strain)  # remove species name column

# Ensure numeric (this should be automatic if data are numeric)
Xa_numeric <- as.data.frame(lapply(Xa, as.numeric))

# Now scale
Xa_z <- scale(Xa_numeric)

# Reassemble with species
dfA_z <- as.data.frame(Xa_z) %>%
  mutate(Strain = species) %>%
  select(Strain, everything())

#=======================================================
# Repeat for matrixB
species_b <- matrixB$Strain
Xb <- matrixB %>% select(-Strain)
Xb_numeric <- as.data.frame(lapply(Xb, as.numeric))
Xb_z <- scale(Xb_numeric)
dfB_z <- as.data.frame(Xb_z) %>%
  mutate(Strain = species_b) %>%
  select(Strain, everything())



# Define substrate groups  
gag_substrates <- c("chondroitin", "heparan", "hyaluronan")
protein_substrates <- c("collagen", "fibronectin", "elastin")

# Helper function: z-score and compute Euclidean distance matrix
compute_dist <- function(df, substrates) {
  # Extract numeric columns by substrate name
  numeric_data <- df %>%
    select(all_of(substrates)) %>%
    mutate_all(as.numeric)
  
  # Z-score columns
  z_data <- scale(numeric_data)
  
  # Euclidean distance matrix (rows = species)
  dist_matrix <- dist(z_data, method = "euclidean")
  
  return(dist_matrix)
}

# Mantel test function
run_mantel <- function(dist1, dist2, permutations = 9999) {
  mantel_res <- mantel(dist1, dist2, method = "spearman", permutations = permutations)
  tibble(
    Mantel_r = mantel_res$statistic,
    p_value = mantel_res$signif,
    permutations = permutations
  )
}

# Prepare data frames - remove Strain column for calculations
dfA_num <- matrixA %>% select(-Strain)
dfB_num <- matrixB %>% select(-Strain)

# All substrates Mantel
distA_all <- compute_dist(matrixA, colnames(dfA_num))
distB_all <- compute_dist(matrixB, colnames(dfB_num))
mantel_all <- run_mantel(distA_all, distB_all)

# GAG-only Mantel
distA_gag <- compute_dist(matrixA, gag_substrates)
distB_gag <- compute_dist(matrixB, gag_substrates)
mantel_gag <- run_mantel(distA_gag, distB_gag)

# Protein-only Mantel
distA_protein <- compute_dist(matrixA, protein_substrates)
distB_protein <- compute_dist(matrixB, protein_substrates)
mantel_protein <- run_mantel(distA_protein, distB_protein)

# Combine & label results
results <- bind_rows(
  mantel_all %>% mutate(Test = "Mantel (All substrates)"),
  mantel_gag %>% mutate(Test = "Mantel (GAG-only)"),
  mantel_protein %>% mutate(Test = "Mantel (Protein-only)")
) %>% select(Test, everything())

print(results)

```
