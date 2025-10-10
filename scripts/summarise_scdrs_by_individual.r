# get summary 
library(data.table)
library(arrow)
library(tidyverse)

# load annotation
df_annot <- fread("resources/cell_type_annot/scClassify_asian_level3_final.csv",
                  key = "cellname")
df_cell_score <- read_parquet("results/aggregate/bioheart975-full.scdrs.cell_score.tsv.parquet.gz")
df_cell_mcp <- read_parquet("results/aggregate/bioheart975-full.scdrs.cell_mcp.tsv.parquet.gz")

setDT(df_cell_score, key = "index")
setDT(df_cell_mcp, key = "index")

# update cell type annotation
df_cell_score[df_annot, cell_type := i.celltype]
df_cell_mcp[df_annot, cell_type := i.celltype]

# calculate proportion of significantly enriched cell
# significantly enriched cell -> mcp < 0.05
traits <- c("dbp", "sbp", "primary_ht")
df_prop_sig <- df_cell_mcp %>% 
  pivot_longer(all_of(traits), names_to = "trait") %>% 
  group_by(trait) %>% 
  summarise(prop_sig = mean(value < 0.05),
            n_cells = n())