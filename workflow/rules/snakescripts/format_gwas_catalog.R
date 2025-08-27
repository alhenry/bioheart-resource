# format GWAS catalog data for use in the workflow

library(data.table)
library(tidyverse)
library(glue)

INPUT <- snakemake@input
PARAMS <- snakemake@params
OUTPUT <- snakemake@output

# INPUT <- list(
#     gwas_file = "resources/gwas/gwas_catalog/sbp/GCST90310294.h.tsv.gz",
#     frq_files = glue("resources/genotypes/bioheart975/chr{i}.frq", i = 1:22)
# )

# PARAMS <- list(
#     max_maf_diff = 0.2
# )

df_gwas <- fread(INPUT$gwas_file)

df_frq <- map_df(INPUT$frq_files, fread)

df_gwas[, maf := pmin(effect_allele_frequency, 1 - effect_allele_frequency)]

df_gwas_format <- df_gwas %>% 
    mutate(SNP = str_replace_all(variant_id, "_", ":")) %>% 
    inner_join(df_frq[,.(SNP, MAF)]) %>% 
    filter(abs(maf - MAF) <= PARAMS$max_maf_diff) %>% 
    select(SNP, CHR = chromosome, BP = base_pair_location, P = p_value, N = n)


fs::dir_create(dirname(OUTPUT$formatted_gwas))
fwrite(df_gwas_format, OUTPUT$formatted_gwas, sep = "\t", quote = FALSE)