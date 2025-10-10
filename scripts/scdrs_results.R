# scdrs exploration

library(data.table)
library(tidyverse)
library(scales)
library(patchwork)
library(arrow)
library(ragg)
library(paletteer)
library(ggforce)
library(ggridges)
library(ggrepel)
library(readxl)
library(glue)
library(geomtextpath)
library(yaml)

# df_cell_map <- fread("resources/misc/cell_map.tsv") %>% 
#   mutate(cell_type = factor(cell_type, unique(cell_type)))
# 
# df_trait_map <- read_excel("resources/metadata/trait_metadata_curated.xlsx") %>% 
#   filter(include) %>% 
#   setDT()
# df_trait_cat <- read_excel("resources/metadata/trait_metadata_curated.xlsx",
#                            sheet = "trait_category_order") %>% 
#   setDT() 

df_trait_map <- fs::dir_ls("resources/metadata/phenotypes") %>% 
  map_df(~yaml::read_yaml(.x))

df_stats <- fread("results/aggregate/bioheart975-full.scdrs.cell_type_stats.tsv")
df_top <- fread("results/aggregate/bioheart975-full.scdrs.cell_type_top.tsv")

df_cell_map <- fread("resources/metadata/cell_map.tsv") %>% 
  mutate(cell_type = factor(cell_type, unique(cell_type)))

df_cells <- read_parquet("results/aggregate/bioheart975-full.scdrs.cell_score.tsv.parquet.gz")

# get cell mcp value
df_mcp <- read_parquet("results/aggregate/bioheart975-full.scdrs.cell_mcp.tsv.parquet.gz") %>% 
  pivot_longer(any_of(df_trait_map$trait_id), names_to = "phenotype", values_to = "assoc_mcp") %>% 
  setDT()

# n_cells <- nrow(df_cells)

# check distribution
# ggplot(df_mcp) +
#   geom_histogram(aes(x = assoc_mcp), bins = 100, fill = "gray", color = "black") +
#   geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
#   theme_bw()
# 
# df_mcp %>%
#   filter(phenotype == "sbp") %>%
#   ggplot() +
#   geom_histogram(aes(x = assoc_mcp), bins = 50, fill = "gray", color = "black", linewidth = 0.1) +
#   theme_bw() +
#   geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
#   facet_wrap(~cell_type)

setDT(df_cells)

# df_cells[df_cell_map,
#          `:=`(cell_label = factor(i.cell_type, df_cell_map$cell_type),
#               major_cell_type = factor(i.major_cell_type, unique(df_cell_map$major_cell_type))),
#          on = c(cell_type = "wg2_scpred_prediction")]
# df_stats[df_cell_map,
#          `:=`(cell_label = factor(i.cell_type, df_cell_map$cell_type),
#               major_cell_type = factor(i.major_cell_type, unique(df_cell_map$major_cell_type))),
#          on = c(cell_type = "wg2_scpred_prediction")]
# df_mcp[df_cell_map,
#        `:=`(cell_label = factor(i.cell_type, df_cell_map$cell_type),
#             major_cell_type = factor(i.major_cell_type, unique(df_cell_map$major_cell_type))),
#        on = c(cell_type = "wg2_scpred_prediction")]

setDT(df_trait_map)
df_stats[df_trait_map, `:=`(pheno_label = i.trait_name),
         on = c(phenotype = "trait_id")]
df_mcp[df_trait_map, `:=`(pheno_label = i.trait_name),
       on = c(phenotype = "trait_id")]

df_mcp_sig <- df_mcp[assoc_mcp < 0.05] 

df_mcp_sig[df_stats, `:=`(ct_sig = i.assoc_mcp < 0.05,
                          n_cell_ct = i.n_cell),
           on = c("phenotype", "cell_type")]

df_mcp_sig_count <- df_mcp_sig %>% 
  group_by(cell_type, ct_sig, n_cell_ct, pheno_label) %>% 
  tally()   %>% 
  group_by(pheno_label) %>% 
  mutate(cell_rank = frank(n, ties.method = "first"),
         prop = n * as.numeric(ct_sig) / sum(as.numeric(ct_sig) * n_cell_ct)) %>% 
  arrange(pheno_label, -prop) %>% 
  mutate(cum_prop = cumsum(prop), lag_cum_prop = lag(cum_prop, default = 0)) %>% 
  setDT()

df_mcp_sig_count[df_stats, `:=`(assoc_mcp = i.assoc_mcp,
                                n_cell = i.n_cell,
                                assoc_mcz = i.assoc_mcz),
                 on = c("pheno_label", "cell_type")]

df_pheno_order <- df_mcp_sig_count %>% 
  group_by(pheno_label) %>% 
  summarise(n = sum(n * as.numeric(assoc_mcp < 0.05) * as.numeric(ct_sig)),
            n_prop = n / sum(as.numeric(ct_sig) * n_cell_ct)) %>% 
  arrange(desc(n_prop))

(p_sig_count <- df_mcp_sig_count %>% 
    mutate(pheno_label = factor(pheno_label, df_pheno_order$pheno_label)) %>%
    ggplot(aes(x = prop, y = pheno_label, group = pheno_label)) +
    theme_void() +
    # geom_col(fill = "#ECEFF1FF", data = ~filter(.x, assoc_mcp >= 0.05),
    #          color = "white", linewidth = 0.2) +
    geom_col(aes(fill = cell_type), width = 0.5,
             data = ~filter(.x, assoc_mcp < 0.05),
             color = "black", linewidth = 0.5) +
    geom_text(aes(label = paste0(" ", percent(x, 1)), x = x), hjust = 0,
              size = 10/.pt,
              data = ~.x %>%
                filter(ct_sig) %>% 
                group_by(pheno_label) %>%
                summarise(x = max(cum_prop))) +
    scale_fill_paletteer_d("ggthemes::Tableau_20",
                          #  breaks = df_cell_map$cell_type,
                           name = NULL,
                           guide = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 2)) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.1)),
                       labels = scales::percent,
                       # limits = c(0, max(df_mcp_sig_count$n) * 1.05),
                       guide = guide_axis(cap = "lower")) +
    coord_cartesian(clip = "off") +
    labs(x = "Average proportion of phenotype-enriched cells") +
    scale_y_discrete(labels = ~str_wrap(.x, 15)) +
    theme(legend.position = "bottom",
          panel.spacing.y = unit(0.1, "lines"),
          panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
          # axis.text.y = element_text(size = 8, hjust = 1),
          axis.text.x = element_text(size = 10, hjust = 0.5),
          axis.line.x = element_line(arrow = arrow(length = unit(0.2, "lines"), ends = "last")),
          axis.ticks.x = element_line(),
          axis.text.y = element_text(size = 10, hjust = 1),
          axis.title.x = element_text(size = 10),
          axis.ticks.length.x = unit(0.25, "lines"),
          # plot.margin = margin(1,1,1,1, unit = "lines"),
          legend.text = element_text(size = 9, margin = margin(l = 2)),
          legend.box.margin = margin(l = -5, unit = "lines"),
          legend.key.size = unit(0.5, "lines"),
          legend.key.spacing.y = unit(0.1, "lines"),
          strip.background = element_blank())
)

ggsave("figures/scdrs/bioheart975-full.scdrs.cell_type_sig_count.png", bg = "white",
       p_sig_count, width = 8, height = 5, device = agg_png, scaling = 1)

# export to supplementary tables
tab_scdrs <- df_stats %>% 
  left_join(df_mcp_sig_count) %>%
  mutate(prop_ct = n / n_cell_ct) %>% 
  mutate(pheno_label = factor(pheno_label, df_pheno_order$pheno_label)) %>% 
  arrange(pheno_label, desc(prop))

# write_gs(tab_scdrs, "scdrs_results")

# Visualize overall cell stats / count
# df_stats_pheno <- df_stats %>% 
#   group_by(pheno_label, phenotype, pheno_cat, supercategory) %>% 
#   mutate(cell_rank = dense_rank(assoc_mcz))

# number of enriched traits across cell types
df_n_traits_by_cells <- df_mcp_sig_count %>% 
  group_by(cell_type) %>% 
  mutate(pheno_length = n(),
         pheno_rank = dense_rank(n))

# trait_cat_col <- select(df_trait_cat, cat_order, color) %>% deframe()

(p_bar_celltype_stats <- df_stats %>% 
  ggplot(aes(x = pheno_label, y = -log10(assoc_mcp))) +
  geom_col(aes(fill = pheno_label, alpha = ifelse(assoc_mcp < 0.05, 1, 0.2)),
           width = 1, colour = "white") +
  scale_alpha_identity() +
  # geom_violin(aes(colour = pheno_cat, fill = after_scale(alpha(colour, 0.5))),
  #             linewidth = 0.1) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "red3") +
  labs(x = NULL, y = bquote(-log[10] ~italic("P")["phenotype enrichment"])) +
  facet_wrap(~cell_type, scales = "free_x", nrow = 4) +
  scale_fill_paletteer_d("ggthemes::Classic_10", name = NULL) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.clip = "off",
        strip.text = element_text(face = "bold", size = 8),
        legend.position = "bottom",
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
)

ggsave("figures/scdrs/bioheart975-full.scdrs.cell_type_stats_bar.png", bg = "white",
       p_bar_celltype_stats, width = 10, height = 8, device = agg_png, scaling = 1.1)

# significance count of phenotypes by cell type
df_pheno_count <- df_mcp_sig_count %>% 
  group_by(cell_type, pheno_label) %>% 
  summarise(sig_count = sum(n * as.numeric(assoc_mcp < 0.05)),
            sig_prop = sum( n / n_cell * as.numeric(assoc_mcp < 0.05))) %>% 
  ungroup() %>% 
  complete(nesting(pheno_label), nesting(cell_type),
           fill = list(sig_count = 0, sig_prop = 0))

# quantify numbers (for presentation / manuscript)
df_stats[assoc_mcp < 0.05,.N, by = cell_type][, .(min = min(N), median = median(N), max = max(N))]

# n cell with enrichment of  ≥1 phenotype
n_enriched_cells <- df_mcp[assoc_mcp < 0.05, length(unique(index))]
n_enriched_cells / df_mcp[,.N]

# enriched phenotypes

df_mcp_sig_count[prop > 0, .(sum_prop = sum(prop, na.rm = T)), by = pheno_label][,median(sum_prop)]


# Crohn's disease example 

# Visualize the cell score UMAP
# (p_cell_type_umap <- ggplot(df_cells, aes(x = umap_1, y = umap_2, color = cell_type)) +
#     theme_bw() +
#     geom_point(size = 0.1, alpha = 0.4) +
#     labs(x = "UMAP 1", y = "UMAP 2") +
#     coord_fixed() +
#     scale_color_manual(values = deframe(df_cell_map[,.(cell_type, color)]),
#                        breaks = df_cell_map$cell_type,
#                        name = NULL,
#                        guide = guide_legend(override.aes = list(size = 3, alpha = 1)))
# )

# Visualize cell score UMAP for crohns
df_cell_long <- df_cells %>% 
  pivot_longer(c(dbp, primary_ht, sbp), names_to = "phenotype", values_to = "scdrs") %>% 
  group_by(cell_type, phenotype) %>% 
  mutate(across(c(umap_1, umap_2), ~ .x %in% boxplot.stats(.x)$out,
         .names = "outlier_{col}"))

setDT(df_cell_long)
df_cell_long[df_trait_map, pheno_label := i.trait_name,
             on = c(phenotype = "trait_id")]
df_cell_long[df_mcp, assoc_mcp := i.assoc_mcp,
             on = c("index", "phenotype", "cell_type")]
df_cell_long[df_stats, ct_sig := i.assoc_mcp < 0.05,
             on = c("cell_type", "phenotype")]

# calculate centroid
centroid <- function(x) mean(x[!x %in% boxplot.stats(x)$out], na.rm = T)

df_labels <- df_cell_long[ct_sig == TRUE, .(umap_1 = centroid(umap_1),
                                         umap_2 = centroid(umap_2)),
                          by = list(cell_type, pheno_label)]

set.seed(100)
(p_scdrs_umap <- df_cell_long %>% 
    group_by(phenotype) %>% 
    # randomise cell order so that overlap is random
    slice_sample(n = nrow(.)) %>%
    # arrange(-assoc_mcp) %>%
    ggplot(aes(x = umap_1, y = umap_2)) +
    theme_void() +
    facet_wrap(~pheno_label) +
    geom_point(aes(colour = ifelse(ct_sig, scdrs, NA)),
               # data = ~filter(.x, assoc_mcp < 0.05 & ct_sig),
               shape = "circle",
               # show.legend = TRUE,
               size = 0.2, alpha = 0.5) +
    geom_label_repel(
      aes(label = cell_type), color = "black", 
      fontface = "bold", force = 1, label.size = NA, label.padding = 0.1,
      size = 11/.pt, fill = alpha("white", 0.8),
      data = df_labels, max.overlaps = Inf
    ) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    coord_fixed() +
    scale_colour_viridis_c(
      name = "scDRS", option = "magma", direction = -1,
      na.value = "gray90") +
    theme(panel.border = element_rect(linewidth = 1, color = "black"),
          legend.position = "right",
          legend.position.inside = c(0.75, 0.02),
          strip.text = element_text(face = "bold", size = 12),
          strip.clip = "off",
          aspect.ratio = 1,
          # legend.justification = c("left", "bottom"),
          axis.title.x = element_text(size = 12, margin = margin(t = 0.5, unit = "lines")),
          axis.title.y = element_text(size = 12, margin = margin(r= 0.5, unit = "lines"),
                                      angle = 90))
)

# umap celltype legend
(p_umap_legend <- df_cells %>%
    ggplot(aes(x = umap_1, y = umap_2, color = cell_type)) +
    theme_void() +
    geom_point(size = 0.1, alpha = 0.4) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    coord_fixed() +
    scale_color_manual(values = deframe(df_cell_map[,.(cell_type, color)]),
                       breaks = df_cell_map$cell_type,
                       name = NULL,
                       guide = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
    theme(legend.position = "right",
          legend.text = element_text(size = 9, margin = margin(l = 2)),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10, angle = 90),
          legend.key.size = unit(0.5, "lines"),
          aspect.ratio = 1,
          panel.border = element_rect(linewidth = 1, color = "black"),
          legend.key.spacing.y = unit(0.2, "lines"))
)

ggsave("figures/scdrs/bioheart975-full.scdrs.crohns.umap.png", bg = "white",
       p_scdrs_umap, width = 12, height = 5, scaling = 1, device = agg_png)

ggsave("figures/scdrs/bioheart975.scdrs.cell_type_umap_legend.png", bg = "white",
       p_umap_legend, width = 8, height = 8, device = agg_png, scaling = 1.5)

# Ridge plot
df_cell_order <- df_stats %>%
  group_by(phenotype) %>% 
  arrange(-assoc_mcz) %>%
  mutate(cell_type = factor(cell_type, unique(cell_type)),
         sig = assoc_mcp < 0.05,
         text_col = ifelse(sig, "red3", "black"),
         text_face = ifelse(sig, "bold", "plain")) %>% 
  setDT()

df_avg <- df_cell_long[,.(mean_scdrs = mean(scdrs, na.rm = TRUE),
                       sd_scdrs = sd(scdrs, na.rm = TRUE),
                       n_cells = .N,
                       ct_sig = !is.na(max(ct_sig)),
                       n_sig = sum(as.numeric(assoc_mcp < 0.05)),
                       prop_sig = sum(as.numeric(assoc_mcp < 0.05)) / .N),
                    by = list(cell_type, pheno_label)]

plusmin <- "±"
(p_ridge <- df_cell_long %>%
    ggplot(aes(x = scdrs, y = fct_rev(cell_type))) +
    theme_bw() +
    geom_density_ridges(
      aes(colour = ct_sig, fill = after_scale(alpha(colour, 0.5))),
      scale = 1.4, linewidth = 0.2, rel_min_height = 0,
    ) +
    facet_grid(cols = vars(pheno_label), scales = "free_y", space = "free") +
    labs(x = "scDRS Crohn's disease", y = NULL) +
    scale_colour_manual(
      values = c("FALSE" = "gray50", "TRUE" = "red3"),
      name = NULL,
      guide = "none") +
    # scale_colour_manual(
    #   values = deframe(df_cell_map[,.(cell_type, color)]),
    #   name = NULL,
    #   guide = "none") +
    geom_text(aes(label = paste(number(mean_scdrs, 0.1), plusmin, number(sd_scdrs,0.1))),
              data = df_avg[ct_sig == TRUE], size = 10/.pt, color = "red3",
              x = Inf, vjust = 0, position = position_nudge(y = 0.1), hjust = 1) +
    geom_text(aes(label = paste(number(mean_scdrs, 0.1), plusmin, number(sd_scdrs,0.1))),
              data = df_avg[ct_sig == FALSE], size = 10/.pt, color = "black",
              x = Inf, vjust  = 0, position = position_nudge(y = 0.1), hjust = 1) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             fontface = "plain", label = paste("Mean", plusmin, "SD"), size = 10/.pt) +
    scale_x_continuous(labels = label_number(accuracy = 0.01),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    theme(axis.title.x = element_text(size = 10),
          axis.text.y = element_text(vjust = 0, size = 10,
                                     face = rev(df_cell_order$text_face),
                                     colour = rev(df_cell_order$text_col)),
          axis.text.x = element_text(size = 10),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0.3, "lines"),
          axis.line.x = element_line(),
          panel.border = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing.x = unit(2, "lines"),
          # plot.margin = margin(l = -5, unit = "lines"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major = element_blank())
)

ggsave("figures/scdrs/bioheart975-full.scdrs.cell_type_distribution.png", bg = "white",
       p_ridge, width = 9, height = 9, device = agg_png, scaling = 0.9)



(p_bar_prop <- df_avg %>%
    ggplot(aes(y = factor(cell_type, rev(df_cell_order$cell_type)))) +
    geom_col(aes(x = n_cells), position = position_nudge(y = 0.25),
             fill = "gray80", width = 0.5,
             color = "black", linewidth = 0.5) +
    geom_col(aes(x = n_sig, fill = "red3"), position = position_nudge(y = 0.25),
             # fill = "#4E4E4EFF",
             colour = "black", width = 0.5, linewidth = 0.5) +
    theme_bw() +
    scale_fill_identity(
      name = NULL, labels = expression(italic("P")["scDRS"] < 0.05),
      guide = guide_legend()
    ) +
    geom_text(aes(x = Inf, color = ct_sig, label = paste0(" ", percent(prop_sig, 0.1))),
              size = 10/.pt, position = position_nudge(y = 0.1),
              hjust = 1, vjust = 0) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red3"),
                       guide =  "none") +
    labs(x = "Number of cells", y = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
             fontface = "plain", label = bquote("% enriched cells"), size = 10/.pt) +
    scale_x_continuous(expand = expansion(mult = c(0.01)),
                       limits = c(0, 1.7e4),
                       breaks = c(0, 5000, 10000),
                       labels = c("0", "5k", "10k"),
                       guide = guide_axis(cap = "both")) +
    layer_scales(p_ridge)$y +
    coord_cartesian(clip = "off") +
    theme(axis.title.x = element_text(size = 11, hjust = 0),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.ticks.y = element_blank(),
          axis.ticks.length.x = unit(0.3, "lines"),
          axis.line.x = element_line(),
          panel.border = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(l = 0.5, unit = "lines"),
          legend.position = "none",
          panel.grid.major.y = element_blank())
)


(p_ridge_bar <- (p_ridge + p_bar_prop)  +
    plot_layout(width = c(1, 0.4)) 
)


# Decile plot
df_q_scdrs <- df_cells %>%
  select(crohns, B2_ALL_eur_leave_23andme, cell_label) %>% 
  pivot_longer(-cell_label, names_to = "phenotype", values_to = "scdrs") %>%
  left_join(df_trait_map %>% select(phenotype = trait_id, pheno_label = name)) %>%
  group_by(pheno_label, phenotype, cell_label) %>%
  mutate(q_scdrs = ntile(scdrs, 10)) %>% 
  group_by(pheno_label, phenotype, cell_label, q_scdrs) %>%
  summarise(mean = mean(scdrs), sd = sd(scdrs), n = n()) %>% 
  setDT()

df_q_scdrs[df_stats, ct_sig := i.assoc_mcp < 0.05, on = c("cell_label", "phenotype")]

dendritic_cells <- c("ASDC", "cDC1", "cDC2", "pDC")

(p_decile <- df_q_scdrs %>% 
    filter(cell_label %in% dendritic_cells) %>% 
    ggplot(aes(x = q_scdrs, y = mean,
               alpha = ifelse(ct_sig, 1, 0.4),
               colour = pheno_label, group = pheno_label)) +
    facet_wrap(~cell_label, axes = "all") +
    # geom_line(linewidth = 0.5, colour = "#ECEFF1FF",
    #           data = ~filter(.x, !(cell_label %in% df_cell_order[sig == TRUE, cell_label]))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 2) +
    geom_line(aes(label = pheno_label), linewidth = 1) +
    scale_alpha_identity() +
    # geom_text_repel(
    #   aes(label = cell_label), size = 10/.pt,
    #   data = ~filter(.x, q_scdrs == 10 &
    #                    cell_label %in% df_cell_order[sig == TRUE, cell_label]),
    #   direction = "y",
    #   colour = "black",
    #   xlim = c(10.2, Inf),
    #   min.segment.length = 0,
    #   segment.size = 0.1, show.legend = FALSE) +
    labs(y = "Mean scDRS",
         x = "scDRS decile") +
    coord_cartesian(clip = "off") +
    scale_x_continuous(breaks = 1:10, labels = 1:10,
                       limits = c(1, 10),
                       expand = expansion(add = c(0.4)),
                       guide = guide_axis(cap = "both")) +
    scale_colour_paletteer_d(
      "ggthemes::Tableau_10",
      # labels = ~str_wrap(.x, 15),
      name = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 12),
          legend.position = "bottom",
          legend.box.spacing = unit(0.2, "lines"),
          plot.margin = margin(l = 1, unit = "lines"),
          # legend.position.inside = c(0.99, 0.01),
          # legend.justification = c("right", "bottom"),
          axis.line = element_line(),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          axis.ticks.length = unit(0.3, "lines"))
)

# combine plots

design <- "
AAAAABBBBB
AAAAABBBBB
AAAAACCCCC
AAAAACCCCC
AAAAACCCCC
AAAAADDDDD
AAAAADDDDD
"
list_p <- list(A = p_sig_count,
               B = p_scdrs_umap_crohns,
               C = p_ridge_bar,
               D = p_decile) %>% 
  map(wrap_elements)

(plots <- wrap_plots(list_p, design = design) +
    plot_annotation(tag_levels = 'a') &
    theme(plot.tag = element_text(size = 16, face = "bold"),
          plot.tag.position = c(0.01, 1))
)

ggsave("figures/scdrs/tenk10k_phase1_sample.scdrs.crohns.all.v6.png",
       plots, bg = "white", width = 12.5, height = 18, device = agg_png, scaling = 1.2)

ggsave("figures/publication_pdf/Fig3 - Polygenic Enrichment.pdf",
       plots, bg = "white", width = 12.5, height = 18, device = cairo_pdf, scale = 1/1.2)

# 
# ggsave("figures/scdrs/tenk10k_phase1_sample.scdrs.all.png",
#        p_cell_type_umap, width = 8, height = 6, device = agg_png)
# 
# ggsave("figures/scdrs/tenk10k_phase1_sample.scdrs.crohns.umap.png", bg = "white",
#        p_scdrs_score_crohns, width = 9, height = 6, scaling = 1.4, device = agg_png)
# 
# ggsave("figures/scdrs/tenk10k_phase1_sample.scdrs.crohns.violin.png", bg = "white",
#        p_violin, width = 10, height = 6, device = agg_png, scaling = 1.8)
# 
# ggsave("figures/scdrs/tenk10k_phase1_sample.scdrs.crohns.ridge.png", bg = "white",
#        p_ridge, width = 9, height = 6, device = agg_png, scaling = 1.8)
# 
# ggsave("figures/scdrs/tenk10k_phase1_sample.scdrs.crohns.decile.png", bg = "white",
#        p_decile, width = 8, height = 6, device = agg_png, scaling = 1.8)

# Scratch
df_pheno_count %>% 
  ggplot(aes(x = pheno_label, y = cell_label)) +
  geom_point(aes(fill = sig_count),
             shape = "circle filled",
             position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~cell_label, scales = "free", nrow = 7) +
  scale_fill_paletteer_d("ggthemes::Classic_10", name = "Phenotype category") +
  labs(y = "Counts of Enriched Phenotypes") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank())

phenotypes <- intersect(df_trait_map[,trait_id], colnames(df_cells))
df_cell_q99 <- df_cells %>% 
  pivot_longer(all_of(phenotypes), names_to = "phenotype") %>% 
  group_by(phenotype) %>% 
  filter(value > quantile(value, 0.99)) %>% 
  group_by(phenotype, cell_label) %>% 
  tally() %>%
  mutate(cell_rank = frank(n, ties.method = "first"),
         prop = n / sum(n)) %>% 
  arrange(phenotype, -n) %>% 
  mutate(cum_prop = cumsum(prop), lag_cum_prop = lag(cum_prop, default = 0)) %>% 
  setDT()

df_cell_q99[df_trait_map, `:=`(pheno_label = i.name, pheno_cat = i.cat_rev,
                               pheno_slim = i.include,
                               supercategory = i.supercategory),
            on = c(phenotype = "trait_id")]

df_cell_q99[df_stats, `:=`(assoc_mcp = i.assoc_mcp),
            on = c("phenotype", "cell_label")]

(p_cell_q99 <- df_cell_q99 %>% 
    group_by(pheno_cat) %>% 
    mutate(pheno_num = dense_rank(pheno_label)) %>% 
    filter(phenotype != "pancreas") %>% 
    ggplot() +
    theme_void() +
    geom_rect(aes(xmin = lag_cum_prop, xmax = cum_prop,
                  ymin = pheno_num - 0.5, ymax = pheno_num + 0.5),
              data = ~filter(.x, prop <= 0.1),
              fill = "#ECEFF1FF", colour = "gray90", linewidth = 0.5) +
    geom_rect(aes(fill = cell_label, xmin = lag_cum_prop, xmax = cum_prop,
                  ymin = pheno_num - 0.5, ymax = pheno_num + 0.5),
              data = ~filter(.x, prop > 0.1),
              colour = "black", linewidth = 0.5) +
    facet_wrap(vars(pheno_cat), dir = "h", ncol = 1,
               space = "free_y", scales = "free") +
    # geom_col(aes(fill = as.factor(cell_rank)), width = 0.5) +
    scale_fill_manual(values = deframe(df_cell_map[,.(cell_type, color)]),
                      breaks = df_cell_map$cell_type,
                      name = NULL,
                      guide = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
    scale_x_continuous(expand = expansion(add = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(0.01)) +
    geom_text(
      aes(label = pheno_label, y = pheno_num), size = 8/.pt,
      x = -Inf, hjust = 1, nudge_x = 0.01,
      data = ~distinct(.x, pheno_num, pheno_label, pheno_cat)) + 
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Proportion of cells in 99th percentile") +
    theme(panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          strip.clip = "off",
          strip.text = element_text(hjust = 0, face = "bold",
                                    margin = margin(l = 0.75, unit = "lines")),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(size = 10),
          legend.position = "right",
          plot.margin = margin(l = 14, r = 1, unit = "lines"),
          panel.grid.minor.x = element_blank())
)

ggsave("figures/scdrs/tenk10k_phase1_sample.scdrs.cell_type_q99.png", bg = "white",
       p_cell_q99, width = 8, height = 16, device = agg_png, scaling = 1)



# visualize mean scDRS per phenotype - cell type
df_mean_sum <- df_cells %>% 
  pivot_longer(all_of(phenotypes), names_to = "phenotype") %>%
  group_by(cell_type, cell_label, phenotype) %>% 
  summarise(mean_scdrs = mean(value, na.rm = TRUE),
            sd_scdrs = sd(value, na.rm = TRUE), n_cells = n())

setDT(df_mean_sum)
df_mean_sum[df_trait_map, `:=`(pheno_label = i.name, pheno_cat = i.cat_rev,
                               pheno_slim = i.include,
                               supercategory = i.supercategory),
            on = c(phenotype = "trait_id")]
df_mean_sum[df_stats, `:=`(assoc_mcp = i.assoc_mcp, assoc_mcz = i.assoc_mcz),
            on = c("phenotype", "cell_label")]

df_mean_sum %>% 
  ggplot(aes(y = assoc_mcz, x = pheno_label)) +
  theme_bw() +
  geom_point(aes(colour = cell_label), data = ~.x[assoc_mcp < 0.05]) +
  facet_wrap(~pheno_cat, space = "free_x", scale = "free_x") +
  scale_colour_manual(values = deframe(df_cell_map[,.(cell_type, color)]),
                      breaks = df_cell_map$cell_type,
                      name = NULL,
                      guide = guide_legend(override.aes = list(size = 3, alpha = 1)))


df_stats_pheno %>% 
  filter(assoc_mcp < 0.05) %>%
  filter(supercategory == "disease", pheno_slim) %>% 
  ggplot(aes(y = assoc_mcz, x = cell_rank)) +
  facet_wrap(~pheno_label) +
  theme_void() +
  # geom_point(color = "#ECEFF1FF", data = ~filter(.x, assoc_mcp >= 0.05)) +
  geom_col(aes(fill = cell_label),
           data = ~filter(.x, assoc_mcp < 0.05)) +
  # coord_polar() +
  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray20") +
  scale_fill_manual(values = deframe(df_cell_map[,.(cell_type, color)]),
                    breaks = df_cell_map$cell_type,
                    name = NULL,
                    guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(strip.background =element_blank(),
        strip.clip = "off",
        panel.border = element_rect(),
        # legend.position = "bottom",
        panel.grid.major.y = element_line(colour = "gray95"),
        axis.text = element_blank())
