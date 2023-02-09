library(tidyverse)
library(rstatix)
library(ComplexHeatmap)
library(circlize)

# Supp Figure 1: qual metabolites in diversity low, medium, high groups

# load in data and subset to patient level ----------------------------

heatmap_cmpds <- read_csv("./data/qual_compounds.csv") 
heatmap_lookup <- read_csv("./data/qual_heatmap_lookup.csv")

metab <- read_csv("./data/LD850.meta.quant.metabolomics.csv") %>% 
  filter(!is.na(firstSample))

# change compound names to proper format ----------------------------------

qual <- read_csv("./data/LD850.qual.metabolomics.csv") %>% 
  filter(! (compound %in% c("indole-3-carboxaldehyde",
                            "indole-3-acetate",
                            "indole-3-propionate",
                            "trans-indole-3-acrylate"))) %>% 
  filter(seq_id %in% metab$seq_id) %>% 
  mutate(compound = case_when(
    compound == "indole3carboxyaldehyde" ~ "indole3carboxaldehyde", 
    compound == "isovaleric-acid" ~ "isovalerate",
    compound == "indole3acetic" ~ "indole3acetate", 
    compound == "indole3acrylicacid" ~ "indole3acrylate", 
    compound == "indole3propionic" ~ "indole3propionate", 
    compound == "indole3lacticacid" ~ "indole3lactate",
    TRUE ~ compound),
    compound = recode(compound,
                      Preq1 = "PreQ1"))

qual.temp <- qual %>% 
  mutate(compound = str_to_title(compound)) %>% 
  filter(compound %in% heatmap_cmpds$compound|is.na(compound)) %>% 
  count(compound) 

namemap <- heatmap_cmpds %>% 
  anti_join(qual.temp) %>% 
  mutate(oldname = tolower(gsub("\\-| ","", compound)))

rename.var <- namemap$compound
names(rename.var) <- namemap$oldname

qualdata <- qual %>% 
  mutate(compound = if_else(compound %in% names(rename.var),
                            rename.var[compound], compound),
         compound = str_to_title(compound),
  ) %>% 
  filter(compound %in% heatmap_cmpds$compound|is.na(compound))  

# select samples with complete cases of all 83 compounds ------------------

comp_val <- qualdata %>% 
  group_by(compound) %>% 
  summarise(median_val = median(value),
            min_val = min(value[value > 0])) %>% 
  mutate(median_val = if_else(median_val == 0, min_val/10, median_val)) 

heatmap_data <- qualdata %>% 
  left_join(comp_val) %>% 
  mutate(heatmap_val = ifelse(log(value/ median_val, base = 2) == -Inf, 
                              0, log(value/ median_val, base = 2))) %>% 
  ungroup() %>% 
  select(-c(value, median_val, min_val)) %>% 
  group_by(seq_id, compound) %>% 
  slice_max(heatmap_val, with_ties = F, n = 1) %>% 
  left_join(metab %>% 
              select(seq_id, invSimp)) %>% 
  left_join(heatmap_lookup) %>% 
  mutate(class = factor(
    class,
    levels = c(
      "Fatty Acid",                   # 1
      "Dicarboxylic Acid",            # 2
      "Amino Acid",                   # 3
      "Bile Acid",                    # 4
      "Indole",                       # 5
      "Phenolic Aromatic",            # 6
      "Kynurine Pathway",             # 7
      "Vitamin"                       # 8
    )
  ),
  subclass = factor(
    subclass,
    levels = c(
      "Short-Chain Fatty Acid",       # 1
      "Branched-Chain Fatty Acid",    # 2
      "Aminated Fatty Acid",          # 3
      "Long-Chain Fatty Acid",        # 4
      "Dicarboxylic Acid",            # 5
      "Amino Acid",                   # 6
      "Primary Bile Acid",            # 7
      "Secondary Bile Acid",          # 8
      "Conjugated Bile Acid",         # 9
      "Indole",                       # 10
      "Phenolic Aromatic",            # 11
      "Kynurine Pathway",             # 12
      "Vitamin"                       # 13
    )
  )) %>% 
  arrange(class, subclass, compound) %>% 
  ungroup() %>% 
  add_count(seq_id, name = "n_metab") %>% 
  filter(n_metab >= nrow(comp_val))

# 237 samples to plot
# heatmap_data %>% 
#   count(ID)

# creat heatmap order for rows and columns --------------------------------

heatmap_column_df <- heatmap_data %>% 
  distinct(seq_id, invSimp) %>% 
  arrange(invSimp) 

cuts <- quantile(heatmap_column_df$invSimp, probs = c(0.33, 0.67))
invSimp_max <- ceiling(max(heatmap_column_df$invSimp))

heatmap_sample_order <- heatmap_column_df %>% 
  mutate(diversity = cut(invSimp, breaks = c(0, cuts, invSimp_max), 
                         labels = c("Low", "Medium","High")))

# ks test for low, med and high inverse simpson groups --------------------

metab_p <- heatmap_data %>% 
  left_join(heatmap_sample_order) %>% 
  group_by(compound) %>%
  rstatix::kruskal_test(heatmap_val~diversity) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  select(compound, statistic, p, p.adj) 

heatmap_metab_order <- heatmap_data %>%
  select(class, subclass, compound) %>% 
  distinct() %>% 
  left_join(metab_p) %>% 
  arrange(class, subclass, p.adj)

# check range of log2 folder change
range(heatmap_data$heatmap_val)

# Heatmap legend color
col_fun <- colorRamp2(breaks = c(-6, 0, 6), colors = c("#00aaad", "white", "#ad003a"))

# Global parameter for annotation 
# ht_opt$COLUMN_ANNO_PADDING <- unit(2.5, "mm")
heatmap_mat <-
  heatmap_data %>% 
  pivot_wider(c(seq_id), names_from = compound, 
              values_from = heatmap_val) %>%
  column_to_rownames(var = "seq_id") %>% 
  as.matrix() %>% 
  t() 

heatmap_mat <- heatmap_mat[heatmap_metab_order$compound,heatmap_sample_order$seq_id]

inv_df <- heatmap_data %>% 
  distinct(seq_id, invSimp) %>% 
  arrange(invSimp) %>% 
  column_to_rownames(var = "seq_id") 

top_anno = HeatmapAnnotation(`Inverse Simpson` = anno_barplot(inv_df$invSimp))

pvalue_col_fun = colorRamp2(c(0, 0.045, 0.2), c("#06d106", "#ebe121","#E3E4E6"))
row_anno = rowAnnotation(`Adjusted KS\nP-values` = heatmap_metab_order$p.adj,
                         col = list(`Adjusted KS\nP-values` = pvalue_col_fun))

draw(row_anno)

gg_metab_heatmap <- Heatmap(heatmap_mat, 
                            name = "Fold Change (log2)",
                            col = col_fun,
                            # na_col = "grey83",
                            na_col = "black",
                            rect_gp = gpar(col = "grey40", lwd = 1.5),
                            column_names_gp = grid::gpar(fontsize = 6.5),
                            column_gap = unit(2.5, "mm"),
                            column_title_gp = gpar(fontsize = 14),
                            column_title_rot = 0,
                            column_split = heatmap_sample_order$diversity,
                            cluster_columns = FALSE,
                            show_column_names = FALSE,
                            show_column_dend = FALSE,
                            row_names_gp = gpar(fontsize = 8),
                            row_gap = unit(2.25, "mm"),
                            row_names_side = c("left"),
                            row_split = heatmap_metab_order$subclass,
                            row_title_rot = 0,
                            cluster_rows = FALSE,
                            show_row_dend = FALSE,
                            top_annotation = top_anno,
                            right_annotation = row_anno,
                            heatmap_height = unit(14, "in"),
                            heatmap_width =  unit(23.5, "in")
)

gg_metab_heatmap

## print out heatmap -------------------------------------------------------

pdf("./results/7.LD237.qual_heatmap.firstsample.pdf", height = 15, width = 28, onefile = F)
gg_metab_heatmap
dev.off()
