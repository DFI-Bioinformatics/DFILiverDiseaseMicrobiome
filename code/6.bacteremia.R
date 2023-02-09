library(tidyverse)
library(EnhancedVolcano)
library(rstatix)
library(microbiomeMarker)
library(phyloseq)

# Figure 6 + Supp Figure 9: analysis on bacteremia vs. No bacteremia 

# load source code for plot -----------------------------------------------

source("./code/getRdpPal.R")
source("./code/utility.source.R")

# load in data ------------------------------------------------------------

metab <- read_csv("./data/LD850.meta.quant.metabolomics.csv") 

mpa <- read_csv("./data/LD850.metaphlan.csv")

qual <- read_csv("./data/LD850.qual.metabolomics.csv") 

ratio <- read_csv("./data/LD850.ratios.quant.metabolomics.csv")

# subset to bacteremia cohort ----------------------------------------------------

meta <- metab %>% 
  filter(!is.na(bacteremia))

mpa <- mpa %>% 
  filter(seq_id %in% meta$seq_id)

qual <- qual %>% 
  filter(seq_id %in% meta$seq_id) %>% 
  left_join(meta %>% 
              select(seq_id, ID, bacteremia))

ratio <- ratio %>% 
  filter(seq_id %in% meta$seq_id) %>% 
  gather("ratios", "value", -seq_id, -ID) %>% 
  left_join(meta %>% 
              select(seq_id, ID, bacteremia))

# tax + sbp sorted by bifido ----------------------------------------------

## taxonomy barplot --------------------------------------------------------

t <- mpa %>%
  add_count(seq_id, wt = relative_abundance, name = "totalAbd") %>%
  mutate(pctseqs = relative_abundance/totalAbd)  %>%
  dplyr::count(seq_id, taxid, Kingdom, Phylum, Class, Order, Family,
               Genus, Species, wt = pctseqs, name = "pctseqs") %>%
  add_count(seq_id, wt = pctseqs, name = "Total") %>%
  mutate(pctseqs = pctseqs/Total)  %>%
  mutate(genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus)) %>%
  left_join(meta %>% select(seq_id, ID)) 

## color palette
pal <- getRdpPal(t)

taxgg <- t %>%
  left_join(meta %>% distinct(ID, bifido_percent = Bifidobacterium, bacteremia)) %>%
  ggplot(aes(x=reorder(ID,-bifido_percent), y=pctseqs)) +
  geom_bar(aes(fill=Genus),stat="identity") +
  # geom_text(aes(y=1-y.text,x=ID,label=tax.label),angle=90,
  #           lineheight=0.6,size=2.2) +
  scale_fill_manual(values=pal) +
  theme_bw() +
  # facet_grid(.~ sbp_diag, scales = "free", space = "free") +
  theme(legend.position = "none",
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_blank(),
        #axis.text.x = element_blank()
        strip.background.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        axis.ticks.x = element_blank()
  ) +
  labs(x="", y= "Bacterial Abund.") +
  geom_hline(aes(yintercept = 0.9), linetype = 2, color = "white") +
  scale_y_continuous(breaks = c(seq(0, 1, 0.25), 0.9))

taxgg

## bacteremia yes vs. no ---------------------------------------------------------

sbpgg <- meta %>%
  ggplot(aes(reorder(ID, -Bifidobacterium), "Bacteremia status")) +
  geom_tile(aes(fill = bacteremia)) +
  labs(x="", y= "") +
  scale_fill_manual(values = c("no bacteremia" = "#09ba38", "bacteremia" = "#db3b3b")) +
  theme_bw() +
  theme(# plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.text.x = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, face = "bold"),
    #axis.text.x = element_blank()
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.02, "lines"),
    axis.ticks.x = element_blank()
  )

## combine all plots -------------------------------------------------------

pdf("./results/6.LD228.Bacteremia.metaphlan.pdf",
    height = 5.25, width = 12.5)
gg.stack(taxgg, sbpgg, 
         heights = c(1,0.13),
         newpage = F)
dev.off()

# lefse analysis between Bacteremia vs. no Bacteremia ----------------------

## transform into matrix ----------

kmat <- mpa %>%
  transmute(seq_id, taxid, 
            numseqs = as.numeric(estimated_number_of_reads_from_the_clade)) %>%
  pivot_wider(names_from =taxid,
              values_from = numseqs,
              values_fn = sum,
              values_fill = 0) %>%
  column_to_rownames(var = "seq_id")

# check dimension
dim(kmat)

## taxonomy table needed for phyloseq object --------
tax_df <- mpa %>%
  distinct(taxid, Kingdom, Phylum, Class, Order, Family, Genus, Species)  %>%
  relocate(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  column_to_rownames(var = "taxid") %>%
  as.matrix()


## make into a phyloseq object ---------------------------------------------

krakenphy <- phyloseq(tax_table(tax_df),
                      otu_table(kmat, taxa_are_rows = F),
                      sample_data(meta %>%
                                    distinct(ID, seq_id, bacteremia) %>%
                                    column_to_rownames(var = "seq_id")
                      ))


## run lefse ---------------------------------------------------------------

{
  phy_lefse <- run_lefse(
    ps = krakenphy,
    group = "bacteremia",
    subgroup = NULL,
    taxa_rank = "all",
    transform = "identity",
    norm = "CPM",
    lda_cutoff = 2.0,
    sample_min = 3,
  )
  
  marker_table(phy_lefse) %>%
    as_tibble() %>% # view
    mutate(plabel = sapply(feature, last2tax),
           ef_lda = if_else(enrich_group == "bacteremia", ef_lda, -ef_lda),
           level = factor(str_extract(plabel, "^[a-z]"),
                          levels = c("k","p","c","o","f","g"),
                          labels = c("Phylum","Class","Order",
                                     "Family","Genus","Species"))
    ) %>%
    filter(!is.na(level)) %>%
    ggplot(aes(x = reorder(plabel, -ef_lda), y = ef_lda, fill = enrich_group)) +
    geom_col(color = "black", alpha = 0.8) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linetype="dashed", color = "black"),
      axis.text = element_text(size = 8.5, color = "black"),
      axis.title = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 10, color = "black"),
      legend.text = element_text(size = 9, color = "black"),
      plot.title = element_text(size = 12, color = "black"),
      strip.text.y = element_text(angle = 0)
    ) +
    paletteer::scale_fill_paletteer_d("ggthemes::Tableau_10") +
    scale_y_continuous(breaks = seq(-4,4,1))+
    coord_flip() +
    labs(y = "LDA Score (log10)", x = "Phylogeny\n", fill = "bacteremia\nStatus",
         title = "LEfSe Analysis: Grouped by bacteremia Status") +
    facet_grid("level", scales = "free", space = "free_y")
  
  ggsave(paste0("./results/6.LD228.bacteremia.",lubridate::today(),".lefse_analysis.pdf"),
         height = 6.5, width = 8.2)
  
}

# quant metabolites comparison between sbp vs. no sbp --------------------------

metab.long <- meta %>%
  select(ID, `3-oxolithocholic acid`:succinate) %>%
  gather("compound","value",-ID) %>%
  right_join(meta %>%
               distinct(ID, bifido_percent = Bifidobacterium,
                        bacteremia, bifido_class)) 

metab.long %>%
  count(ID, sort = T)

## boxplot -----------------------------------------------------------------

subtext <- meta %>%
  count(bacteremia) %>%
  mutate(lab = paste0(bacteremia,"=", n)) %>%
  pull(lab) %>%
  paste0(collapse = "; ")

metab.long %>%
  mutate(compound = factor(compound,
                           levels = c("acetate", "butyrate", "propionate", "succinate",
                                      "cholic acid","glycocholic acid", "taurocholic acid",
                                      "lithocholic acid", "deoxycholic acid", "isodeoxycholic acid",
                                      "3-oxolithocholic acid", "alloisolithocholic acid")),
         modified_diagnosis = factor(bacteremia, 
                                     levels = c("no bacteremia", "bacteremia"))) %>%
  ggplot(aes(modified_diagnosis, value, color = modified_diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  ggpubr::stat_compare_means() +
  facet_wrap("compound", scales = "free") +
  scale_color_manual(values =  paletteer::paletteer_d("ggthemes::Tableau_10")[c(2,1)] ) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(face = "bold", size = 13.5),
        strip.background.x = element_blank(),
        axis.text = element_text(size = 11.5),
        axis.title = element_text(size = 12.5)) +
  labs(x = "Bacteremia Status", y = "Metabolites Conc.\n(log transformed)",
       caption = subtext) +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.1))

ggsave(file.path("./results/6.LD228.bacteremia.metabolites.quant.pdf"), 
       width = 12.5, height = 9.3)

# quant bile acid ratios --------------------------------------------------

ratio %>% 
  mutate(modified_diagnosis = factor(bacteremia, 
                                     levels = c("no bacteremia", "bacteremia"))) %>% 
  ggplot(aes(modified_diagnosis, value, fill = modified_diagnosis)) +
  geom_boxplot(outlier.shape = NA,  alpha = 0.75) +
  geom_jitter(shape = 21, alpha = 0.65) +
  scale_fill_manual(values =  paletteer::paletteer_d("ggthemes::Tableau_10")[c(2,1)] ) +
  facet_wrap(. ~ ratios, scales = "free",
             ncol = 3) +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.001),
                     labels = scales::scientific) +
  theme_bw() +
  ggpubr::stat_compare_means() +
  theme(legend.position = "none",
        strip.text.x = element_text(face = "bold", size = 13.5),
        strip.background.x = element_blank(),
        axis.text = element_text(size = 11.5),
        axis.title = element_text(size = 12.5)) +
  labs(x = "", y = "Bile Acid ratio\n(log transformed)",
       caption = subtext, 
       title = paste0("Bile Acid ratio", " comparison ",
                      "between bacteremia vs. non-bacteremia"))

ggsave("./results/6.LD228.bacteremia_BAratios.boxplot.pdf", 
       width = 11.5, height = 9.5) 

# qual volcano plot -------------------------------------------------------

### no bacteremia as the reference group
bacteremia_log2fc <- qual %>%
  group_by(compound, bacteremia) %>%
  summarise(med = median(value)) %>%
  mutate(med = ifelse(med == 0, 0.0001, med)) %>% 
  spread(bacteremia, med) %>%
  mutate(log2fc = log2(bacteremia/`no bacteremia`))

bacteremia_pval <- qual %>%
  group_by(compound) %>%
  rstatix::wilcox_test(value ~ bacteremia) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

bacteremia_tot <- left_join(bacteremia_log2fc, bacteremia_pval) %>%
  column_to_rownames(var = "compound")

# Volcano Plot
set.seed(123456)

bacteremia_tot %>%
  filter(abs(log2fc) > 1)

xylims = ceiling(max(bacteremia_tot$log2fc))
plims <- ceiling(max(-log10(bacteremia_tot$p)))

volcano <-
  EnhancedVolcano(bacteremia_tot,
                  lab = rownames(bacteremia_tot),
                  title = 'Bacteremia versus No Bacteremia',
                  y = "p",
                  x = "log2fc",
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 3.0,
                  labSize = 2.0,
                  xlim = c(-xylims,xylims),
                  ylim = c(0, plims + 1),
                  col=paletteer::paletteer_d("ggthemes::fivethirtyeight")[c(3,2,5,4)],
                  colAlpha = 1,
                  legendPosition = "right",
                  legendLabels = c(expression(p > 0.05*";" ~ Log[2] ~ FC < "\u00B1"*1),
                                   expression(p > 0.05*";" ~ Log[2] ~ FC >= "\u00B1"*1),
                                   expression(p <= 0.05*";" ~ Log[2] ~ FC < "\u00B1"*1),
                                   expression(p <= 0.05*";" ~ Log[2] ~ FC >= "\u00B1"*1)),
                  legendLabSize = 14,
                  drawConnectors = T,
                  widthConnectors = 0.5,
                  arrowheads = F,
                  gridlines.minor = F,
                  gridlines.major = F) +
  labs(subtitle = NULL,
       y = expression( -Log[10] ~ P)) +
  annotate("text", x = 0.65*xylims, y = plims + 0.5, label = "Bacteremia",
           size = 6, color = paletteer::paletteer_d("ggthemes::Tableau_10")[1]) +
  annotate("rect", xmin = 1.05, xmax = Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = paletteer::paletteer_d("ggthemes::Tableau_10")[1]) +
  annotate("text", x = -0.65*xylims, y = plims + 0.5, label = "No Bacteremia",
           size = 6, color = paletteer::paletteer_d("ggthemes::Tableau_10")[2]) +
  annotate("rect", xmin = -1.05, xmax = -Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = paletteer::paletteer_d("ggthemes::Tableau_10")[2]) +
  guides(color = guide_legend(nrow = 4),
         shape = guide_legend(nrow = 4))

volcano

ggsave(plot = volcano, file.path("./results/6.LD228.volcano_metab_bacteremia.qual.pdf"),
       width = 12, height = 6, units = "in")
