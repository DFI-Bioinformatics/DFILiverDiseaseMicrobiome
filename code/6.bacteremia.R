library(tidyverse)
library(EnhancedVolcano)
library(rstatix)
library(microbiomeMarker)
library(phyloseq)
library(ggpubr)

# analysis on bacteremia vs. no bacteremia 
## Figure 6:
### E: LEfSe Analysis
### F: logistic regression adjusting for other confounding factors
### G: SCFA or Bile acids comparison
### H: Bile acid ratio comparison
## Supp Fig 10:
### A: taxonomy plot
### B: volcano plot of normalized metabolites
### C, D: additional metabolite comparison 

# load source code for plot -----------------------------------------------

source("./code/getRdpPal.R")
source("./code/utility.source.R")

# load in data ------------------------------------------------------------

metab <- read_csv("./data/LD847_HD22.quant.meta.csv") %>% 
  filter(!is.na(bacteremia)) %>% 
  mutate(bacteremia = factor(bacteremia,
                             levels = c(1,0),
                             labels = c("bacteremia","no bacteremia")))

mpa <- read_csv("./data/LD847_HD22.metaphlan.csv") %>% 
  filter(seq_id %in% metab$seq_id)

qual <- read_csv("./data/LD847_HD22.qual.metabolomics.csv") %>% 
  filter(seq_id %in% metab$seq_id) %>% 
  left_join(metab %>% 
              select(seq_id, ID, bacteremia))

ratio <- read_csv("./data/LD847.ratios.quant.metabolomics.csv") %>% 
  filter(seq_id %in% metab$seq_id) %>% 
  gather("ratios", "value", -seq_id, -ID) %>% 
  left_join(metab %>% 
              select(seq_id, ID, bacteremia))

# Fig 6, E: LEfSe ---------------------------------------------------------

kmat <- mpa %>%
  transmute(seq_id, taxid, 
            numseqs = relative_abundance) %>%
  pivot_wider(names_from =taxid,
              values_from = numseqs,
              values_fn = sum,
              values_fill = 0) %>%
  column_to_rownames(var = "seq_id")

## taxonomy table needed for phyloseq object --------
tax_df <- mpa %>%
  distinct(taxid, Kingdom, Phylum, Class, Order, Family, Genus, Species)  %>%
  relocate(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  column_to_rownames(var = "taxid") %>%
  as.matrix()


## make into a phyloseq object ---------------------------------------------

krakenphy <- phyloseq(tax_table(tax_df),
                      otu_table(kmat, taxa_are_rows = F),
                      sample_data(metab %>%
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
  
  ggsave(paste0("results/6.E.LD246.bacteremia.lefse_analysis.pdf"),
         height = 6.5, width = 8.2)
  
}

# Fig 6, F: logistic regression -------------------------------------------

bacteremia.logit <- glm(bacteremia ~ bifido.cat + bacteremia.PPI + bacteremia.meldna_final,
                 data = metab %>% 
                   mutate(bifido.cat = factor(bifido.cat, 
                                              levels = c("No expansion", "Expansion")),
                          bacteremia = if_else(bacteremia == "bacteremia",
                                               1, 0)), family = "binomial")

exp(cbind(OR = coef(bacteremia.logit), confint(bacteremia.logit))) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>% 
  left_join(tidy(bacteremia.logit)) %>% 
  write_csv("results/6.F.bacteremia.meld.PPI.adjust.bifido.or.logit.csv")

# Fig 6, G: SCFA or Bile acids comparison ---------------------------------

metab.long <- metab %>%
  select(ID, `3-oxolithocholic acid`:succinate) %>%
  gather("compound","value",-ID) %>%
  right_join(metab %>%
               distinct(ID, bifido_percent = Bifidobacterium,
                        bacteremia)) 

sbp.quant.pval <- metab.long %>% 
  group_by(compound) %>% 
  wilcox_test(value ~ bacteremia) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  add_significance() %>% 
  add_y_position() %>% 
  mutate(y.position = log10(y.position) * 1.25,
         compound = factor(compound,
                           levels = c("acetate", "butyrate", "propionate", "succinate",
                                      "cholic acid","glycocholic acid", "taurocholic acid",
                                      "lithocholic acid", "deoxycholic acid", "isodeoxycholic acid",
                                      "3-oxolithocholic acid", "alloisolithocholic acid"))) 

subtext <- metab %>%
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
  ggpubr::stat_pvalue_manual(sbp.quant.pval, label = "p") +
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
  scale_y_log10(expand = expansion(mult = c(0,0.1)))

ggsave(file.path("results/6.G.LD246.bacteremia.metabolites.quant.pdf"), 
       width = 12.5, height = 9.3)

# Fig 6, H: Bile acid ratios --------------------------------------------------

ratio.pval <- ratio %>% 
  group_by(ratios) %>% 
  wilcox_test(value ~ bacteremia) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  add_significance() %>% 
  add_y_position() %>% 
  mutate(y.position = log10(y.position) * 1.25) 

ratio %>% 
  mutate(modified_diagnosis = factor(bacteremia, 
                                     levels = c("no bacteremia", "bacteremia"))) %>% 
  ggplot(aes(modified_diagnosis, value, color = modified_diagnosis)) +
  geom_boxplot(outlier.shape = NA,  alpha = 0.75) +
  geom_jitter(alpha = 0.65) +
  scale_color_manual(values =  paletteer::paletteer_d("ggthemes::Tableau_10")[c(2,1)] ) +
  facet_wrap(. ~ ratios, scales = "free",
             ncol = 3) +
  scale_y_log10(expand = expansion(mult = c(0,0.1))) +
  ggpubr::stat_pvalue_manual(ratio.pval, label = "p") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(face = "bold", size = 13.5),
        strip.background.x = element_blank(),
        axis.text = element_text(size = 11.5),
        axis.title = element_text(size = 12.5)) +
  labs(x = "", y = "Bile Acid ratio\n(log transformed)",
       caption = subtext, 
       title = paste0("Bile Acid ratio", " comparison ",
                      "between bacteremia vs. non-bacteremia"))

ggsave("results/6.H.LD246.bacteremia_BAratios.boxplot.pdf", 
       width = 11.5, height = 9.5) 

# Supp Fig 10, A: taxonomy barplot ---------------------------------------------

t <- mpa %>%
  add_count(seq_id, wt = relative_abundance, name = "totalAbd") %>%
  mutate(pctseqs = relative_abundance/totalAbd)  %>%
  dplyr::count(seq_id, taxid, Kingdom, Phylum, Class, Order, Family,
               Genus, Species, wt = pctseqs, name = "pctseqs") %>%
  add_count(seq_id, wt = pctseqs, name = "Total") %>%
  mutate(pctseqs = pctseqs/Total)  %>%
  mutate(genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus)) %>%
  left_join(metab %>% select(seq_id, ID)) 

## stacked barplot ---------------------------------------------------------

pal <- getRdpPal(t)

taxgg <- t %>%
  left_join(metab %>% distinct(seq_id, ID, 
                              bifido_percent = Bifidobacterium, bacteremia)) %>%
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
  scale_y_continuous(breaks = c(seq(0, 1, 0.25), 0.9),
                     expand = expansion(mult = c(0.005,0.005)))

taxgg

sbpgg <- metab %>%
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

pdf("results/10S.A.LD246.bacteremia.metaphlan.pdf",
    height = 5.25, width = 12.5)
gg.stack(taxgg, sbpgg, 
         heights = c(1,0.13),
         newpage = F)
dev.off()

# Supp Fig 10, B: volcano plot --------------------------------------------

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

ggsave(plot = volcano, file.path("results","10S.B.LD246.volcano_metab_bacteremia.qual.pdf"),
       width = 12, height = 6, units = "in")
