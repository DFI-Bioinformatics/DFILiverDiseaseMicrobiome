library(tidyverse)
library(EnhancedVolcano)
library(rstatix)
library(microbiomeMarker)
library(phyloseq)
library(ggpubr)

# analysis on SBP vs. No SBP 
## Figure 6:
### A: LEfSe Analysis
### B: logistic regression adjusting for other confounding factors
### C: SCFA or Bile acids comparison
### D: Bile acid ratio comparison
## Supp Fig 9:
### A: taxonomy plot
### B: volcano plot of normalized metabolites
### C, D: additional metabolite comparison 

# load source code for plot -----------------------------------------------

source("./code/getRdpPal.R")
source("./code/utility.source.R")

# load in data ------------------------------------------------------------

metab <- read_csv("./data/LD847_HD22.quant.meta.csv") %>% 
  filter(!is.na(sbp))

mpa <- read_csv("./data/LD847_HD22.metaphlan.csv") %>% 
  filter(seq_id %in% metab$seq_id)

qual <- read_csv("./data/LD847_HD22.qual.metabolomics.csv") %>% 
  filter(seq_id %in% metab$seq_id) %>% 
  left_join(metab %>% 
              select(seq_id, ID, sbp_diag))

ratio <- read_csv("./data/LD847.ratios.quant.metabolomics.csv") %>% 
  filter(seq_id %in% metab$seq_id) %>% 
  gather("ratios", "value", -seq_id, -ID) %>% 
  left_join(metab %>% 
              select(seq_id, ID, sbp_diag))

# Fig 6, A: LEfSe ------------------------------------------------------------

## transform into matrix ----------

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
                                    distinct(ID, seq_id, sbp_diag) %>%
                                    column_to_rownames(var = "seq_id")
                      ))


## run lefse ---------------------------------------------------------------

{
  phy_lefse <- run_lefse(
    ps = krakenphy,
    group = "sbp_diag",
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
           ef_lda = if_else(enrich_group == "sbp", ef_lda, -ef_lda),
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
    paletteer::scale_fill_paletteer_d("rcartocolor::Bold") +
    scale_y_continuous(breaks = seq(-4,4,1))+
    coord_flip() +
    labs(y = "LDA Score (log10)", x = "Phylogeny\n", fill = "sbp\nStatus",
         title = "LEfSe Analysis: Grouped by sbp Status") +
    facet_grid("level", scales = "free", space = "free_y")
  
  ggsave(paste0("results/6.A.LD122.sbp.lefse_analysis.pdf"),
         height = 8.5, width = 8.2)
  
}

# Fig 6, B: logistic regression  -------------------------------------------

sbp.logit <- glm(sbp ~ bifido.cat + sbp.PPI + sbp.meldna_final,
                 data = metab %>% mutate(bifido.cat = factor(bifido.cat, levels = c("No expansion", "Expansion"))), family = "binomial")

exp(cbind(OR = coef(sbp.logit), confint(sbp.logit))) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>% 
  left_join(tidy(sbp.logit)) %>% 
  write_csv("results/6.B.sbp.meld.PPI.adjust.bifido.or.logit.csv")

# Fig 6, C: SCFA or Bile acids comparison ----------------------------------

metab.long <- metab %>%
  select(ID, `3-oxolithocholic acid`:succinate) %>%
  gather("compound","value",-ID) %>%
  right_join(metab %>%
               distinct(ID, bifido_percent = Bifidobacterium,
                        sbp_diag)) 

metab.long %>%
  count(ID, sort = T)

sbp.quant.pval <- metab.long %>% 
  group_by(compound) %>% 
  wilcox_test(value ~ sbp_diag) %>% 
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
  count(sbp_diag) %>%
  mutate(lab = paste0(sbp_diag,"=", n)) %>%
  pull(lab) %>%
  paste0(collapse = "; ")

metab.long %>%
  mutate(compound = factor(compound,
                           levels = c("acetate", "butyrate", "propionate", "succinate",
                                      "cholic acid","glycocholic acid", "taurocholic acid",
                                      "lithocholic acid", "deoxycholic acid", "isodeoxycholic acid",
                                      "3-oxolithocholic acid", "alloisolithocholic acid"))) %>%
  ggplot(aes(sbp_diag, value, color = sbp_diag)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_wrap("compound", scales = "free") +
  paletteer::scale_color_paletteer_d("rcartocolor::Bold") +
  theme_bw() +
  stat_pvalue_manual(sbp.quant.pval, label = "p") +
  theme(legend.position = "none",
        strip.text.x = element_text(face = "bold", size = 13.5),
        strip.background.x = element_blank(),
        axis.text = element_text(size = 11.5),
        axis.title = element_text(size = 12.5)) +
  labs(x = "", y = "Metabolites Conc.\n(log transformed)",
       caption = subtext) +
  scale_y_log10(expand = expansion(mult = c(0,0.1)))

ggsave("results/6.C.LD122.sbp.metabolites.quant.pdf", width = 12.5, height = 9.3)

# Fig 6, D: bile acid ratios --------------------------------------------------

ratio.pval <- ratio %>% 
  group_by(ratios) %>% 
  wilcox_test(value ~ sbp_diag) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  add_significance() %>% 
  add_y_position() %>% 
  mutate(y.position = log10(y.position) * 1.25) 

ratio %>% 
  ggplot(aes(sbp_diag, value, color = sbp_diag)) +
  geom_boxplot(outlier.shape = NA,  alpha = 0.75) +
  geom_jitter(alpha = 0.65) +
  paletteer::scale_color_paletteer_d("rcartocolor::Bold") +
  facet_wrap(. ~ ratios, scales = "free",
             ncol = 3) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(face = "bold", size = 13.5),
        strip.background.x = element_blank(),
        axis.text = element_text(size = 11.5),
        axis.title = element_text(size = 12.5)) +
  scale_y_log10(expand = expansion(mult = c(0,0.1))) +
  stat_pvalue_manual(ratio.pval, label = "p") +
  labs(x = "", y = "Bile Acid ratio\n(log transformed)",
       caption = subtext, 
       title = paste0("Bile Acid ratio", " comparison ",
                      "between SBP vs. non-SBP"))

ggsave("results/6.D.LD122.sbp_BAratios.boxplot.pdf", 
       width = 11.5, height = 9.5) 


# Supp Fig 9, A: taxonomy barplot ---------------------------------------------

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

# checks
t %>% 
  count(seq_id, wt = pctseqs)

## stacked tax plot --------------------------------------------------------

pal <- getRdpPal(t)

taxgg <- t %>%
  left_join(metab %>% distinct(seq_id, ID, bifido_percent = Bifidobacterium, sbp_diag)) %>%
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

## sbp yes vs. no ---------------------------------------------------------

sbpgg <- metab %>%
  ggplot(aes(reorder(ID, -Bifidobacterium), "SBP status")) +
  geom_tile(aes(fill = sbp_diag)) +
  # facet_grid(. ~ sbp_diag, scale = "free", space = "free_x") +
  labs(x="", y= "") +
  scale_fill_manual(values = c("no sbp" = "#09ba38", "sbp" = "#db3b3b")) +
  scale_y_discrete(expand = expansion(mult = c(0.005,0.005)))+
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

pdf("results/9S.A.LD122.sbp.metaphlan.sbp.pdf",
    height = 5.25, width = 12.5)
gg.stack(taxgg, sbpgg, 
         heights = c(1,0.08),
         newpage = F)
dev.off()

# Supp Fig 9, B: volcano plot ---------------------------------------------

### no sbp as the reference group
sbp_log2fc <- qual %>%
  group_by(compound, sbp_diag) %>%
  summarise(med = median(value)) %>% 
  mutate(med = ifelse(med == 0, 0.0001, med)) %>%
  spread(sbp_diag, med) %>%
  mutate(log2fc = log2(sbp/`no sbp`))

sbp_pval <- qual %>%
  group_by(compound) %>%
  rstatix::wilcox_test(value ~ sbp_diag) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

sbp_tot <- left_join(sbp_log2fc, sbp_pval) %>%
  column_to_rownames(var = "compound")

# Volcano Plot
set.seed(123456)

xylims = ceiling(max(sbp_tot$log2fc))
plims <- ceiling(max(-log10(sbp_tot$p)))

volcano <-
  EnhancedVolcano(sbp_tot,
                  lab = rownames(sbp_tot),
                  title = 'SBP versus No SBP',
                  y = "p",
                  x = "log2fc",
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 3.0,
                  labSize = 3.5,
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
  annotate("text", x = (xylims-0.5), y = plims + 1, label = "SBP",
           size = 6, color = paletteer::paletteer_d("rcartocolor::Bold")[2]) +
  annotate("rect", xmin = 1.05, xmax = Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = paletteer::paletteer_d("rcartocolor::Bold")[2]) +
  annotate("text", x = -c(xylims-0.5), y = plims + 1, label = "No SBP",
           size = 6, color = paletteer::paletteer_d("rcartocolor::Bold")[1]) +
  annotate("rect", xmin = -1.05, xmax = -Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = paletteer::paletteer_d("rcartocolor::Bold")[1]) +
  guides(color = guide_legend(nrow = 4),
         shape = guide_legend(nrow = 4))

volcano

ggsave(plot = volcano, file.path("results","9S.B.LD122.volcano_metab_sbp.qual.pdf"),
       width = 12, height = 6, units = "in")
