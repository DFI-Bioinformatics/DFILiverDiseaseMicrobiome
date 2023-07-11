library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
library(paletteer)
library(EnhancedVolcano)

# Figure 4: Bifido expansion vs. no expansion comparison in the presence of lactulose 
## A. volcano plot to show differentially abundant normalized metabolites
## C. boxplot of three chosen SCFAs and three bile acids
## D. boxplot of bile acids ratios

# Figure S4: Bifido expansion vs. no expansion comparison in the absence of lactulose 
## A. volcano plot to show differentially abundant normalized metabolites
## B. boxplot of three chosen SCFAs and three bile acids
## C. boxplot of bile acids ratios

# load data ---------------------------------------------------------------

metab <- read_csv("./data/LD847_HD22.quant.meta.csv") %>% 
  filter(fig3.lactu == "Yes", cohort == "Liver\nDisease") %>% 
  mutate(diversity = factor(diversity, 
                            levels = c("Low", "Medium", "High", "Healthy\nDonor")))

ratios <- read_csv("./data/LD847.ratios.quant.metabolomics.csv") %>% 
  filter(seq_id %in% metab$seq_id) %>% 
  select(seq_id, ID, ratio_CA_GCA, ratio_DCA_CA,ratio_DCA_GCA) %>% 
  gather("ratios", "value", -seq_id, -ID) %>% 
  left_join(metab %>% 
              transmute(seq_id, 
                        bifido_cat = if_else(Bifidobacterium < 0.1,
                                                     "Less Bifido\n(<10%)",
                                                     "More Bifido\n(>=10%)"), 
                        lactulose)) %>% 
  mutate(ratios = factor(ratios,
                         levels = c("ratio_CA_GCA", "ratio_DCA_CA","ratio_DCA_GCA"),
                         labels = c("CA:GCA ratio",
                                    "DCA:CA ratio",
                                    "DCA:GCA ratio")))

qual <- read_csv("./data/LD847_HD22.qual.metabolomics.csv") %>% 
  right_join(metab %>% 
               transmute(ID, seq_id, bifido_cat = if_else(Bifidobacterium < 0.1,
                                                          "Less Bifido\n(<10%)",
                                                          "More Bifido\n(>=10%)"), lactulose)) 

# set up hard-coded variables ---------------------------------------------

sel.comps <- c("acetate","butyrate","propionate",
               "taurocholic acid",
               "lithocholic acid",
               "alloisolithocholic acid")

bifido_pal = c("#5050FFFF","#CE3D32FF")

# Fig 4, in the presence of lactulose -------------------------------------

## A. volcano plot  -------------------------------------------------------
### less Bifido as the reference group
lactyes_log2fc <- qual %>%
  filter(lactulose == "Yes") %>%
  group_by(compound, bifido_cat) %>%
  summarise(med = median(value)) %>%
  mutate(med = ifelse(med == 0, 0.0001, med)) %>%
  spread(bifido_cat, med) %>%
  mutate(log2fc = log2(`More Bifido\n(>=10%)`/`Less Bifido\n(<10%)`))

lactyes_pval <- qual %>%
  filter(lactulose == "Yes") %>%
  group_by(compound) %>%
  rstatix::wilcox_test(value ~ bifido_cat) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

lactyes_tot <- left_join(lactyes_log2fc, lactyes_pval) %>%
  column_to_rownames(var = "compound")

#### p-adjusted -------
set.seed(123456)

xylims = ceiling(max(lactyes_tot$log2fc))
plims <- ceiling(max(-log10(lactyes_tot$p)))

volcano <-
  EnhancedVolcano(lactyes_tot,
                  lab = rownames(lactyes_tot),
                  title = 'FDR corrected: More versus Less Bifido (cutoff = 10%)',
                  y = "p.adj",
                  x = "log2fc",
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 3.0,
                  labSize = 2.5,
                  xlim = c(-xylims,xylims),
                  ylim = c(0, plims + 1),
                  col=paletteer::paletteer_d("ggthemes::fivethirtyeight")[c(3,2,5,4)],
                  colAlpha = 0.65,
                  legendPosition = "right",
                  legendLabels = c(expression(p > 0.05*";" ~ Log[2] ~ FC < "\u00B1"*1),
                                   expression(p > 0.05*";" ~ Log[2] ~ FC >= "\u00B1"*1),
                                   expression(p <= 0.05*";" ~ Log[2] ~ FC < "\u00B1"*1),
                                   expression(p <= 0.05*";" ~ Log[2] ~ FC >= "\u00B1"*1)),
                  legendLabSize = 14,
                  drawConnectors = T,
                  widthConnectors = 0.25,
                  maxoverlapsConnectors = Inf,
                  arrowheads = F,
                  gridlines.minor = F,
                  gridlines.major = F) +
  labs(subtitle = "Lactuolose == Yes",
       y = expression( -Log[10] ~ P.adj)) +
  annotate("text", x = 0.65*xylims, y = plims + 0.5, label = "More Bifido",
           size = 6, color = ggsci::pal_igv("default")(2)[2]) +
  annotate("rect", xmin = 1.05, xmax = Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = ggsci::pal_igv("default")(2)[2]) +
  annotate("text", x = -0.65*xylims, y = plims + 0.5, label = "Less Bifido",
           size = 6, color = ggsci::pal_igv("default")(2)[1]) +
  annotate("rect", xmin = -1.05, xmax = -Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = ggsci::pal_igv("default")(2)[1]) +
  guides(color = guide_legend(nrow = 4),
         shape = guide_legend(nrow = 4))

volcano

ggsave(plot = volcano,
       file.path("results","4.A.LD262.volcano_metab_bifido.lactuloseYes.qual_padjusted.pdf"),
       width = 16, height = 8, units = "in")

## C. boxplot of chosen SCFAs and bile acids -------------------------------

sel.comps.df <- metab %>% 
  select(seq_id, all_of(sel.comps)) %>% 
  gather("compound", "value", -seq_id) %>% 
  filter(compound %in% sel.comps) %>% 
  mutate(compound = factor(compound, levels = sel.comps),
         value = if_else(value < 0, 0, value)) %>% 
  left_join(metab %>% 
              transmute(seq_id, 
                        bifido_cat = if_else(Bifidobacterium < 0.1,
                                                  "Less Bifido\n(<10%)",
                                                  "More Bifido\n(>=10%)"), lactulose)) %>% 
  filter(lactulose == "Yes")

sel.comps.sub <- sel.comps.df %>% 
  distinct(seq_id, bifido_cat) %>% 
  group_by(bifido_cat) %>% 
  summarise(n=n()) %>% 
  mutate(bifido_cat = gsub("\n"," ",bifido_cat),
         label = paste0(bifido_cat, "=",n)) 

sel.comps.str <- paste0(sel.comps.sub$label, collapse = "\n")

sel.metab.pdf <- sel.comps.df %>% 
  group_by(compound) %>% 
  wilcox_test(value ~ bifido_cat) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  add_significance(p.col = "p.adj") %>% 
  add_xy_position("bifido_cat") %>% 
  mutate(y.position = log10(y.position)*1.25)

sel.comps.df %>%
  ggplot(aes(x=bifido_cat, y=value)) +
  geom_boxplot(aes(fill = bifido_cat), 
               width = 0.65, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(aes(color = bifido_cat), width = 0.25,alpha = 0.55) +
  scale_fill_manual(values = bifido_pal) +
  scale_color_manual(values = bifido_pal) +
  theme_bw() +
  facet_wrap(. ~ compound, scales = "free", ncol = 3) +
  # stat_compare_means(label = "p.signif") +
  stat_pvalue_manual(sel.metab.pdf, label = "p.adj.signif", tip.length = 0.015) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  labs( x= "", y = "SCFA & Bile Acid Conc.\n(log transformed)",
        caption = sel.comps.str,
        subtitle = "Lactulose = Yes") +
  scale_y_log10(
    labels = scales::scientific,
    breaks = c(0, 0.0001, 0.001 , 0.01, 0.1, 1, 10, 100, 1000),
    expand = expansion(mult = c(0., 0.1))
  )
# scale_y_log10(expand = expansion(mult = c(0,0.1)),
#                    breaks = trans_breaks("log10", function(x) 10^x),
#                    labels = trans_format("log10", math_format(10^.x))
#                    ) # +
# annotation_logticks() 

ggsave("results/4.C.LD262.lactuloseYes.sel.metab.boxplot.pdf",
       height = 6.22, width = 9)

## D. boxplot of bile acids ratios ----------------------------------------

ratio.sub <- ratios %>% 
  filter(lactulose == "Yes") %>% 
  distinct(seq_id, bifido_cat) %>% 
  group_by(bifido_cat) %>% 
  summarise(n=n()) %>% 
  mutate(bifido_cat = gsub("\n"," ",bifido_cat),
         label = paste0(bifido_cat, "=",n)) 

ratio.str <- paste0(ratio.sub$label, collapse = "\n")

ratio.pdf <- ratios %>% 
  filter(lactulose == "Yes") %>% 
  group_by(ratios) %>% 
  wilcox_test(value ~ bifido_cat) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  add_significance(p.col = "p.adj") %>% 
  add_xy_position("bifido_cat") %>% 
  mutate(y.position = log10(y.position)*1.25)

ratios %>% 
  filter(lactulose == "Yes") %>% 
  ggplot(aes(x=bifido_cat, y=value)) +
  geom_boxplot(aes(fill = bifido_cat), 
               width = 0.65, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(aes(color = bifido_cat), width = 0.25,alpha = 0.55) +
  scale_fill_manual(values = bifido_pal) +
  scale_color_manual(values = bifido_pal) +
  theme_bw() +
  facet_wrap(. ~ ratios, scales = "free", ncol = 3) +
  # stat_compare_means(label = "p.signif") +
  stat_pvalue_manual(ratio.pdf, label = "p.adj.signif", tip.length = 0.015) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  labs( x= "", y = "Bile Acid ratios\n(log transformed)",
        caption = ratio.str,
        subtitle = "Lactulose = Yes") +
  scale_y_log10(
    labels = scales::scientific,
    # breaks = c(0, 0.0001, 0.001 , 0.01, 0.1, 1, 10, 100, 1000),
    expand = expansion(mult = c(0., 0.1))
  )

ggsave("results/4.D.LD262.lactuloseYes.bifido.ratio.pdf",
       width = 9.75, height = 4.75)

# Fig S4, in the absence of lactulose -------------------------------------

## A. volcano plot  -------------------------------------------------------
### less Bifido as the reference group
lactno_log2fc <- qual %>%
  filter(lactulose == "No") %>%
  group_by(compound, bifido_cat) %>%
  summarise(med = median(value)) %>%
  mutate(med = ifelse(med == 0, 0.0001, med)) %>%
  spread(bifido_cat, med) %>%
  mutate(log2fc = log2(`More Bifido\n(>=10%)`/`Less Bifido\n(<10%)`))

lactno_pval <- qual %>%
  filter(lactulose == "No") %>%
  group_by(compound) %>%
  rstatix::wilcox_test(value ~ bifido_cat) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

lactno_tot <- left_join(lactno_log2fc, lactno_pval) %>%
  column_to_rownames(var = "compound")

# Volcano Plot
set.seed(123456)

xylims = ceiling(max(lactno_tot$log2fc))
plims <- ceiling(max(-log10(lactno_tot$p)))

#### p-adjusted -------
volcano <-
  EnhancedVolcano(lactno_tot,
                  lab = rownames(lactno_tot),
                  # title = 'FDR corrected: More versus Less Bifido (cutoff = 10%)',
                  title = 'FDR corrected: More versus Less Bifido (cutoff = 10%)',
                  y = "p.adj",
                  x = "log2fc",
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 3.0,
                  labSize = 2.5,
                  xlim = c(-xylims,xylims),
                  ylim = c(0, plims + 1),
                  col=paletteer::paletteer_d("ggthemes::fivethirtyeight")[c(3,2,5,4)],
                  colAlpha = 0.65,
                  legendPosition = "right",
                  legendLabels = c(expression(p > 0.05*";" ~ Log[2] ~ FC < "\u00B1"*1),
                                   expression(p > 0.05*";" ~ Log[2] ~ FC >= "\u00B1"*1),
                                   expression(p <= 0.05*";" ~ Log[2] ~ FC < "\u00B1"*1),
                                   expression(p <= 0.05*";" ~ Log[2] ~ FC >= "\u00B1"*1)),
                  legendLabSize = 14,
                  drawConnectors = T,
                  widthConnectors = 0.25,
                  arrowheads = F,
                  maxoverlapsConnectors = Inf,
                  gridlines.minor = F,
                  gridlines.major = F) +
  labs(subtitle = "Lactuolose == No",
       y = expression( -Log[10] ~ P.adj)) +
  annotate("text", x = 0.65*xylims, y = plims + 0.5, label = "More Bifido",
           size = 6, color = ggsci::pal_igv("default")(2)[2]) +
  annotate("rect", xmin = 1.05, xmax = Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = ggsci::pal_igv("default")(2)[2]) +
  annotate("text", x = -0.65*xylims, y = plims + 0.5, label = "Less Bifido",
           size = 6, color = ggsci::pal_igv("default")(2)[1]) +
  annotate("rect", xmin = -1.05, xmax = -Inf, ymin = -log(0.05, base = 10),
           ymax = Inf,
           alpha = .1, fill = ggsci::pal_igv("default")(2)[1]) +
  guides(color = guide_legend(nrow = 4),
         shape = guide_legend(nrow = 4))

volcano

ggsave(plot = volcano,
       file.path("results","4S.A.LD262.volcano_metab_bifido.lactuloseNo.qual_padjusted.pdf"),
       width = 16, height = 8, units = "in")

## B. boxplot of chosen SCFAs and bile acids -------------------------------

sel.comps.df <- metab %>% 
  select(seq_id, all_of(sel.comps)) %>% 
  gather("compound", "value", -seq_id) %>% 
  filter(compound %in% sel.comps) %>% 
  mutate(compound = factor(compound, levels = sel.comps),
         value = if_else(value < 0, 0, value)) %>% 
  left_join(metab %>% 
              transmute(seq_id, bifido_cat = if_else(Bifidobacterium < 0.1,
                                                     "Less Bifido\n(<10%)",
                                                     "More Bifido\n(>=10%)"), lactulose)) %>% 
  filter(lactulose == "No")

sel.comps.sub <- sel.comps.df %>% 
  distinct(seq_id, bifido_cat) %>% 
  group_by(bifido_cat) %>% 
  summarise(n=n()) %>% 
  mutate(bifido_cat = gsub("\n"," ",bifido_cat),
         label = paste0(bifido_cat, "=",n)) 

sel.comps.str <- paste0(sel.comps.sub$label, collapse = "\n")

sel.metab.pdf <- sel.comps.df %>% 
  group_by(compound) %>% 
  wilcox_test(value ~ bifido_cat) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  add_significance(p.col = "p.adj") %>% 
  add_xy_position("bifido_cat") %>% 
  mutate(y.position = log10(y.position)*1.25)

sel.comps.df %>%
  ggplot(aes(x=bifido_cat, y=value)) +
  geom_boxplot(aes(fill = bifido_cat), 
               width = 0.65, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(aes(color = bifido_cat), width = 0.25,alpha = 0.55) +
  scale_fill_manual(values = bifido_pal) +
  scale_color_manual(values = bifido_pal) +
  theme_bw() +
  facet_wrap(. ~ compound, scales = "free", ncol = 3) +
  # stat_compare_means(label = "p.signif") +
  stat_pvalue_manual(sel.metab.pdf, label = "p.adj.signif", tip.length = 0.015) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  labs( x= "", y = "SCFA & Bile Acid Conc.\n(log transformed)",
        caption = sel.comps.str,
        subtitle = "Lactulose = No") +
  scale_y_log10(
    labels = scales::scientific,
    breaks = c(0, 0.0001, 0.001 , 0.01, 0.1, 1, 10, 100, 1000),
    expand = expansion(mult = c(0., 0.1))
  )
# scale_y_log10(expand = expansion(mult = c(0,0.1)),
#                    breaks = trans_breaks("log10", function(x) 10^x),
#                    labels = trans_format("log10", math_format(10^.x))
#                    ) # +
# annotation_logticks() 

ggsave("results/4S.B.LD262.lactuloseNo.sel.metab.boxplot.pdf",
       height = 6.22, width = 9)

## C. boxplot of bile acids ratios ----------------------------------------

ratio.sub <- ratios %>% 
  filter(lactulose == "No") %>% 
  distinct(seq_id, bifido_cat) %>% 
  group_by(bifido_cat) %>% 
  summarise(n=n()) %>% 
  mutate(bifido_cat = gsub("\n"," ",bifido_cat),
         label = paste0(bifido_cat, "=",n)) 

ratio.str <- paste0(ratio.sub$label, collapse = "\n")

ratio.pdf <- ratios %>% 
  filter(lactulose == "No") %>% 
  group_by(ratios) %>% 
  wilcox_test(value ~ bifido_cat) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  add_significance(p.col = "p.adj") %>% 
  add_xy_position("bifido_cat") %>% 
  mutate(y.position = log10(y.position)*1.25)

ratios %>% 
  filter(lactulose == "No") %>% 
  ggplot(aes(x=bifido_cat, y=value)) +
  geom_boxplot(aes(fill = bifido_cat), 
               width = 0.65, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(aes(color = bifido_cat), width = 0.25,alpha = 0.55) +
  scale_fill_manual(values = bifido_pal) +
  scale_color_manual(values = bifido_pal) +
  theme_bw() +
  facet_wrap(. ~ ratios, scales = "free", ncol = 3) +
  # stat_compare_means(label = "p.signif") +
  stat_pvalue_manual(ratio.pdf, label = "p.adj.signif", tip.length = 0.015) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  labs( x= "", y = "Bile Acid ratios\n(log transformed)",
        caption = ratio.str,
        subtitle = "Lactulose = No") +
  scale_y_log10(
    labels = scales::scientific,
    # breaks = c(0, 0.0001, 0.001 , 0.01, 0.1, 1, 10, 100, 1000),
    expand = expansion(mult = c(0., 0.1))
  )

ggsave("results/4S.C.LD262.lactuloseNo.bifido.ratio.pdf",
       width = 9.75, height = 4.75)
