library(tidyverse)
library(EnhancedVolcano)
library(rstatix)

# Figure 4: metabolite comparison between Less Bifido vs. More Bifido ----------
## under condition lactulose / no lactulose

## load in quant data ------------------------------------------------------------

metab <- read_csv("./data/LD850.meta.quant.metabolomics.csv") 

## prepare long quant data format -------------------------------------------

metab.long <- metab %>%
  select(ID, `3-oxolithocholic acid`:succinate) %>%
  gather("compound","value",-ID) %>%
  left_join(metab %>% 
              select(ID, 
                     bifido_class, 
                     lactulose)) 

### in lactulose == "yes" -------------------------------------------

tarcomps <- c("acetate","butyrate","propionate",
              "taurocholic acid", 
              "lithocholic acid",
              "alloisolithocholic acid")

### boxplot 

subtext <- metab.long %>%
  filter(lactulose == "yes",
         compound %in% tarcomps) %>% 
  distinct(ID,bifido_class) %>% 
  count(bifido_class) %>% 
  mutate(lactulose = paste0("lactulose=", bifido_class),
         lab = paste0(gsub("\n","",bifido_class),"=", n)
  ) %>%
  pull(lab) %>%
  paste0(collapse = "\n")

metab.long %>%
  filter(lactulose == "yes",
         compound %in% tarcomps) %>% 
  mutate(compound = factor(compound, levels = tarcomps)) %>% 
  ggplot(aes(bifido_class, value, fill = bifido_class)) +
  geom_boxplot(outlier.shape = NA,  alpha = 0.75) +
  geom_jitter(shape = 21, alpha = 0.65) +
  scale_fill_manual(values = ggsci::pal_igv("default")(2)) +
  facet_wrap("compound", scales = "free") +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.001),
                     labels = scales::scientific) +
  theme_bw() +
  ggpubr::stat_compare_means() +
  theme(legend.position = "none",
        strip.text.x = element_text(face = "bold", size = 13.5),
        strip.background.x = element_blank(),
        axis.text = element_text(size = 11.5),
        axis.title = element_text(size = 12.5)) +
  labs(x = "", y = "SCFA & Bile Acid Conc.\n(log transformed)",
       subtitle = "Lactulose == Yes",
       caption = subtext)

ggsave("./results/4.LD850.metab_bifido.lactuloseYes.pdf",
       width = 12.5, height = 8.5)

### in lactulose == "no" -------------------------------------------

### boxplot 

subtext <- metab.long %>%
  filter(lactulose == "no",
         compound %in% tarcomps) %>% 
  distinct(ID,bifido_class) %>% 
  count(bifido_class) %>% 
  mutate(lactulose = paste0("lactulose=", bifido_class),
         lab = paste0(gsub("\n","",bifido_class),"=", n)
  ) %>%
  pull(lab) %>%
  paste0(collapse = "\n") 

metab.long %>%
  filter(lactulose == "no",
         compound %in% tarcomps) %>% 
  mutate(compound = factor(compound, levels = tarcomps)) %>% 
  ggplot(aes(bifido_class, value, fill = bifido_class)) +
  geom_boxplot(outlier.shape = NA,  alpha = 0.75) +
  geom_jitter(shape = 21, alpha = 0.65) +
  scale_fill_manual(values = ggsci::pal_igv("default")(2)) +
  facet_wrap("compound", scales = "free") +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.001),
                     labels = scales::scientific) +
  theme_bw() +
  ggpubr::stat_compare_means() +
  theme(legend.position = "none",
        strip.text.x = element_text(face = "bold", size = 13.5),
        strip.background.x = element_blank(),
        axis.text = element_text(size = 11.5),
        axis.title = element_text(size = 12.5)) +
  labs(x = "", y = "SCFA & Bile Acid Conc.\n(log transformed)",
       subtitle = "Lactulose == No",
       caption = subtext)

ggsave("./results/4.LD850.metab_bifido.lactuloseNo.pdf",
       width = 12.5, height = 8.5)

## quant ratio boxplot -----------------------------------------------------

ratio <- read_csv("./data/LD850.ratios.quant.metabolomics.csv") %>% 
  left_join(metab %>% 
              select(seq_id, ID, bifido_class, lactulose))

ratio_long <- ratio %>% 
  gather("ratios", "value", -seq_id, -ID, 
         -lactulose,-bifido_class)

## with lactulose -----------------------------------------------

subtext <- ratio_long %>% 
  distinct(seq_id, bifido_class,lactulose) %>% 
  count(bifido_class,lactulose) %>% 
  mutate(lactulose = paste0("lactulose=", lactulose),
         lab = paste0(lactulose, "-", gsub("\n","",bifido_class),"=", n)
  ) %>%
  pull(lab) %>%
  paste0(collapse = "\n") 

ratio_long %>% 
  mutate(lactulose = paste0("lactulose=", lactulose)) %>% 
  ggplot(aes(bifido_class, value, fill = bifido_class)) +
  geom_boxplot(outlier.shape = NA,  alpha = 0.75) +
  geom_jitter(shape = 21, alpha = 0.65) +
  scale_fill_manual(values = ggsci::pal_igv("default")(2)) +
  facet_wrap(. ~ lactulose + ratios, scales = "free",
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
       title = paste0("Bile Acid ratio", " comparison"))

ggsave("./results/4.LD850.bifido.lactulose_BAratios.boxplot.pdf", 
       width = 10.5, height = 17.5) 

## load in qual data -------------------------------------------------------

qual <- read_csv("./data/LD850.qual.metabolomics.csv") %>% 
  left_join(metab %>% 
              select(ID, seq_id, bifido_class, lactulose))

### volcano plot of qual data -----------------------------------------------

#### in lactulose == "yes" -------------------------------------------

### less Bifido as the reference group
lactyes_log2fc <- qual %>%
  filter(lactulose == "yes") %>%
  group_by(compound, bifido_class) %>%
  summarise(med = median(value)) %>%
  mutate(med = ifelse(med == 0, 0.0001, med)) %>%
  spread(bifido_class, med) %>%
  mutate(log2fc = log2(`More Bifido\n(>=10%)`/`Less Bifido\n(<10%)`))

lactyes_pval <- qual %>%
  filter(lactulose == "yes") %>%
  group_by(compound) %>%
  rstatix::wilcox_test(value ~ bifido_class) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

lactyes_tot <- left_join(lactyes_log2fc, lactyes_pval) %>%
  column_to_rownames(var = "compound")

# Volcano Plot
set.seed(123456)

xylims = ceiling(max(lactyes_tot$log2fc))
plims <- ceiling(max(-log10(lactyes_tot$p)))

#### not p-adjusted -------
volcano <-
  EnhancedVolcano(lactyes_tot,
                  lab = rownames(lactyes_tot),
                  title = 'More versus Less Bifido (cutoff = 10%)',
                  y = "p",
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
       y = expression( -Log[10] ~ P)) +
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
       file.path("./results/4.LD850.volcano_metab_bifido.lactuloseYes.qual.pdf"),
       width = 16, height = 8, units = "in")

#### p-adjusted -------
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
       file.path("./results/4.LD850.volcano_metab_bifido.lactuloseYes.qual_padjusted.pdf"),
       width = 16, height = 8, units = "in")

### in lactulose == "no" -------------------------------------------

### less Bifido as the reference group
lactno_log2fc <- qual %>%
  filter(lactulose == "no") %>%
  group_by(compound, bifido_class) %>%
  summarise(med = median(value)) %>%
  mutate(med = ifelse(med == 0, 0.0001, med)) %>%
  spread(bifido_class, med) %>%
  mutate(log2fc = log2(`More Bifido\n(>=10%)`/`Less Bifido\n(<10%)`))

lactno_pval <- qual %>%
  filter(lactulose == "no") %>%
  group_by(compound) %>%
  rstatix::wilcox_test(value ~ bifido_class) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

lactno_tot <- left_join(lactno_log2fc, lactno_pval) %>%
  column_to_rownames(var = "compound")

# Volcano Plot
set.seed(123456)

xylims = ceiling(max(lactno_tot$log2fc))
plims <- ceiling(max(-log10(lactno_tot$p)))

#### not p-adjusted -------
volcano <-
  EnhancedVolcano(lactno_tot,
                  lab = rownames(lactno_tot),
                  # title = 'FDR corrected: More versus Less Bifido (cutoff = 10%)',
                  title = 'More versus Less Bifido (cutoff = 10%)',
                  y = "p",
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
       y = expression( -Log[10] ~ P)) +
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
       file.path("./results/4.LD850.volcano_metab_bifido.lactuloseNo.qual.pdf"),
       width = 16, height = 8, units = "in")

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
       file.path("./results/4.LD850.volcano_metab_bifido.lactuloseNo.qual_padjusted.pdf"),
       width = 16, height = 8, units = "in")


# Figure 5: Proteobacteria, Enterococcus abundance comparison between less and more Bifido groups ------
## under lactulose / no lactulose contions

metab %>% 
  select(seq_id, ID, Enterococcus, Proteobacteria) %>% 
  gather("taxon", "n", -seq_id, -ID) %>% 
  left_join(metab %>% 
              select(seq_id, ID, lactulose, bifido_class)) %>% 
  mutate(lactulose = paste0("lactulose status = ", lactulose)) %>%
  ggplot(aes(bifido_class, n, color = bifido_class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.65) +
  facet_grid(. ~ lactulose + taxon) +
  ggpubr::stat_compare_means(label.y = 7.15) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(face = "bold", size = 13.5),
        strip.background.x = element_blank(),
        axis.text = element_text(size = 11.5),
        axis.title = element_text(size = 12.5))  +
  labs(x = "Bifido", y = "Bacteria Abund.\n(log transformed)") +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.01))

ggsave(file.path("./results/4.LD850.Bifido_Entero_Proteo.lactulose.pdf"),
       width = 12.85, height = 5.3)
