library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
library(paletteer)

# Figure 1: 
## A. taxonomy + metabolite plots
## B. selected bacteria boxplot
## C. selected metabolite boxplot
# - between Low, Medium, High vs. Healthy Donor

# load source code for plot -----------------------------------------------

source("./code/getRdpPal.R")

# load in data ------------------------------------------------------------

metab <- read_csv("./data/LD847_HD22.quant.meta.csv") %>% 
  mutate(diversity = factor(diversity, 
                            levels = c("Low", "Medium", "High", "Healthy\nDonor")))

mpa <- read_csv("./data/LD847_HD22.metaphlan.csv")

# A. taxonomy + metabolite plots -------------------------------------------

bk <- mpa %>%
  add_count(seq_id, wt = relative_abundance, name = "totalAbd") %>%
  mutate(pctseqs = relative_abundance/totalAbd)  %>%
  dplyr::count(seq_id, taxid, Kingdom, Phylum, Class, Order, Family,
               Genus, wt = pctseqs, name = "pctseqs") %>%
  add_count(seq_id, wt = pctseqs, name = "Total") %>%
  mutate(pctseqs = pctseqs/Total)  %>%
  mutate(genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus))

## taxonomy plot ---------------------------------------------------------

pal <- getRdpPal(bk)

ggdf <- bk %>%
  group_by(seq_id, Kingdom, Phylum,Class,Order,Family,Genus,genLab) %>%
  summarize(pctseqs=sum(pctseqs)) %>%
  ungroup() %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus,genLab) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus))) %>%
  group_by(seq_id) %>%
  arrange(Genus) %>%
  mutate(cum.pct = cumsum(pctseqs),
         y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>%
  ungroup() %>%
  dplyr::select(-cum.pct)

taxgg <- ggdf %>%
  left_join(metab %>% select(seq_id, cohort, invSimp, diversity)) %>% 
  ggplot(aes(x=reorder(seq_id,invSimp), y=pctseqs)) +
  geom_bar(aes(fill=Genus),stat="identity") +
  # geom_text(aes(y=1-y.text,x=ID,label=tax.label),angle=90,
  #           lineheight=0.6,size=2.2) +
  scale_fill_manual(values=pal) +
  theme_bw() +
  facet_grid(.~ diversity, scales = "free", space = "free_x")+
  theme(legend.position = "none",
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5),
        axis.text.x = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        strip.background = element_rect(fill = "white", colour="white"),
        strip.text.x = element_text(face = "bold", size = 12),
        axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.005,0.005)))+
  labs(x="", y= "Bacterial Abund.")

taxgg

## scfa panel --------------------------------------------------------------

scfametab <- c("acetate","butyrate","propionate")

scfagg <- metab %>%
  select(seq_id, ID, cohort, invSimp, diversity, all_of(scfametab)) %>% 
  gather("compound", "value", -c(seq_id, ID, cohort, invSimp, diversity)) %>% 
  mutate(compound = factor(compound, levels = scfametab),
         biletype=case_when(
           compound %in% c("invSimp") ~ "Inverse Simpson",
           compound %in% c("cholic acid") ~ "primary",
           compound %in% c("lithocholic acid","deoxycholic acid","isodeoxycholic acid") ~ "secondary",
           compound %in% c("glycocholic acid","taurocholic acid") ~ "conjugated",
           compound %in% c("propionate","acetate","butyrate") ~ "short chain fatty acid",
           compound %in% c("succinate") ~ "succinate",
           TRUE ~ "tertiary")) %>%
  ggplot(aes(x=reorder(seq_id, invSimp), y=value)) +
  geom_col(aes(fill = biletype)) +
  scale_fill_manual(values = "#00468BFF") +
  facet_grid(compound ~ diversity, scale = "free", space = "free_x") +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5.5),
        strip.text.y = element_text(angle = 0, face = "bold", hjust = 0),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0.005,0.005))) +
  labs(x = "", y = "Short Chain Fatty Acid Conc.")

scfagg

## bile acid panel ---------------------------------------------------------

bilemetab <- c("cholic acid",
               "glycocholic acid","taurocholic acid",
               "lithocholic acid","deoxycholic acid",
               "alloisolithocholic acid", "3-oxolithocholic acid")

bikegg <-  metab %>%
  select(seq_id, ID, cohort, invSimp, diversity, all_of(bilemetab)) %>% 
  gather("compound", "value", -c(seq_id, ID, cohort, invSimp, diversity)) %>% 
  mutate(compound = factor(compound, levels = bilemetab),
         biletype=case_when(
           compound %in% c("invSimp") ~ "Inverse Simpson",
           compound %in% c("cholic acid") ~ "primary",
           compound %in% c("lithocholic acid","deoxycholic acid","isodeoxycholic acid") ~ "secondary",
           compound %in% c("glycocholic acid","taurocholic acid") ~ "conjugated",
           compound %in% c("propionate","acetate","butyrate") ~ "short chain fatty acid",
           compound %in% c("succinate") ~ "succinate",
           TRUE ~ "tertiary"),
         value = if_else(value <=0, 0, value)) %>%
  ggplot(aes(x=reorder(seq_id,invSimp), y=value)) +
  geom_col(aes(fill = biletype)) +
  facet_grid(compound ~ diversity, scale = "free", space = "free_x") +
  theme_bw()+
  scale_fill_manual(values = c("#ED0000FF", "#42B540FF", "#0099B4FF", "#FDAF91FF")) +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.1),
                     expand = expansion(mult = c(0.005,0.005))) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5.5),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold", hjust = 0),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        strip.background = element_blank()) +
  labs(x = "", y = "Bile Acid Conc.\n(log10 transformed)")

bikegg

## inverse Simpson ---------------------------------------------------------

invgg <- metab %>%
  ggplot(aes(x = reorder(seq_id, invSimp), invSimp)) +
  geom_col(alpha = 0.65) +
  facet_grid(. ~ diversity, scale = "free_x", space = "free_x") +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5.5),
        axis.ticks.x = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold", hjust = 0),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        strip.background = element_blank()) +
  labs(x = "", y = "Inverse\nSimpson") +
  geom_hline(aes(yintercept = 5), linetype = 2, color = "red") +
  scale_y_continuous(breaks = c(0,5,10,20, 30),
                     expand = expansion(mult = c(0.005,0.005)))

## combine all plots -------------------------------------------------------

pdf(file.path("./results/1.A.LD847_HD22.sepmetab.metaphlan.pdf"), 
    height = 8.25, width = 13.5)
gg.stack(taxgg, invgg, scfagg, bikegg,
         heights = c(1.7, 0.35, 0.85,2.25),
         newpage = F)
dev.off()

# B. selected bacteria boxplot ---------------------------------------------

my.comps <- list(c("High","Healthy\nDonor"),
                 c("Medium","Healthy\nDonor"),
                 c("Low","Healthy\nDonor"))

sel.bact <- c("Bacteroidetes", "Bifidobacterium",
              "Lachnospiraceae", "Oscillospiraceae",
              "Enterococcus", "Proteobacteria")

bactabd <- metab %>%
  select(seq_id, all_of(sel.bact)) 

sel.bact.df <- bactabd %>% 
  gather("taxon", "pctseqs", -seq_id) %>% 
  left_join(metab %>% select(seq_id, invSimp, diversity)) %>% 
  mutate(taxon = factor(taxon, 
                        levels=sel.bact) ) 

sel.bact.pdf <- sel.bact.df %>% 
  group_by(taxon) %>% 
  wilcox_test(pctseqs ~ diversity,
              comparisons = my.comps,
              p.adjust.method = "BH") %>% 
  mutate(m.padj = p.adjust(p, method = "BH")) %>% 
  add_significance(p.col = "m.padj") %>% 
  view

sel.bact.df %>% 
  ggplot(aes(diversity, pctseqs)) + 
  geom_boxplot(aes(fill = diversity), width = 0.65, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(aes(color = diversity), width = 0.25,alpha = 0.55) +
  paletteer::scale_color_paletteer_d("rockthemes::facelift", 
                                     direction = -1) +
  paletteer::scale_fill_paletteer_d("rockthemes::facelift",
                                    direction = -1) +
  theme_bw() +
  facet_wrap(. ~ taxon, scales = "free", ncol = 2) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             comparisons = my.comps, size = 3.5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  labs( x= "Diveristy", y = "Bacteria Abundance") +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))

ggsave("./results/1.B.LD847_HD22.bacteriumAbd.bwplt.pdf",
       height = 8.15, width = 5.25)

# C. selected metabolite boxplot ------------------------------------------

sel.comps <- c("acetate","butyrate","propionate",
               "taurocholic acid",
               "lithocholic acid",
               "alloisolithocholic acid")

sel.metab.df <- metab %>%
  select(seq_id, ID, cohort, invSimp, diversity, all_of(sel.comps)) %>% 
  gather("compound", "value", -c(seq_id, ID, cohort, invSimp, diversity)) %>% 
  mutate(compound = factor(compound, levels = sel.comps),
         value = if_else(value < 0, 0, value))

sel.metab.pdf <- sel.metab.df %>% 
  group_by(compound) %>% 
  wilcox_test(value ~ diversity,
              comparisons = my.comps,
              p.adjust.method = "BH") %>% 
  mutate(m.padj = p.adjust(p, method = "BH")) %>% 
  add_significance(p.col = "m.padj") %>% 
  view

sel.metab.df %>%
  ggplot(aes(x=diversity, y=value)) +
  geom_boxplot(aes(fill = diversity), width = 0.65, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(aes(color = diversity), width = 0.25,alpha = 0.55) +
  paletteer::scale_color_paletteer_d("rockthemes::facelift", 
                                     direction = -1) +
  paletteer::scale_fill_paletteer_d("rockthemes::facelift",
                                    direction = -1) +
  theme_bw() +
  facet_wrap(. ~ compound, scales = "free", ncol = 2, dir="v") +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             comparisons = my.comps, size = 3.5,
                             p.adjust.methods = "BH") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  labs( x= "Diveristy", y = "Metabolite Conc.") +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(1),
                     # labels = trans_format("log10", math_format(10^.x)),
                     labels = function(x) format(x, scientific = TRUE),
                     expand = expansion(mult = c(0,0.1)))
# scale_y_log10(expand = expansion(mult = c(0,0.1)),
#                    breaks = trans_breaks("log10", function(x) 10^x),
#                    labels = trans_format("log10", math_format(10^.x))
#                    ) # +
# annotation_logticks() 

ggsave("results/1.C.LD847_HD22.sel.metab.boxplot.pdf",
       height = 8.15, width = 5.25)
