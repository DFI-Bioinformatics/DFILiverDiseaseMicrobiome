library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
library(paletteer)

# Figure 3: 
## A. lactulose taxonomy plot + clinical variables
## B. bacterium comparison
# Figure 5A
## Enterococcus and Proteobacteria comparison between Bifido expanded vs. no expansion.
# Table S1 
## Bifidobacteria relative abundance as response variable
## Enterococcus relative abundance  as response variable

# load source code for plot -----------------------------------------------

source("./code/getRdpPal.R")

# load in data ------------------------------------------------------------

metab <- read_csv("./data/LD847_HD22.quant.meta.csv") %>% 
  filter(fig3.lactu == "Yes") %>% 
  mutate(diversity = factor(diversity, 
                            levels = c("Low", "Medium", "High", "Healthy\nDonor")))

mpa <- read_csv("./data/LD847_HD22.metaphlan.csv") %>% 
  filter(seq_id %in% metab$seq_id)

# Fig 3, A. taxonomy plot -----------------------------------------------------

bk <- mpa %>%
  dplyr::count(seq_id, taxid, Kingdom, Phylum, Class, Order, Family,
               Genus, Species, wt = relative_abundance, name = "pctseqs") %>%
  add_count(seq_id, wt = pctseqs, name = "Total") %>%
  mutate(pctseqs = pctseqs/Total)  %>%
  mutate(genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus))

t <- bk %>%
  left_join(metab %>%
              distinct(seq_id,
                       ID))

pal <- getRdpPal(t)

ggdf <- t %>%
  group_by(ID, Kingdom, Phylum,Class,Order,Family,Genus,genLab) %>%
  summarize(pctseqs=sum(pctseqs)) %>%
  ungroup() %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus,genLab) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus))) %>%
  group_by(ID) %>%
  arrange(Genus) %>%
  mutate(cum.pct = cumsum(pctseqs),
         y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>%
  ungroup() %>%
  dplyr::select(-cum.pct) %>%
  mutate(#tax.label=ifelse(Species=="unclassified",paste(Family,Genus,Species,sep="\n"),paste(Genus,Species,sep="\n")),
    tax.label= ifelse(grepl("unclassified",genLab), Family,genLab),
    tax.label = ifelse(pctseqs >= .1, as.character(tax.label), ""))

sorted.taxgg <- ggdf %>%
  left_join(metab %>%
              distinct(seq_id,
                       ID, lactulose, Bifidobacterium)) %>% 
  ggplot(aes(x=reorder(ID, -Bifidobacterium), y=pctseqs)) +
  geom_bar(aes(fill=Genus),stat="identity") +
  # geom_text(aes(y=1-y.text,x=ID,label=tax.label),angle=90,
  #           lineheight=0.6,size=2.2) +
  scale_fill_manual(values=pal) +
  theme_bw() +
  theme(legend.position = "none",
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5),
        #axis.text.x = element_blank()
        strip.background.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        axis.ticks.x = element_blank()
  ) +
  facet_grid(.~ lactulose, scales = "free", space = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0.005,0.005)))+
  labs(x="", y= "Bacterial Abund.")

sorted.taxgg

## clinical variables ------------------------------------------------------

clin.gg <- metab %>%
  select(ID, abx, laxative, rifaximin, PPI,
         bifido_percent = Bifidobacterium,
         portal_htn, acute_or_chronic) %>%
  gather("var","value", -ID,-bifido_percent) %>%
  mutate(var = factor(var, levels = rev(c("rifaximin", "abx",
                                          "PPI", "laxative",
                                          "acute_or_chronic",
                                          "portal_htn")))) %>%
  left_join(metab %>% select(ID, lactulose)) %>%
  ggplot(aes(reorder(ID,-bifido_percent), var)) +
  geom_tile(aes(fill = value)) +
  theme_bw() +
  theme(legend.position = "right",
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5),
        #axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_text(size = 8.5),
        axis.ticks.x = element_blank(),
        legend.key.height= unit(0.5, 'lines'),
        legend.key.width= unit(0.5, 'lines'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9.5)
  ) +
  facet_grid( . ~ lactulose, scales = "free", space = "free") +
  labs(x="", y= "", fill = "administered") +
  scale_y_discrete(expand = expansion(mult = c(0.005,0.005)))+
  scale_fill_manual(values = c("No" = "#09ba38", "Yes" = "#db3b3b", 
                               "Healthy Donor" = "#09ba38",
                               "Acute" = "#DD5129FF",
                               "Chronic" = "#FAB255FF",
                               "Cirrhosis (Compensated)" = "#43B284FF",
                               "Cirrhosis (Decompensated)" = "#0F7BA2FF"))

clin.gg

## stool cons --------------------------------------------------------------

stool.pal <- paletteer::paletteer_d("nbapalettes::warriors_city2")[c(1,3,2)]

stool.gg <- metab %>%
  select(ID, stool_consistency,
         bifido_percent = Bifidobacterium) %>%
  left_join(metab %>% select(ID, lactulose)) %>%
  ggplot(aes(reorder(ID,-bifido_percent), 1)) +
  geom_tile(aes(fill = stool_consistency)) +
  theme_bw() +
  theme(legend.position = "right",
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5),
        #axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.height= unit(0.5, 'lines'),
        legend.key.width= unit(0.5, 'lines'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9.5)
  ) +
  facet_grid( . ~ lactulose, scales = "free", space = "free") +
  labs(x="", y= "", fill = "stool\nconsistency") +
  scale_y_continuous(expand = expansion(mult = c(0.005,0.005)),
                     breaks = c(1), labels = "Stool\nConsistency") +
  scale_fill_manual(values = stool.pal)

stool.gg

## lactulose dose --------------------------------------------------------------

lact.gg <- metab %>%
  select(ID, oral_lactulose,
         bifido_percent = Bifidobacterium) %>%
  gather("var","value", -ID,-bifido_percent) %>%
  left_join(metab %>% select(ID, lactulose)) %>%
  ggplot(aes(reorder(ID,-bifido_percent), var)) +
  geom_tile(aes(fill = value)) +
  theme_bw() +
  theme(legend.position = "right",
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5),
        #axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.height= unit(0.5, 'lines'),
        legend.key.width= unit(0.5, 'lines'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9.5)
  ) +
  facet_grid( . ~ lactulose, scales = "free", space = "free") +
  labs(x="", y= "", fill = "lactolose dose") +
  scale_y_discrete(expand = expansion(mult = c(0.005,0.005))) +
  scale_fill_binned(trans = yingtools2::log_epsilon_trans(10),
                    type = "viridis",
                    labels = function(x) format(x, scientific = TRUE)) 

lact.gg

## vanA --------------------------------------------------------------------

sorted.toxingg <- metab %>%
              select(seq_id, ID, lactulose,vanA,
                     bifido_percent = Bifidobacterium) %>%
  mutate(log10.var = log10(vanA),
         log10.var = if_else(is.infinite(log10.var),
                             0, log10.var),
         var = "vanA") %>% 
  # summary(log10.var)
  ggplot(aes(reorder(ID,-bifido_percent), var)) +
  geom_tile(aes(fill = log10.var), stat = "identity") +
  scale_fill_gradient2(low = "white",
                       mid= "cyan1",
                       high = "#0c0970",
                       midpoint = 1.5,
                       breaks = seq(-1, 4, by = 1),
                       # labels = c("1", "2", "3", "4", "5", "6", "7"),
                       labels = scales::label_math(),
                       guide = "legend",
                       limits = c(-1, 4)
  ) +
  guides(fill = guide_legend(title = "RPKM",
                             reverse = TRUE)) +
  theme_bw() +
  facet_grid( . ~ lactulose, scales = "free", space = "free") +
  scale_y_discrete(expand = expansion(mult = c(0.005,0.005))) +
  theme(legend.position = "right",
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8.5),
        #axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.height= unit(0.5, 'lines'),
        legend.key.width= unit(0.5, 'lines'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9.5),
  ) +
  labs(x="", y= "", fill = "RPKM")

sorted.toxingg

## combine everything together ---------------------------------------------

pdf(file.path("results","3.A.LD262_HD22.lactulose.pdf"),
    height = 7.75, width = 14.5)
gg.stack(sorted.taxgg, sorted.toxingg, lact.gg, clin.gg, stool.gg, 
         heights = c(1,0.12, 0.12, 0.60, 0.12),
         newpage = F)
dev.off()

# Fig 3, B. bacterium comparison -------------------------------------------------

## lactulose vs no lactulose w/o abx + rifaximin ----------------------------

no_abx_rif.df <- metab %>% 
  filter(abx == "No", rifaximin == "No", cohort == "Liver\nDisease") %>% 
  select(seq_id, lactulose, Bifidobacterium, Proteobacteria, Enterococcus) %>% 
  gather("taxon","pctseqs",-c(seq_id, lactulose)) %>% 
  mutate(lactulose = factor(lactulose,
                            levels = c("Yes","No"),
                            labels = c("Lactulose","No lactulose"))) %>%
  replace_na(list(pctseqs = 0)) 

no_abx_rif.pdf <- no_abx_rif.df %>% 
  group_by(taxon) %>% 
  wilcox_test(pctseqs ~ lactulose) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  add_significance(p.col = "p.adj") %>% 
  add_xy_position("lactulose")

no_abx_rif.subs <- no_abx_rif.df %>% 
  distinct(seq_id,lactulose) %>% 
  group_by(lactulose) %>% 
  summarise(n = n()) %>% 
  mutate(lab = paste0(lactulose,":",n))

no_abx_rif.subs.str <- paste(no_abx_rif.subs$lab,collapse = "; ")

no_abx_rif.df %>% 
  ggplot(aes(lactulose, pctseqs)) +
  geom_boxplot(aes(fill = lactulose), width = 0.65, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(aes(color = lactulose), width = 0.25,alpha = 0.55) +
  paletteer::scale_color_paletteer_d("rockthemes::deelite") +
  paletteer::scale_fill_paletteer_d("rockthemes::deelite") +
  theme_bw() +
  stat_pvalue_manual(no_abx_rif.pdf) +
  facet_wrap(. ~ taxon, scales = "free", ncol = 3) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  labs( x= "", y = "Taxon Abund.",
        title = "Comparison of lactulose status in samples without antibiotics nor rifaximin exposure",
        caption = no_abx_rif.subs.str) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))

ggsave("results/3.B.LD117.no_abx_no_rifaximin.boxplot.pdf", 
       width = 7.75, height = 3.85)

# Fig 5, Bifido expansion vs. no expansion --------------------------------

bifido_pal = c("#5050FFFF","#CE3D32FF")

bact.df <- metab %>% 
  filter(cohort == "Liver\nDisease") %>% 
  select(seq_id, lactulose, Enterococcus, Proteobacteria, Bifidobacterium) %>% 
  mutate(bifido_cat = if_else(Bifidobacterium < 0.1,
                              "Less Bifido\n(<10%)",
                              "More Bifido\n(>=10%)")) %>% 
  select(-Bifidobacterium) %>% 
  gather("taxon","pctseqs",-c(bifido_cat, seq_id, lactulose)) 

bact.p <- bact.df %>% 
  group_by(lactulose, taxon) %>% 
  t_test(pctseqs ~ bifido_cat) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  add_y_position() %>% 
  add_significance() 

bact.df %>% 
  ggplot(aes(bifido_cat, pctseqs)) +
  geom_boxplot(aes(fill = bifido_cat), 
               width = 0.65, alpha = 0.55, outlier.shape = NA) +
  geom_jitter(aes(color = bifido_cat), width = 0.25,alpha = 0.55) +
  scale_fill_manual(values = bifido_pal) +
  scale_color_manual(values = bifido_pal) +
  theme_bw() +
  facet_wrap(. ~ lactulose + taxon, ncol = 4) +
  # stat_compare_means(label = "p.signif") +
  stat_pvalue_manual(bact.p, label = "p.adj.signif", tip.length = 0.015) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  labs( x= "", y = "Bacteria Abund.\n(log transformed)") +
  scale_y_log10(
    # labels = scales::scientific,
    breaks = c(0, 0.001, 0.01, 0.1, 1),
    expand = expansion(mult = c(0., 0.1))
  )

ggsave("results/5.A.lactulose.Proteo_Entero.Bifido.boxplot.pdf",
       height = 4.5, width = 9)

# Table S1 -----------------------------------------------------------------
# Exclude healthy donor from this analysis

## Bifidobacterium ---------------------------------------------------------
lm.bifido.ld.biLact <- lm(Bifidobacterium ~ 
                            lactulose + PPI + stool_consistency + abx,
                          data = metab %>% 
                            filter(cohort != "Healthy\nDonor") )

tidy(lm.bifido.ld.biLact) %>% 
  write_csv("./results/3.TableS1.lm.LD262.Bifido.csv")

## Enterococcus -----------------------------------------------------------

lm.entero.ld.biLact <- lm(Enterococcus ~ 
                            lactulose + PPI + stool_consistency + abx,
                          data = metab %>% 
                            filter(cohort != "Healthy\nDonor") )

tidy(lm.entero.ld.biLact) %>% 
  write_csv("./results/3.TableS1.lm.LD262.Entero.csv")
