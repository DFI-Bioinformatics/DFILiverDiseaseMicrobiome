library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
library(paletteer)

# Figure 2: UMAP of LD847 samples + various metabolite overlay
# Figure S2: SCFA level comparison between dominant taxon groups
# FIgure S3: Bile acid level comparison between dominant taxon groups

# load in data ------------------------------------------------------------

metab <- read_csv("./data/LD847_HD22.quant.meta.csv") %>% 
  mutate(diversity = factor(diversity, 
                            levels = c("Low", "Medium", "High", "Healthy\nDonor")))

mpa <- read_csv("./data/LD847_HD22.metaphlan.csv")

taxumapcoord <- read_csv("./data/taxumap_embedding.LD847.n375.csv") %>% 
  rename(seq_id = index_column)

n.neigh = 375
key.fn = "2.B.LD847"

# A. overlay dominant taxon to umap ------------------------------------------

bk <- mpa %>%
  add_count(seq_id, wt = relative_abundance, name = "totalAbd") %>%
  mutate(pctseqs = relative_abundance/totalAbd)  %>%
  dplyr::count(seq_id, taxid, Kingdom, Phylum, Class, Order, Family,
               Genus, wt = pctseqs, name = "pctseqs") %>%
  mutate(genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus))

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
  dplyr::select(-cum.pct) %>%
  mutate(#tax.label=ifelse(Species=="unclassified",paste(Family,Genus,Species,sep="\n"),paste(Genus,Species,sep="\n")),
    tax.label= ifelse(grepl("unclassified",genLab), Family,genLab),
    tax.label = ifelse(pctseqs >= .1, as.character(tax.label), ""))

## find the most abundant taxon ----------------------------------------
domi <- ggdf %>%
  group_by(seq_id) %>%
  slice_max(pctseqs, n = 1, with_ties = F) %>%
  left_join(metab %>% select(ID, seq_id)) %>%
  select(-c(y.text,tax.label)) %>%
  mutate(Genus = as.character(Genus),
         Genus = if_else(pctseqs < 0.05, "Other", Genus))

g2tax <- domi %>%
  ungroup() %>%
  count(Phylum, Order, Family, Genus, genLab, sort = T) %>%
  mutate(taxgrp = case_when(
    genLab %in% c("Enterococcus","Streptococcus","Staphylococcus", "Bifidobacterium") ~ genLab,
    Family %in% c("Lachnospiraceae","Ruminococcaceae","Oscillospiraceae","Erysipelotrichaceae",
                  "Lactobacillaceae") ~ Family,
    Order %in% c("Clostridiales") ~ Order,
    Phylum %in% c("Proteobacteria","Actinobacteria","Bacteroidetes") ~ Phylum,
    TRUE ~ "Others"
  )) %>%
  count(Genus, taxgrp, wt = n)

## assign colors ------------------------------------------------------------

palnames <- g2tax %>%
  count(taxgrp, wt = n, sort = T) %>%
  arrange(taxgrp) %>%
  pull(taxgrp)

palcols <- c("#B181A3","#51AB9B","#926184","#129246",
             "orange","#EC9B96","#3b51a3","#9AAE73",
             "grey","red","#f1eb25","#9FB846")

names(palcols) <- palnames

palcols <- palcols[c(1:8,10:12,9)]

ch.palcols <- palcols

ch.palcols["Actinobacteria"] <- "grey"

g2tax

tax.gg <- taxumapcoord %>%
  left_join(domi) %>%
  left_join(g2tax) %>% 
  ggplot(aes(taxumap1, taxumap2)) +
  geom_point(aes(color = taxgrp),alpha = 0.45, size = 5.5) +
  scale_color_manual(values = ch.palcols) +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = "", y = "",
       title = paste0("taxUMAP based on shotgun taxonomy"),
       caption = paste0("n_neighbors=",n.neigh),
       subtitle = paste0("colored by ", "top taxon"),
       color = "Dominant\nTaxon") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold.italic", size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 12))

tax.gg

ggsave("./results/2.A.LD847.umap.domination.recolor.pdf",
       width = 8.95, height = 7.55)


# B. overlay metabolite conc to umap -----------------------------------------

my_comps <- list(c("Bifidobacterium","Lactobacillaceae"),
                 c("Bifidobacterium","Bacteroidetes"),
                 c("Bifidobacterium","Lachnospiraceae"),
                 c("Bifidobacterium","Enterococcus"),
                 c("Bifidobacterium","Proteobacteria"))

## metabolomics long format ------------------------------------------------

metablong <- metab %>%
  filter(cohort == "Liver\nDisease") %>% 
  select(ID,seq_id, `3-oxolithocholic acid`:succinate) %>%
  gather("compound","value", -seq_id, -ID)

metablong %>%
  count(compound)

## Fig S2, scfa ---------------------------------------------------

abreaks = c(10, 20, 50, 100, 150)
bbreaks = c(5, 10, 20, 40, 60)
pbreaks = c(2.5, 5, 10, 45,  90)

allbreaks = list("acetate" = abreaks,
                 "propionate" = pbreaks,
                 "butyrate" = bbreaks)

for (scfa in c("acetate","propionate","butyrate")){
  
  tempmetab <- metablong %>%
    filter(compound %in% scfa) %>%
    left_join(domi) %>%
    left_join(g2tax) %>%
    left_join(taxumapcoord) %>%
    mutate(taxgrp = fct_infreq(taxgrp),
           taxgrp = fct_relevel(taxgrp, "Others", after = Inf))
  
  clustgg <- tempmetab %>%
    ggplot(aes(taxgrp, value, color = taxgrp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.55) +
    scale_color_manual(values = palcols) +
    ggpubr::stat_compare_means(comparisons = my_comps,
                               p.adjust.methods = "BH") +
    theme_bw() +
    theme(legend.position = "none",
          strip.text.x = element_text(face = "bold", size = 12),
          strip.background.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
    labs(x = "Dominant Taxon", y = "Metabolite Conc. (mM)") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))
  
  clustgg
  
  scfagg <- tempmetab %>%
    ggplot(aes(taxumap1, taxumap2, color = value)) +
    geom_point(alpha = 0.35, size = 2.5) +
    theme_bw()+
    scale_color_binned(direction = -1, type = "viridis",
                       n.breaks = 5,
                       breaks = allbreaks[[scfa]]) +
    labs(x = "taxUmap1", y = "taxUmap2", color = "")
  
  scfagg +
    labs(title = scfa) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold.italic", size = 10),
          strip.background = element_blank(),
          axis.text = element_text(size = 6.5),
          axis.title = element_text(size = 7.5),
          legend.key.size = unit(0.75, "lines"),
          legend.text = element_text(size = 6),
          plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 8.5)) 
  
  ggsave(file.path("results",paste0(key.fn,".",gsub(" ","_",scfa),".umap.lessCat.pdf")),
         width = 4.25, height = 3.35)
  
  cplt <- cowplot::plot_grid(scfagg, clustgg, rel_widths = c(1.15, 1.25), align = "h")
  
  title <- cowplot::ggdraw() +
    cowplot::draw_label(paste0(scfa), fontface='bold')
  
  cowplot::plot_grid(title, cplt, ncol=1, rel_heights=c(0.1, 1))
  
  ggsave(file.path("results",paste0(key.fn,".",gsub(" ","_",scfa),".combined.toptax.lessCat.pdf")),
         width = 9.5, height = 4.35)
}

## Fig S3, bile acids -----------------------------------------------------

metablong %>%
  filter(! compound %in% c("acetate","propionate","butyrate","succinate")) %>%
  left_join(domi) %>%
  left_join(g2tax) %>%
  left_join(taxumapcoord) %>%
  mutate(taxgrp = fct_infreq(taxgrp),
         taxgrp = fct_relevel(taxgrp, "Others", after = Inf),
         value = if_else(value <0, 0, value)) %>%
  ggplot(aes(taxgrp, value, color = taxgrp)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.55) +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.1),
                     labels = function(x) format(x, scientific = T)
  ) +
  facet_wrap("compound", scales = "free", ncol = 4)+
  scale_color_manual(values = palcols) +
  ggpubr::stat_compare_means(comparisons = my_comps, p.adjust.methods = "BH") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(face = "bold", size = 12),
        strip.background.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(x = "Dominant Taxon", y = "Metabolite Conc. (ug/ml)")

ggsave(file.path("results", paste0(key.fn,".toptaxon.bile.pdf")),
       width = 17.5, height = 9.35)

### bile acid umap overlay ------------------------------------------------------

bilemetab <- c("cholic acid",
               "glycocholic acid","taurocholic acid",
               "lithocholic acid","deoxycholic acid",
               "alloisolithocholic acid", "3-oxolithocholic acid")

for (scfa in bilemetab){
  
  tempmetab <- metablong %>%
    filter(compound %in% scfa) %>%
    left_join(domi) %>%
    left_join(g2tax) %>%
    left_join(taxumapcoord) %>%
    mutate(taxgrp = fct_infreq(taxgrp),
           taxgrp = fct_relevel(taxgrp, "Others", after = Inf))
  
  clustgg <- tempmetab %>%
    mutate(value = if_else(value < 0, 0, value)) %>% 
    ggplot(aes(taxgrp, value, color = taxgrp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.55) +
    scale_color_manual(values = palcols) +
    scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.1),
                       expand = expansion(mult = c(0,0.1)),
                       labels = function(x) format(x, scientific = TRUE)) +
    ggpubr::stat_compare_means(comparisons = my_comps,
                               p.adjust.methods = "BH") +
    theme_bw() +
    theme(legend.position = "none",
          strip.text.x = element_text(face = "bold", size = 12),
          strip.background.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
    labs(x = "Dominant Taxon", y = "Metabolite Conc. (ug/mL)")
  
  clustgg
  
  scfagg <- tempmetab %>%
    ggplot(aes(taxumap1, taxumap2, color = value)) +
    geom_point(alpha = 0.35, size = 2.5) +
    theme_bw()+
    scale_color_binned(type = "viridis",
                       direction = -1,
                       trans = yingtools2::log_epsilon_trans(0.1)) +
    labs(x = "taxUmap1", y = "taxUmap2", color = "")
  
  scfagg +
    labs(title = scfa) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold.italic", size = 10),
          strip.background = element_blank(),
          axis.text = element_text(size = 6.5),
          axis.title = element_text(size = 7.5),
          legend.key.size = unit(0.75, "lines"),
          legend.text = element_text(size = 6),
          plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 8.5)) 
  
  ggsave(file.path("results",paste0(key.fn,".",gsub(" ","_",scfa),".umap.lessCat.pdf")),
         width = 4.25, height = 3.35)
  
  cplt <- cowplot::plot_grid(scfagg, clustgg, rel_widths = c(1.25, 1.30), align = "h")
  
  title <- cowplot::ggdraw() +
    cowplot::draw_label(paste0(scfa), fontface='bold')
  
  cowplot::plot_grid(title, cplt, ncol=1, rel_heights=c(0.1, 1))
  
  ggsave(file.path("results",paste0(key.fn,".",gsub(" ","_",scfa),".combined.toptax.pdf")),
         width = 9.5, height = 4.35)
}

