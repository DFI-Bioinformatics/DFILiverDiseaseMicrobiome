library(ggpmisc)

# Figure 2: TaxUmap on 850 samples + statistical analysis on top taxon

## TaxUmap is performed on the commendline, the resultant coordinates are saved to use
## Github repo: https://github.com/jsevo/taxumap

# load source code for plot -----------------------------------------------

source("./code/getRdpPal.R")

# load data ---------------------------------------------------------------

taxumapcoord <- read_csv(("./data/embedding.LD850.csv")) %>%
  rename(seq_id = index_column)

metab <- read_csv("./data/LD850.meta.quant.metabolomics.csv") 

mpa <- read_csv("./data/LD850.metaphlan.csv")

# transform data into plot-able dataframe ----------------------------------

bk <- mpa %>%
  add_count(seq_id, wt = relative_abundance, name = "totalAbd") %>%
  mutate(pctseqs = relative_abundance/totalAbd)  %>%
  dplyr::count(seq_id, taxid, Kingdom, Phylum, Class, Order, Family,
               Genus, wt = pctseqs, name = "pctseqs") %>%
  add_count(seq_id, wt = pctseqs, name = "Total") %>%
  mutate(pctseqs = pctseqs/Total)  %>%
  mutate(genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus))

pal <- getRdpPal(bk)
pal <- c(pal, "Other" = "#99adcf")

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

# determine top taxon -----------------------------------------------------
## if the most abundant taxon is less than 0.05, we call it "Other"
domi <- ggdf %>%
  group_by(seq_id) %>%
  slice_max(pctseqs, n = 1, with_ties = F) %>%
  left_join(metab %>% select(ID, seq_id)) %>%
  select(-c(y.text,tax.label)) %>%
  mutate(Genus = as.character(Genus),
         Genus = if_else(pctseqs < 0.05, "Other", Genus))

# plot scatter plot of taxumap coordinates --------------------------------

taxumapcoord %>%
  left_join(domi) %>%
  ggplot(aes(taxumap2, taxumap1, color = Genus)) +
  geom_point(alpha = 0.55, size = 3.5) +
  scale_color_manual(values = pal, guide = "none") +
  theme_bw() +
  # theme(legend.position = "none") +
  labs(x = "", y = "",
       title = paste0("taxUMAP based on shotgun taxonomy"),
       caption = paste0("n_neighbors=",350),
       subtitle = paste0("colored by ", "top taxon") )

ggsave("./results/2.LD850.umap.domination.pdf",
       width = 9.5, height = 8.2)

# overlay metabolites conc. and statistical tests -------------------------

## top taxon as group ------------------------------------------------------------------

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

## define color for top taxon groups ---------------------------------------

palnames <- g2tax %>%
  count(taxgrp, wt = n, sort = T) %>%
  arrange(taxgrp) %>%
  pull(taxgrp)

palcols <- c("#B181A3","#51AB9B","#926184","#129246",
                      "orange","#EC9B96","#3b51a3","#9AAE73",
                      "grey","red","#f1eb25","#9FB846")
                      
names(palcols) <- palnames

palcols <- palcols[c(1:8,10:12,9)]

g2tax

## define comparison groups ------------------------------------------

my_comps <- list(c("Bifidobacterium","Lactobacillaceae"),
                 c("Bifidobacterium","Bacteroidetes"),
                 c("Bifidobacterium","Lachnospiraceae"),
                 c("Bifidobacterium","Enterococcus"),
                 c("Bifidobacterium","Proteobacteria"))

# transform metabolomics data into long form ------------------------------

metablong <- metab %>%
  select(ID,seq_id, `3-oxolithocholic acid`:succinate) %>%
  gather("compound","value", -seq_id, -ID)

# scfa --------------------------------------------------------------------

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
    labs(x = "Dominant Taxon", y = "Metabolite Conc. (mM)")
  
  clustgg
  
  scfagg <- tempmetab %>%
    ggplot(aes(taxumap2, taxumap1, color = value)) +
    geom_point(alpha = 0.45, size = 2.5) +
    theme_bw()+
    scale_color_binned(direction = -1, type = "viridis",
                       n.breaks = 5,
                       breaks = allbreaks[[scfa]]) +
    labs(x = "taxUmap1", y = "taxUmap2", color = "")
  
  scfagg
  
  cplt <- cowplot::plot_grid(scfagg, clustgg, rel_widths = c(1.15, 1.25), align = "h")
  
  title <- cowplot::ggdraw() +
    cowplot::draw_label(paste0(scfa), fontface='bold')
  
  cowplot::plot_grid(title, cplt, ncol=1, rel_heights=c(0.1, 1))
  
  ggsave(paste0("./results/","2.",gsub(" ","_",scfa),".combined.toptax.lessCat.pdf"),
         width = 9.5, height = 4.35)
}

# bile acids --------------------------------------------------------------

metablong %>%
  filter(! compound %in% c("acetate","propionate","butyrate","succinate")) %>%
  left_join(domi) %>%
  left_join(g2tax) %>%
  left_join(taxumapcoord) %>%
  mutate(taxgrp = fct_infreq(taxgrp),
         taxgrp = fct_relevel(taxgrp, "Others", after = Inf)) %>%
  ggplot(aes(taxgrp, value, color = taxgrp)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.55) +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.1),
                     # labels = scales::math_format(10^.log10(x))
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

ggsave(paste0("./results/2.","LD850.toptaxon.bile.pdf"),
       width = 17.5, height = 9.35)

# bile acid umap overlay ------------------------------------------------------

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
    ggplot(aes(taxgrp, value, color = taxgrp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.55) +
    scale_color_manual(values = palcols) +
    scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.1)) +
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
    ggplot(aes(taxumap2, taxumap1, color = value)) +
    geom_point(alpha = 0.45, size = 2.5) +
    theme_bw()+
    scale_color_binned(type = "viridis",
                       direction = -1,
                       trans = yingtools2::log_epsilon_trans(0.1)) +
    labs(x = "taxUmap1", y = "taxUmap2", color = "")
  
  scfagg
  
  cplt <- cowplot::plot_grid(scfagg, clustgg, rel_widths = c(1.25, 1.30), align = "h")
  
  title <- cowplot::ggdraw() +
    cowplot::draw_label(paste0(scfa), fontface='bold')
  
  cowplot::plot_grid(title, cplt, ncol=1, rel_heights=c(0.1, 1))
  
  ggsave(paste0("./results/2.",gsub(" ","_",scfa),".combined.toptax.pdf"),
         width = 9.5, height = 4.35)
}


# linear regression between bifido vs. scfa -------------------------------
## in more Bifido group

metablong %>%
  filter(compound %in% c("acetate","propionate","butyrate")) %>% 
  left_join(metab %>% 
              select(seq_id, Bifidobacterium)) %>% 
  filter(Bifidobacterium > 0.1) %>% 
  ggplot(aes(Bifidobacterium, value)) +
  geom_point(alpha = 0.45, color = "darkblue") +
  facet_grid(compound ~ ., scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(face = "bold", size = 10.5),
        strip.background = element_blank()) +
  geom_smooth( color = "red", method = "lm", linetype = 2) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label),
                                 after_stat(p.value.label),
                                 sep = "*\", \"*"))) +
  labs( x= "Bifidobacterium Abundance",
        y = "Metabolite Conc. (mM)")

ggsave(paste0("./results/2.","scfa.Bifido.more.pdf"),
       height = 9.5, width = 4.35)
  
