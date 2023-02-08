library(tidyverse)
library(vegan)
library(ggdendro)
library(cowplot)

# Figure 1: Hcluster + metabolites

# load source code for plot -----------------------------------------------

source("./code/dendro_source.R")
source("./code/getRdpPal.R")

# load in data ------------------------------------------------------------

metab <- read_csv("./data/LD850.meta.quant.metabolomics.csv") 

mpa <- read_csv("./data/LD850.metaphlan.csv")

# form matrix ---------------------------------------------------------

kmat <- mpa %>%
  select(seq_id, taxid,
            relative_abundance) %>%
  pivot_wider(names_from =taxid,
              values_from = relative_abundance,
              values_fn = sum,
              values_fill = 0) %>%
  column_to_rownames(var = "seq_id")

# hcluster ----------------------------------------------------------------

## transformation before Hcluster ---------------------------------------
abdmat <- kmat/rowSums(kmat)*100

bc <- vegdist(abdmat, "bray")

cl <- hclust(bc, "ward.D2")

## cut dendrogram into 8 clusters ------------------------------------------

p1_dendro = dendro_data_k(cl, k = 8)

p1 <- plot_ggdendro(p1_dendro,
                    direction   = "tb",
                    label.size  = 0,
                    branch.size = 0.5,
                    nudge.label = 0,
                    expand.y    = 0.2)

ggtr <- p1 + theme_void()  +
  coord_cartesian(xlim = c(0, nrow(metab) + 1),
                  ylim = c(-0.05, max(p1_dendro$segments$y) + 0.3), expand = F) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# tax barplot dataframe ---------------------------------------------------

bk <- mpa %>%
  add_count(seq_id, wt = relative_abundance, name = "totalAbd") %>%
  mutate(pctseqs = relative_abundance/totalAbd)  %>%
  dplyr::count(seq_id, taxid, Kingdom, Phylum, Class, Order, Family,
               Genus, wt = pctseqs, name = "pctseqs") %>%
  add_count(seq_id, wt = pctseqs, name = "Total") %>%
  mutate(pctseqs = pctseqs/Total)  %>%
  mutate(genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus))

## generate customized color palette ---------------------------------------

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
  left_join(p1_dendro$labels, by = c("seq_id" = "label")) %>%
  ggplot(aes(x=reorder(seq_id,x), y=pctseqs)) +
  geom_bar(aes(fill=Genus),stat="identity") +
  # geom_text(aes(y=1-y.text,x=ID,label=tax.label),angle=90,
  #           lineheight=0.6,size=2.2) +
  scale_fill_manual(values=pal) +
  theme_bw() +
  coord_cartesian(xlim = c(0, nrow(p1_dendro$labels) + 1), expand = F) +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5),
        axis.text.x = element_blank()) +
  labs(x="", y= "Bacterial Abund.")

# metabolites -------------------------------------------------------------

metabdf <- metab %>%
  select(seq_id, ID, `3-oxolithocholic acid`:succinate) %>%
  gather("compound", "value", -seq_id, -ID)

## scfa panel --------------------------------------------------------------

scfametab <- c("acetate","butyrate","propionate")

scfagg <- metabdf %>%
  filter(compound %in% scfametab) %>%
  mutate(compound = factor(compound, levels = scfametab),
         biletype=case_when(
           compound %in% c("invSimp") ~ "Inverse Simpson",
           compound %in% c("cholic acid") ~ "primary",
           compound %in% c("lithocholic acid","deoxycholic acid","isodeoxycholic acid") ~ "secondary",
           compound %in% c("glycocholic acid","taurocholic acid") ~ "conjugated",
           compound %in% c("propionate","acetate","butyrate") ~ "short chain fatty acid",
           compound %in% c("succinate") ~ "succinate",
           TRUE ~ "tertiary")) %>%
  left_join(p1_dendro$labels, by = c("seq_id" = "label")) %>%
  ggplot(aes(x=reorder(seq_id,x), y=value)) +
  geom_col(aes(fill = biletype)) +
  scale_fill_manual(values = "#00468BFF") +
  facet_grid(compound ~., scale = "free") +
  theme_bw()+
  coord_cartesian(xlim = c(0, nrow(p1_dendro$labels) + 1), expand = F) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5.5),
        strip.text.y = element_text(angle = 0, face = "bold", hjust = 0),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.02, "lines"),
        strip.background.y = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "", y = "Short Chain Fatty Acid Conc.")

## bile acid panel ---------------------------------------------------------

bilemetab <- c("cholic acid",
               "glycocholic acid","taurocholic acid",
               "lithocholic acid","deoxycholic acid",
               "alloisolithocholic acid", "3-oxolithocholic acid")

bikegg <-  metabdf %>%
  filter(compound %in% bilemetab) %>%
  mutate(compound = factor(compound, levels = bilemetab),
         biletype=case_when(
           compound %in% c("invSimp") ~ "Inverse Simpson",
           compound %in% c("cholic acid") ~ "primary",
           compound %in% c("lithocholic acid","deoxycholic acid","isodeoxycholic acid") ~ "secondary",
           compound %in% c("glycocholic acid","taurocholic acid") ~ "conjugated",
           compound %in% c("propionate","acetate","butyrate") ~ "short chain fatty acid",
           compound %in% c("succinate") ~ "succinate",
           TRUE ~ "tertiary")) %>%
  left_join(p1_dendro$labels, by = c("seq_id" = "label")) %>%
  ggplot(aes(x=reorder(seq_id,x), y=value)) +
  geom_col(aes(fill = biletype)) +
  facet_grid(compound ~., scale = "free") +
  theme_bw()+
  scale_fill_manual(values = c("#ED0000FF", "#42B540FF", "#0099B4FF", "#FDAF91FF")) +
  scale_y_continuous(trans = yingtools2::log_epsilon_trans(0.1)) +
  coord_cartesian(xlim = c(0, nrow(p1_dendro$labels) + 1), expand = F) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5.5),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold", hjust = 0),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.02, "lines"),
        strip.background.y = element_blank()) +
  labs(x = "", y = "Bile Acid Conc.\n(log10 transformed)")

# combine all plots -------------------------------------------------------

pdf(file.path("./results/1.LD850.hclustAbd.sepmetab.metaphlan.pdf"), 
    height = 12.5, width = 13.5)
plot_grid(ggtr, taxgg, scfagg, bikegg,
          ncol = 1,
          align = "hv", axis = "tblr",
          rel_heights = c(1,1.7, 0.85,1.95))
dev.off()
