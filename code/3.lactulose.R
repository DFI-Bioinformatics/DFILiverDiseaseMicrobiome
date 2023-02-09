library(tidyverse)
library(yingtools2)

# Figure 3: taxonomy barplot between lactulose vs. no lactulose

# load source code for plot -----------------------------------------------

source("./code/getRdpPal.R")

# load in data ------------------------------------------------------------

metab <- read_csv("./data/LD850.meta.quant.metabolomics.csv") 

mpa <- read_csv("./data/LD850.metaphlan.csv")

# calculate rel abund -----------------------------------------------------

bk <- mpa %>%
  dplyr::count(seq_id, taxid, Kingdom, Phylum, Class, Order, Family,
               Genus, wt = relative_abundance, name = "pctseqs") %>%
  add_count(seq_id, wt = pctseqs, name = "Total") %>%
  mutate(pctseqs = pctseqs/Total)  %>%
  mutate(genLab = Genus,
         Genus = paste0(Phylum,"-",Order,"-", Family, "-",Genus))

# clinical variables ------------------------------------------------------

clin.gg <- metab %>%
  select(ID, `systemic antibiotics`, 
         rifaximin, lactulose,
         bifido_percent = Bifidobacterium) %>%
  gather("var","value", -ID,-bifido_percent) %>%
  mutate(var = factor(var, levels = c("systemic antibiotics",
                                      "rifaximin","lactulose"),
                      labels = c("systemic antibiotics",
                                 "rifaximin","lactulose"))) %>%
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
        axis.ticks.x = element_blank()
  ) +
  facet_grid( . ~ lactulose, scales = "free", space = "free") +
  labs(x="", y= "", fill = "administered") +
  scale_fill_manual(values = c("no" = "#09ba38", "yes" = "#db3b3b"))

clin.gg

# toxins ------------------------------------------------------------------

sorted.toxingg <- metab %>% 
  select(seq_id, vanA, vanB) %>%
  gather("var","value", -seq_id) %>%
  replace_na(list(value = 0)) %>% # view
  left_join(metab %>% select(seq_id, ID, bifido_percent = Bifidobacterium, lactulose)) %>%
  ggplot(aes(reorder(ID,-bifido_percent), var)) +
  geom_tile(aes(fill = value)) +
  scale_fill_binned(trans = yingtools2::log_epsilon_trans(10),
                    type = "viridis",
                    # direction = -1
  ) +
  theme_bw() +
  facet_grid( . ~ lactulose, scales = "free", space = "free") +
  theme(legend.position = "right",
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8.5),
        #axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank()
  ) +
  labs(x="", y= "", fill = "RPKM")

sorted.toxingg

# taxonomy barplot sorted by bifido ---------------------------------------

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
              select(ID, seq_id, 
                     bifido_percent = Bifidobacterium,
                     lactulose_before_stool_final = lactulose)) %>%
  ggplot(aes(x=reorder(ID,-bifido_percent), y=pctseqs)) +
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
  facet_grid( . ~ lactulose_before_stool_final, scales = "free", space = "free") +
  labs(x="", y= "Bacterial Abund.")

sorted.taxgg

## combine all plots -------------------------------------------------------

pdf(file.path("./results/3.LD850.bifido_clinical_toxin.selected.sorted.lactulose.pdf"),
    height = 7.15, width = 14.5)
gg.stack(sorted.taxgg, clin.gg, sorted.toxingg,
         heights = c(1,0.28,0.19),
         newpage = F)
dev.off()
