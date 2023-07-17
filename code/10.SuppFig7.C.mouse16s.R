library(tidyverse)
library(phyloseq)
library(yingtools2)
source("R_Functions.R")
source("getRdpPal.R")


# Load in the phyceberg
phy <- readRDS("../data/mouse_consortia_phyloseq.rds")


#Table for 16s rRNA sequencing data
seq_table <- phyloseq_table(phy) %>% 
  type_convert() %>% 
  filter(pctseqs > 0) %>% 
  mutate(lactulose = factor(lactulose, levels = c("no lactulose", "lactulose"))) 

#Table for metabolomic quantification data
metabolomic_quant_table <- read_csv("../data/mouse_consortia_metab_quant.csv") %>% 
  mutate(lactulose = factor(lactulose, levels = c("no lactulose", "lactulose"))) %>% 
  mutate(sampleid = factor(sampleid, levels = str_sort(unique(sampleid), numeric = TRUE)))



#Generate color palette for stacked barplot
pal <- getRdpPal(seq_table)




gg1 <- seq_table %>%
  group_by(sample, sampleid, Kingdom, Phylum, Class, Order, Family, Genus,
           sample.type, mouse.number, day, lactulose, consortia) %>%
  summarise(pctseqs=sum(pctseqs)) %>%
  ungroup() %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus) %>% 
  mutate(Genus = factor(Genus, levels = unique(Genus))) %>% 
  group_by(sample, lactulose, consortia, day) %>% 
  arrange(Genus) %>% 
  mutate(cum.pct = cumsum(pctseqs), 
         y.text = (cum.pct + c(0, cum.pct[-length(cum.pct)]))/2) %>% 
  ungroup() %>%
  dplyr::select(-cum.pct) %>% 
  mutate(tax.label= if_else(grepl("unclassified",Genus), 
                            "unclassified",
                            as.character(Genus))) %>%
  mutate(tax.label = if_else(pctseqs >= .1, tax.label, "")) %>%
  ggplot(aes(x=sampleid,y=pctseqs)) +
  geom_bar(aes(fill=Genus),stat="identity") +
  geom_text(aes(y=1-y.text,x=sampleid,label=tax.label), angle=90,
            lineheight=0.6,size=2.5) +
  scale_fill_manual(values=pal) +
  facet_grid(. ~ consortia+lactulose, scales = "free",space = "free") +
  ylab("16S Relative Abundance") +
  xlab("") +
  theme_bw() +
  theme(strip.text.x=element_text(size=15),
        legend.position="none",
        axis.text.x=element_text(angle=90),
        axis.text.x.bottom = element_blank()) 


gg2 <- metabolomic_quant_table %>% 
  ggplot(aes(x=sampleid, y=Concentration_mM)) +
  geom_bar(aes(fill = Compound), stat="identity") +
  facet_grid(Compound ~ consortia+lactulose, scales = "free", space = "free_x") +
  ylab("Concentration (mM)") +
  xlab("") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12 ),
        strip.text.y = element_text(angle=0),
        strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        axis.text.x.bottom = element_blank()) +
  scale_fill_manual(values = c("#7F8BA2", 
                               "#F1757F"))



pdf(file = "../results/7S.C.mouse_16s.pdf", height = 25, width = 30)
ggstack <- gg.stack(gg1,gg2,heights=c(4,3), newpage = F)
dev.off()




