library(phyloseq)
library(readxl)
library(yingtools2)

# Supp Fig 6, D: 16S copy number of Bifidobacterium between lactulose- vs. water- treated GF mice.

# load data ---------------------------------------------------------------

qpcr <- readxl::read_xlsx("./data/qPCR_Calculations_MattO_220708 GF.xlsx", skip = 6) %>% 
  dplyr::rename(copy_number = `Original copies/ul`,
                sample = ID) %>% 
  mutate(copy_number = as.numeric(copy_number))

tempphy <- readRDS("./data/MMF.16S.263_MattOdenwald.finalPhy_rdp.rds")

t <- get.otu.melt(tempphy) %>%
  replace_na(list(Species="unclassified",
                  Genus="unclassified",
                  Family="unclassified",
                  Order="unclassified",
                  Class="unclassified",
                  Phylum="unclassified")) %>%
  mutate(Genus=ifelse(Genus=="unclassified",
                      paste(Family,Genus,sep="\n"),
                      as.character(Genus)))

# Supp Fig 6, D ---------------------------------------------------
  
t2 <- t %>% 
    mutate(day = str_extract(string = experiment.sample.name, pattern = "day.+"),
           day = str_to_title(string = day),
           day_number = str_extract(string = day, pattern = "[0-9]+"),
           day_number = as.numeric(day_number),
           treatment = gsub(pattern = "bifido", replacement = "", x = treatment, ignore.case = T),
           treatment = gsub(pattern = "\\s", replacement = "", x = treatment)) %>% 
    filter(grepl(pattern = "Bifido", x = Genus)) %>% 
    arrange(day_number) 
  
  
stats <- 
    t2 %>%
    mutate(sample = sampleid,
           sample = gsub(sample, pattern = "\\.", replacement = "_")) %>% 
    group_by(group, sample, Kingdom,Phylum,Class,Order,Family,Genus,day,mouse.number,treatment,day_number) %>%
    summarise(pctseqs=sum(pctseqs)) %>%
    ungroup() %>%
    arrange(Kingdom, Phylum, Class, Order, Family, Genus) %>% 
    mutate(Genus = factor(Genus, levels = unique(Genus))) %>% 
    group_by(group, sample, day, mouse.number, treatment,day_number) %>% 
    arrange(Genus) %>% 
    full_join(qpcr,
              by = "sample") %>%
    ungroup() %>%
    mutate(tax.label= ifelse(Genus=="unclassified",
                             paste(Family,Genus,sep="\n"),
                             as.character(Genus)),
           day = ifelse(is.na(day), "Day 0", day),
           day_number = ifelse(is.na(day_number), 0, day_number),
           mouse.number = ifelse(is.na(mouse.number), str_extract(string = sample, pattern = "40[0-9]"), mouse.number),
           pctseqs = ifelse(is.na(pctseqs), NA, pctseqs),
           treatment = case_when(is.na(treatment) & mouse.number %in% as.character(seq(401, 404, 1))~"Lactulose",
                                 is.na(treatment) & mouse.number %in% as.character(seq(405, 408, 1))~"Water",
                                 TRUE~treatment)) %>% 
    ungroup() %>% 
    mutate(tax.label = if_else(pctseqs >= .1, tax.label, ""),
           treatment = factor(treatment, levels = c("Water", "Lactulose")),
           pctseqs = 100*pctseqs) %>% 
    group_by(day) %>% 
    rstatix::wilcox_test(copy_number ~ treatment)
  
  
  t2 %>%
    mutate(sample = sampleid,
           sample = gsub(sample, pattern = "\\.", replacement = "_")) %>% 
    group_by(group, sample, Kingdom,Phylum,Class,Order,Family,Genus,day,mouse.number,treatment,day_number) %>%
    summarise(pctseqs=sum(pctseqs)) %>%
    ungroup() %>%
    arrange(Kingdom, Phylum, Class, Order, Family, Genus) %>% 
    mutate(Genus = factor(Genus, levels = unique(Genus))) %>% 
    group_by(group, sample, day, mouse.number, treatment,day_number) %>% 
    arrange(Genus) %>% 
    full_join(qpcr,
              by = "sample") %>%
    ungroup() %>%
    mutate(tax.label= ifelse(Genus=="unclassified",
                             paste(Family,Genus,sep="\n"),
                             as.character(Genus)),
           day = ifelse(is.na(day), "Day 0", day),
           day_number = ifelse(is.na(day_number), 0, day_number),
           mouse.number = ifelse(is.na(mouse.number), str_extract(string = sample, pattern = "40[0-9]"), mouse.number),
           pctseqs = ifelse(is.na(pctseqs), NA, pctseqs),
           treatment = case_when(is.na(treatment) & mouse.number %in% as.character(seq(401, 404, 1))~"Lactulose",
                                 is.na(treatment) & mouse.number %in% as.character(seq(405, 408, 1))~"Water",
                                 TRUE~treatment)) %>% 
    ungroup() %>% 
    mutate(tax.label = if_else(pctseqs >= .1, tax.label, ""),
           treatment = factor(treatment, levels = c("Water", "Lactulose")),
           pctseqs = 100*pctseqs) %>%
    ggplot(aes(x = reorder(day, day_number), y = copy_number , color = mouse.number, size = pctseqs, group = mouse.number)) +
    geom_point(alpha = 0.35) +
    geom_line(size = 0.5) +
    facet_wrap(~treatment) +
    ylab("16S Copy Number\n") +
    xlab("") +
    theme_bw() +
    theme(
      axis.text.x=element_text(angle=60, vjust = 1, hjust = 1, color = "black", size = 12),
      axis.text.y=element_text(color = "black", size = 12),
      axis.title.y=element_text(color = "black", size = 12),
      strip.text.x=element_text(size=13),
      panel.grid = element_blank(),
      legend.position="right",
      legend.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("#8B0000", "#FFD700", "#EE7621", "#FF4040", "#0000CD", "#1E90FF", "#2E8B57", "#43CD80"))+
    scale_size_continuous("Bifidobacterium\nRelative Abundance (%)",
                          limits = c(0,100)) +
    scale_y_log10(
      breaks=c(1,10,100,1000,10000,100000,1000000,10000000,10000000,10000000,100000000),
      labels=scales::scientific(c(1,10,100,1000,10000,100000,1000000,10000000,10000000,10000000,100000000)),
      limits = c(1,100000000)) +
    guides(color = guide_legend("Mouse", ncol = 2))
  
ggsave("results/6S.D.bifido_16s.pdf", height = 6, width = 8)  
