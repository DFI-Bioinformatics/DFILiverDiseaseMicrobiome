
#######To clean up and output tibble from phyloseq input

phyloseq_table <- function (phy, meta_labels = NULL) {
  
  #Extract phyloseq file data
  sample_table <- data.frame(sample_data(phy))
  sample_table <- sample_table %>% 
    mutate(sample = rownames(sample_table))
  sample_table <- sample_table[, colSums(is.na(sample_table)) != nrow(sample_table)] #remove columns with no data
  
  tax_table <- data.frame(tax_table(phy))
  tax_table <- tax_table %>% 
    mutate(taxonID = rownames(tax_table)) %>%
    replace_na(list(Species="unclassified",
                    Genus="unclassified",
                    Family="unclassified",
                    Order="unclassified",
                    Class="unclassified",
                    Phylum="unclassified")) %>%
    mutate(Genus=ifelse(Genus=="unclassified",
                        paste(Family,Genus,sep="\n"),
                        as.character(Genus))) #Change na in taxonomy to unclassified
  
  otu_table <- as.data.frame.table(otu_table(phy))
  otu_table <- otu_table %>% 
    dplyr::rename(sample = Var1, taxonID = Var2, numseqs = Freq) %>% 
    group_by(sample) %>% 
    mutate(nseqs = sum(numseqs)) %>% 
    ungroup() %>% 
    mutate(pctseqs = numseqs / nseqs)
  
  refseq_table <- data.frame(refseq(phy, FALSE))
  refseq_table <- refseq_table %>% 
    mutate(taxonID = rownames(refseq_table)) %>% 
    dplyr::rename(seq = refseq.phy..FALSE.)
  
  #Making the compiled table
  if (is.null(meta_labels) == TRUE) {
    t <- sample_table %>% 
      right_join(otu_table, multiple = "all") %>% 
      left_join(tax_table) %>% 
      left_join(refseq_table)
  } else {
    t <- sample_table %>% 
      right_join(otu_table, multiple = "all") %>% 
      left_join(tax_table) %>% 
      left_join(refseq_table) %>% 
      left_join(meta_labels)
  }
  
  t <- tibble(t)
  
return(t)  
}



################### Color Palette Stuff

get.yt.rdp <- function (tax) 
{
  if (class(tax)[1] %in% c("phyloseq", "taxonomyTable")) {
    tax <- get.tax(tax)
  }
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
             "Genus")
  if (!all(ranks %in% names(tax))) {
    stop("YTError: need to have taxon levels: Kingdom, Phylum, Class, Order, Family, Genus, Species")
  }
  tax.dict <- tax[, ranks] %>% distinct()
  tax.dict$color <- rep(shades("gray", variation = 0.25), 
                        length.out = nrow(tax.dict))
  proteo <- tax.dict$Phylum == "Proteobacteria"
  tax.dict$color[proteo] <- rep(shades("red", variation = 0.4), 
                                length.out = sum(proteo))
  actino <- tax.dict$Phylum == "Actinobacteria"
  tax.dict$color[actino] <- rep(shades("#A77097", variation = 0.25), 
                                length.out = sum(actino))
  bacteroidetes <- tax.dict$Phylum == "Bacteroidetes"
  tax.dict$color[bacteroidetes] <- rep(shades("#51AB9B", variation = 0.25), 
                                       length.out = sum(bacteroidetes))
  clost <- tax.dict$Order == "Clostridiales"
  tax.dict$color[clost] <- rep(shades("#9C854E", variation = 0.25), 
                               length.out = sum(clost))
  lachno <- tax.dict$Family == "Lachnospiraceae"
  tax.dict$color[lachno] <- rep(shades("#EC9B96", variation = 0.25), 
                                length.out = sum(lachno))
  rumino <- tax.dict$Family == "Ruminococcaceae"
  tax.dict$color[rumino] <- rep(shades("#9AAE73", variation = 0.25), 
                                length.out = sum(rumino))
  erys <- tax.dict$Family == "Erysipelotrichaceae"
  tax.dict$color[erys] <- rep(shades("orange", variation = 0.25), 
                              length.out = sum(erys))
  cid.colors.new <- c(Enterococcus = "#129246", Streptococcus = "#9FB846", 
                      Staphylococcus = "#f1eb25", Lactobacillus = "#3b51a3")
  cid <- cid.colors.new[match(tax.dict$Genus, names(cid.colors.new))]
  tax.dict$color <- ifelse(is.na(cid), tax.dict$color, cid)
  tax.palette <- structure(tax.dict$color, names = as.character(tax.dict$Genus))
  tax.palette
}

######
extract_legend<-function(ggplot_object){
  tmp <- ggplot_gtable(ggplot_build(ggplot_object))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}



######### Join two phyloseq files

join_physeq <- function (phy1, phy2) {
  
  #get the sequences and taxonIDs of the first phy file
  phy1_ASVs <- phyloseq_table(phy1) %>% 
    select(taxonID, seq) %>% 
    dplyr::rename(refseq.phy2..FALSE. = seq) %>% 
    distinct() 
  
  #match ASVs and fill in the NAs of phy1 taxonID
  matched_ASVs <- data.frame(refseq(phy2, FALSE))
  matched_ASVs <- matched_ASVs %>% 
    mutate(Var2 = rownames(matched_ASVs)) %>% 
    left_join(phy1_ASVs) %>% 
    mutate(taxonID = ifelse(is.na(taxonID), Var2, taxonID)) %>% 
    select(-refseq.phy2..FALSE.)
  
  
  #Getting the Sample Data
  sample_table <- sample_data(phy2)
  
  
  #Replace tax_table ASVs of phy2 with phy1
  tax_table <- data.frame(tax_table(phy2))
  tax_table <- tax_table %>% 
    mutate(Var2 = rownames(tax_table)) %>% 
    left_join(matched_ASVs)
  rownames(tax_table) <- tax_table$taxonID
  tax_table <- tax_table %>% 
    select(-Var2, -taxonID)
  tax_table <- as.matrix(tax_table)
  tax_table <- tax_table(tax_table)
  
  #Replace otu_table ASVs of phy2 with phy1
  otu_table <- as.data.frame.table(otu_table(phy2))
  otu_table <- otu_table %>% 
    left_join(matched_ASVs) %>% 
    mutate(Var2 = taxonID) %>% 
    select(-taxonID) %>% 
    pivot_wider(names_from = Var1, values_from = Freq)
  otu_table <- as.data.frame(otu_table)
  rownames(otu_table) <- otu_table$Var2
  otu_table <- otu_table %>% 
    select(-Var2)
  otu_table <- as.matrix(t(otu_table))
  otu_table <- otu_table(otu_table, taxa_are_rows = FALSE) 
  
  #Replace refseq_table ASVs of phy2 with phy1
  refseq_table <- data.frame(refseq(phy2, FALSE))
  refseq_table <- refseq_table %>% 
    mutate(Var2 = rownames(refseq_table)) %>% 
    left_join(matched_ASVs)
  refseq <- DNAStringSet(refseq_table$refseq.phy2..FALSE., use.names = TRUE)
  names(refseq) <- refseq_table$taxonID
  
  
  #Bringing it all together
  new_phy2 <- phyloseq(sample_table, tax_table, otu_table, refseq)
  phy <- merge_phyloseq(phy1, new_phy2)
  
  return(phy)  
}






#######Remove all of the ; in ASV names for phylogenetic Newick Files

fix_ASV_names <- function (phy) {
  
  #Get the sample_table from phy ASV
  sample_table <- sample_data(phy)
  
  #Clean the ; in tax_table ASV
  tax_table <- data.frame(tax_table(phy))
  tax_table <- tax_table %>% 
    mutate(taxonID = rownames(tax_table)) %>%
    replace_na(list(Species="unclassified",
                    Genus="unclassified",
                    Family="unclassified",
                    Order="unclassified",
                    Class="unclassified",
                    Phylum="unclassified")) %>%
    mutate(Genus=ifelse(Genus=="unclassified",
                        paste(Family,Genus,sep="\n"),
                        as.character(Genus))) %>% #Change na in taxonomy to unclassified
    mutate(taxonID = str_replace_all(taxonID, ";", "_"))
  rownames(tax_table) <- tax_table$taxonID
  tax_table <- tax_table %>% 
    select(-taxonID)
  tax_table <- as.matrix(tax_table)
  tax_table <- tax_table(tax_table)
  
  
  #Clean the ; in otu_table ASV
  otu_table <- as.data.frame.table(otu_table(phy))
  otu_table <- otu_table %>% 
    mutate(Var2 = str_replace_all(Var2, ";", "_")) %>% 
    pivot_wider(names_from = Var1, values_from = Freq)
  otu_table <- as.data.frame(otu_table)
  rownames(otu_table) <- otu_table$Var2
  otu_table <- otu_table %>% 
    select(-Var2)
  otu_table <- as.matrix(t(otu_table))
  otu_table <- otu_table(otu_table, taxa_are_rows = FALSE)  
  
  #Clean the ; in refseq ASV
  refseq_table <- data.frame(refseq(phy, FALSE))
  refseq_table <- refseq_table %>% 
    mutate(Var2 = rownames(refseq_table)) %>% 
    mutate(Var2 = str_replace_all(Var2, ";", "_"))
  refseq <- DNAStringSet(refseq_table$refseq.phy..FALSE., use.names = TRUE)
  names(refseq) <- refseq_table$Var2
  
  #Bringing it all together
  new_phy <- phyloseq(sample_table, tax_table, otu_table, refseq)
  
  
  return(new_phy) 
}


