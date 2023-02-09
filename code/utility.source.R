last2tax <- function(x) {
  sstr <- str_split(x,"\\|")[[1]]
  sstr <- gsub("_+$"," sp.",sstr)
  return(paste0(sstr[length(sstr)-1], ":",sstr[length(sstr)]))
}