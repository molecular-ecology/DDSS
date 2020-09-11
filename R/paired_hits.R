# Function paired_hits.R restricts hits to paired reads with the same species identified. Any improvements or corrections, please contact David Andow. 
# Input “best.csv” file, vectors of overlap and per.id, and the name of the library
# Output a dataframe with the number of reads for each species match, including those matching multiple hit species
paired.hits <- function(df.in, ov, id, name) {
ov.index <- seq_len(length(ov))
id.index <- seq_len(length(id))

# drop hit.match, sort hit.match.sp alphabetically (involves splitting hit.match.sp so that each species is separated, sorting this list for each hit, pasting it back together into one string, and unlisting), pull out identifiers from Numquery.
df.in <- df.in %>% as.data.frame %>% select(-hit.match) %>% select(-matches("per.id.")) %>% 
   mutate(new = strsplit(hit.match.sp, fixed(" \\| ")) %>% lapply(.,'sort') %>% 
   lapply(. %>% paste(., collapse =' | ')) %>% unlist(.) ) %>% 
   mutate(outcome = str_sub(map_chr(Numquery, function(s) paste0(strsplit(s, ":")[[1]][5:7], collapse=":")), end=-3))

#split by overlap
df.in <- split(df.in, df.in$overlap)
df.ov <- list()
for(i in ov.index) {
df.id <- list()
for(j in id.index) {

# Select rows with paired reads and the same hit.match.sp (now recorded in new), summarize by hit species (new) and read number
df.id[[j]] <- df.in[[i]] %>% filter(overlap == ov[i]) %>% filter(per.id >= id[j]) %>% group_by(outcome, new) %>% filter(n()>1) %>% ungroup %>% count(new) %>% mutate(threshold = paste0(ov[i],"_",id[j]))
} # end for(j in id.index)
df.ov[[i]] <- bind_rows(df.id, .id=NULL)
} # end for(i in ov.index)

# bind rows and convert from long to wide for output
df.out <- bind_rows(df.ov, .id=NULL) %>% mutate(library = name) %>% pivot_wider(   names_from  = c(new),  values_from = c(n)   ) %>% mutate_all(~replace(., is.na(.), 0))
}

