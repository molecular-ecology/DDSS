## Function best_matches.R finds the unique matching hit(s) with the highest percent identity for a read (Numquery), for a given overlap length and marks the ones greater than a percent identify threshold.
# Function inputs are df.in (data frame) and per, a list of percent identity thresholds. df.in has column_names <- c("Numquery", "hit", "hit.index", "hsp.index", "initial.hit.pos", "STRAND", "blast.align.length", "num.mis", "per.id", "hit.start.loc", "overlap")
#Function will "clean up" the hit names so that they are just genus species names, and put these names in a new column. This version will fix the misspellings in the hit names for C sanguinea and H tredecimpunctata. Any improvements or corrections, please contact David Andow.
#It will output a dataframe and a tibble. The dataframe has the hits with the highest percent identity for each overlap length sorted by overlap length, eliminating conditions with NA for hit matches. It will create a new column for each of the percent identity thresholds, indicating if the read exceeds the threshold (1) or not (0). It will have the read id (Numquery), overlap length, original percent identity, indeces for per thresholds and the hit match(es) with the original hit and just genus species. If there is more than one best hit for a read, they are concatenated and separated by a pipe.
#   The tibble will have a summary of the data indicating the number of reads associated with each overlap length, percent identity threshold, and hit species.

best.match <- function (df.in, per) {
# Sort the dataframe by overlap, Numquery, hit.index, hsp.index
df.in <- df.in[with(df.in, order(overlap, Numquery, hit.index, hsp.index)), ] 

#Fix names for Cycloneda sanguinea and Hippodamia tredecimpunctata
df.in$hit <- str_replace_all(df.in$hit, "sanguinesa", " sanguinea")
df.in$hit <- str_replace_all(df.in$hit, "tredecempunctata","tredecimpunctata")

#Clean up hit names so that they only “genus species” names
hitsp <- str_replace_all(df.in$hit, "UNVERIFIED: ", "")
hitsp <- str_replace_all(hitsp, "_", " ")
hitsp <- str_replace_all(hitsp, "  ", " ")
hitsp <- str_extract(hitsp, "[^ ]+ [^ ]+")
df.in <-cbind(df.in, hitsp)
df.in$hitsp <- as.character(df.in$hitsp)

#Remove all records with NA for hit, no hits for a given overlap length and identity threshold
df.in <- df.in[!is.na(df.in$hit),]

#Create a dummy variable to hold unique combinations of overlap and Numquery and split the dataframe
df.in <- cbind(df.in, "dum" = paste(df.in$overlap, df.in$Numquery, sep="-"))
dum.uni <- unique(df.in$dum)
dfnew <- split(df.in, df.in$dum)
cols <- c("Numquery", "overlap", "per.id", "hit.match", "hit.match.sp")
mylist <- c()  # NULL list

i.index <- seq_len(length(dum.uni))

for(i in i.index) {  
	# record indeces for dfnew[[i]] for the highest per.id and keep unique names
	keep.list <- c(which(dfnew[[i]]$per.id == max(dfnew[[i]]$per.id)))  # record indeces with the maximum per.id for the read and overlap 
	hit.match <- paste(unique(c(unlist(dfnew[[i]]$hit[keep.list]))), collapse=" | ")
	hit.match.sp <- paste(unique(c(unlist(dfnew[[i]]$hitsp[keep.list]))), collapse=" | ")
	mylist[[i]] <-c(dfnew[[i]]$Numquery[[1]], dfnew[[i]]$overlap[[1]], max(dfnew[[i]]$per.id), hit.match, hit.match.sp)
} # end for i, read
# compile mylist into a dataframe
df.out <- do.call("rbind.data.frame", mylist) 
colnames(df.out) <- cols 

# convert variable types from factors
df.out[,cols] <- apply(df.out[,cols], 2, function(x) as.character(x))
df.out$overlap <- as.integer(df.out$overlap)
df.out$per.id <- as.numeric(df.out$per.id)

# sort the dataframe by overlap and remove na for matches
df.out <- df.out %>% arrange(overlap) %>% filter(!is.na(hit.match.sp))

# mark the reads that exceed the percent identity thresholds
namevector <- paste("per.id.", per, sep = "")
for(i in seq_len(length(namevector))) {
	blank <- rep(0, length(df.out$per.id)) # list of 0
	ok <- which(df.out$per.id >= per[i])   # indeces of Numquerys with per.id greater than or equal to a threshold
	blank[ok] <-1
	assign(namevector[i], blank)
} # end for i namevector
pi <- as.data.frame(do.call(cbind, mget(namevector))) # mget is needed to get the values in the names
df.out <- cbind(df.out,pi)
df.out[is.na(df.out)] <- 0  # Replacing NA values with 0

# summarize the data
catagories <- df.out %>% arrange(overlap) %>% group_by(overlap, hit.match.sp) #dplyr, sorts, and groups the data

result <- c()  #Null list
for(i in seq_len(length(per))) {
	temp <- catagories %>% filter(get(namevector[i]) == 1) %>% summarize(id.thres = mean(per.id), num.reads = n()) 
	temp$id.thres <- per[i] 
	result[[i]] <- temp 
} #end for i in per
result <- bind_rows(result)

list(df.out, result)

}
