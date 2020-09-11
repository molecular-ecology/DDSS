## Function percent_identity.R calculates the highest percent identity in a query for a given hit.index and hsp.index, here uniquely identified as Numquery, hit.index, and hsp.index, for specified overlap lengths.
#Inputs are a df and a list of overlap lengths. The df has the structure of Numquery, hit, hit.index, hsp.index, STRAND, hit.POS, query.POS, blast.query.from, blast.query.to, blast.align.length, num.mis, and mis.loc.
#Will remove NAs for percent identity and output the highest percent identity, and starting location in the query and hit sequence of the overlap with the highest percent identity, and the overlap length. Calculates the hit position (in the reference) of the first bp in the alignment (initial.hit.pos). It retains Numquery, hit, hit.index, hsp.index, initial.hit.pos, STRAND, blast.align.length, and num.mis. Any improvements or corrections, please contact David Andow.

per.id <- function (df.in, ov) {
# Sort the combined dataframe by Numquery, hit.index, hsp.index, query.POS
df.in <- df.in[with(df.in, order(Numquery, hit.index, hsp.index, query.POS)), ] 
#Create a dummy vaiable to hold unique combinations of Numquery, hit.index and hsp.index
df.in <- cbind(df.in, "dum" = paste(df.in$Numquery, df.in$hit.index, df.in$hsp.index, sep="-"))
dum.uni <- unique(df.in$dum)
mylist <-c() #NULL list
dfnew <- split(df.in, df.in$dum)
#Create column names for the output data frame
cols <- c("Numquery", "hit", "hit.index", "hsp.index", "initial.hit.pos", "STRAND", "blast.align.length", "num.mis", "per.id", "hit.start.loc", "overlap")
#Loop over all i of the unique combinations of Numquery, hit.index and hsp.index
for(i in seq_len(length(dum.uni))) {
	out.list = c() #Null list
	##calculate the initial hit position, which depends on STRAND
	if (dfnew[[i]]$STRAND[[1]] == "+") {
		initial.hit.pos <- dfnew[[i]]$hit.POS[[1]] - (dfnew[[i]]$query.POS[[1]] + dfnew[[i]]$blast.query.from[[1]])} else { 
		initial.hit.pos <- dfnew[[i]]$hit.POS[[1]] + (dfnew[[i]]$query.POS[[1]] - dfnew[[i]]$blast.query.from[[1]])}
	## For j overlap length in ov
	for(j in seq_len(length(ov))) {
		#calculate the output vector for ov[[j]]
		out.list[[j]] <- c(dfnew[[i]]$Numquery[[1]], dfnew[[i]]$hit[[1]], dfnew[[i]]$hit.index[[1]], dfnew[[i]]$hsp.index[[1]], initial.hit.pos, dfnew[[i]]$STRAND[[1]], dfnew[[i]]$blast.align.length[[1]], sum(dfnew[[i]]$num.mis),
		overlap(ov[[j]], 
			dfnew[[i]]$blast.query.from[[1]], 
			dfnew[[i]]$blast.align.length[[1]]+dfnew[[i]]$blast.query.from[[1]]-1, 
			unlist(dfnew[[i]]$mis.loc), 
			initial.hit.pos,
			dfnew[[i]]$STRAND[[1]]) #Function overlap returns a list of per.id, hit.start, and ov[[j]]
		) # end concatenate for out.list[[j]]
		} #end for j in ov
	mylist[[i]] <- do.call("rbind.data.frame",out.list) #combine all vectors in out.list into a data frame
	colnames(mylist[[i]]) <- cols 
# Change the type of each column from Factor to character
mylist[[i]][,cols] <- apply(mylist[[i]][,cols], 2, function(x) as.character(x))
	} #end for i, dum.uni
df.out <- bind_rows(mylist) #combine all data frames into a single data frame using dplyr

# Change the type of each column from Factor to the appropriate type
cols <- c("hit.index", "hsp.index", "initial.hit.pos", "blast.align.length", "num.mis", "hit.start.loc", "overlap")
df.out[,cols] <- apply(df.out[,cols], 2, function(x) as.integer(x))
df.out$per.id <- as.numeric(df.out$per.id)
df.out <- df.out[with(df.out, order(overlap, Numquery, hit.index, hsp.index)), ] 
#remove NAs in per.id
df.out <- df.out[!is.na(df.out$per.id),]
df.out
}
