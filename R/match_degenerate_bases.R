## Function match_degenerate_bases.R takes the custom output from Roberto’s Java script, checks each record for degenerate bases in the hit sequence from the reference database, and eliminates records for which the query matches all of the degenerate bases in a record. Output is the same data structure as the original with two additional columns. num.mis is the number of mismatched base pairs in the record. This will be used to calculate the percent identity of a query sequence to a hit sequence for a given overlap length. mis.loc is a vector of integers for the record that indicates where the mismatches occurred with the value 1 being the same position as query.POS. This will be used to determine percent identity for a given overlap length in the alignment of the query and hit sequence. Any improvements or corrections, please contact David Andow.
# Input to mat.degen.base is a dataframe constructed from the custom output.
# Output from mat.degen.base is a new dataframe with the same structure as the input dataframe, eliminating records with no mismatches after considering matches between degenerate bases in hit sequences and bases in query sequence, reducing the number of mismatches in the record for matches with degenerate bases, and keeping only locations with mismatches in the record.

mat.degen.base <- function(df1) {
degen.mat.list <- list(R= c('A', 'G'), Y= c('C', 'T'), N= c('A', 'T', 'C', 'G'), W= c('A', 'T'), S= c('C', 'G'), M= c('A', 'C'), K= c('T', 'G'), B= c('T', 'C', 'G'), H= c('A', 'T', 'C'), D= c('A', 'T', 'G'), V= c('A', 'C', 'G') )
degen.list <- (list('R', 'Y', 'N', 'W', 'S', 'M', 'K', 'B', 'H', 'D', 'V'))
patR <- 'R'; patY <- 'Y'; patN <- 'N'; patW <- 'W'; patS <- 'S'; patM <- 'M'; patK <- 'K'; patB <- 'B'; patH <- 'H'; patD <- 'D'; patV <- 'V';
subhit <- df1$REF.hit.
# creates a logical vector with length = nrows of bio; TRUE for any record with a degenerate base, FALSE if all are normal bases
degen <- (grepl(patR, subhit)|grepl(patY, subhit)|grepl(patN, subhit)|grepl(patW, subhit)|grepl(patS, subhit)|grepl(patM, subhit)|grepl(patK, subhit)|grepl(patB, subhit)|grepl(patH, subhit)|grepl(patD, subhit)|grepl(patV, subhit) )
df2 <- df1[!degen,] #selects the rows with normal bases (and therefore, mismatches)
degen.data <- df1[degen,] # selects the rows with degenerate bases
mismat.degen.data <-data.frame() #empty dataframe to store the rows of degen.data with mismatches
for(i in seq_len(nrow(degen.data))) # processes degen.data line by line (O(n2) slow)
{
subhit <- degen.data[i, 'REF.hit.']
subque <- degen.data[i, 'ALT.query.']
nch <- nchar(subhit)
ncq <- nchar(subque)
if(nch>1)
{
subhit <- substring(subhit, 1:nch, 1:nch) # split REF(hit) into a vector of single characters
subque <- substring(subque, 1:ncq, 1:ncq) # split ALT(query) into a vector of single characters
}
# check for correct matches to the degenerate bases
degen.hit <- which(subhit %in% degen.list) #locations in subhit with degenerate base
mat.list <- logical(length=length(degen.hit)) # initializes logical vector to FALSE to keep mismatches
for(j in seq_len(length(degen.hit))) # 
{
degen.tar <- which(degen.list %in% subhit[[ degen.hit[[j]] ]]) #identity of degenerate base
if(degen.hit[[j]] <= length(subque)) 
mat <- subque[[degen.hit[[j]] ]] %in% degen.mat.list[[degen.tar]] #base match
else {mat <- FALSE} # length difference means there is no match
if(mat) 
{degen.data[i, 'num.mis'] <- degen.data[i, 'num.mis'] - 1 #reduce mismatch count
mat.list[[j]] <- mat} # indexes in degen.hit with matches
} #end for j
degen.data[[i, 'mis.loc']] <- degen.data[[i, 'mis.loc']][!mat.list] #drops indexes with matches
if(degen.data[i, 'num.mis']>0) mismat.degen.data <- rbind(mismat.degen.data,degen.data[i,]) #some do not match, including normal bases
} #end for i
df2 <- rbind(df2, mismat.degen.data)
df2[with(df2, order(Numquery, hit.index, hsp.index, query.POS)), ]; # sort df2
}
