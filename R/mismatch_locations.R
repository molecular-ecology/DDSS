## Function mismatch_locations.R converts mismatches in bio$blast.mid.var that are designated by “.” and converts them to a number vector related to positions in the read (Numquery) and counts the number of mismatches. This is done prior to checking if there are true matches to degenerate bases. Any improvements or corrections, please contact David Andow.
mismat.loc <- function (df1) {
mis.loc <- (mapply(function(x,y){(y-1)+unlist(gregexpr("[.]", x))},x=df1$blast.mid.var,y=df1$query.POS))
mis.loc <- as.data.frame(t(t(mis.loc))) 
colnames(mis.loc) <- c('mis.loc')
num.mis <- str_count(df1$blast.mid.var, "[.]") # counts the number of mismatches
df2 <- cbind(df1, num.mis, mis.loc)
}
