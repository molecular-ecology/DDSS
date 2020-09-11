## Function correct_text.R coverts all commas to space in the specified character strings, so that they do not interfere with CSV files later on, and changes the names of the dataframe to a consistent naming system. Any improvements or corrections, please contact David Andow.
correct.text <- function (df1) {
library(stringr)
#change names 
cnames <- c("Numquery", "hit", "hit.index", "hsp.index", "query.POS", "hit.POS", "STRAND", "REF.hit.", "ALT.query.", "blast.align.length", "blast.perc.identity", "blast.hit.var", "blast.query.var", "blast.mid.var", "blast.query.from", "blast.query.to")
colnames(df1) <- mapply(function(x,y) {x <- y}, x =colnames(df1), y =cnames)
#change commas
df1$Numquery <- str_replace_all(df1$Numquery, ",", " ")
df1$hit <- str_replace_all(df1$hit, ",", " ")
df1
}
