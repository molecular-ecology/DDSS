## Best hit
## Main script to process data from BlastNToSnp with custom script to acquire 100% identity matches. It corrects text in characters, checks degenerate bases and eliminates records with matches, corrects the number of mismatches, calculates best percent identity for a given overlap length on each alignment. It gives the hit with the highest percent identity for each overlap length, and summarizes the data by giving the number of reads for each overlap length, percent identity threshold and hit species. This runs with time complexity O(n), about 1,000,000 lines in 2500 secs. The functions are in a subdirectory of the working directory named “R”. Any improvements or corrections, please contact David Andow.
#set environment
library(stringr)
library(dplyr)
library(stringi)
library(purrr)
library(tibble)
library(tidyr)
setwd("your directory")  #### CHANGE
source("correct_text.R")
source("match_degenerate_bases.R")
source("mismatch_locations.R")
source("overlap.R")
source("percent_identity.R")
source("best_matches.R")
source("paired_hits.R")

## read and process data
ov <- c(your list of desired overlap length(s)) ####  CHANGE
per <- c(your list of desired percent identity(ies)) #### CHANGE
short.name <- #your designated short name for input file #### CHANGE
bio <-  read.delim(filename, header = TRUE, na.strings= "xyz", stringsAsFactors=FALSE) #### CHANGE filename
bio <- correct.text(bio) # replace commas, standardize column names
bio <- mismat.loc(bio) #determines the locations and number of potential mismatches
newbio <- mat.degen.base(bio) #determines if a potential mismatch to a degenerate base is really a mismatch, and removes all false mismatches
mat <- per.id(newbio, ov) #determines best percent identity in an alignment for overlap length specified in ov
write.csv(mat, file = paste0(short.name,"mat.csv"),row.names=FALSE)
best <- best.match(mat, per)
match.list <- best[1]
write.csv(match.list, file = paste0(short.name, "_best.csv"),row.names=FALSE)
result <- best[2]
write.csv(result, file = paste0(short.name, "_result.csv"),row.names=FALSE)
pair <- paired.hits(match.list, ov, per, name)
write.csv(pair, file = paste0(short.name, "_pair_result.csv"),row.names=FALSE)

