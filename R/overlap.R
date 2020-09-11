## Function overlap.R, takes the output from mat.degen.base, and for each unique combination of Numquery, hit.index and hsp.index, finds the highest percent identity (per.id), a starting position on the query sequence (start.loc) and the corresponding starting position on the hit sequence (hit.start) for a given overlap length (ovrlap). It returns a list with per.id, hit.start, and ovrlap. It requires input of ovrlap, and several variables passed in: q.from (starting location of alignment on query sequence), q.to (ending location of alignment on query sequence), mml (list of mismatch locations) and hit.start (location of the start of alignment on hit sequence). Any improvements or corrections, please contact David Andow.

overlap <- function(ovrlap, from, to, ml, hs, strand) {
out <- list()
qlen <- to - from +1
if(ovrlap > qlen) {return(c(NA, NA, ovrlap))} #missing value because read is too small
else {if(ml[1] == -1) {return(c(100.00, hs, ovrlap))} #perfect matches
     else {
	qm <- rep(0, qlen)
	for(loc in seq_len(length(ml))) {qm[[ ml[[loc]]-from+1 ]] <- 1}  #vector with 1 at each mismatch location
	ovrlen <- rep(1, ovrlap) # vector of 1s equal to the overlap length
	num.misses <- numeric(0)
	start.loc <- numeric(0)
	# look from the beginning (from)
	x <- head(c(ovrlen, rep(0,qlen)), qlen) # pad the back with 0 to the length of qmm
	start.loc <- c(start.loc, from) #start location of the first overlap section
	num.misses  <- c(num.misses, x %*% qm) #number of misses in the first ovrlap bps
	#look from the end (to)
	x <- tail(c(rep(0,qlen), ovrlen), qlen) # pad the front with 0 to the length of qmm
	start.loc <- c(start.loc, to - ovrlap) #start location of the last overlap section
	num.misses  <- c(num.misses, x %*% qm) #number of misses in the last ovrlap bps
	for(miss in seq_len(length(ml))) #look relative to each ml[[miss]]
		{
		if(ml[[miss]] <= (to-ovrlap)) # look toward the end from ml[[miss]]
			{
			x <- tail(c(rep(0,qlen), ovrlen), ovrlap + ml[[miss]]-from+1) # pad the front with 0
			x <- head(c(x, rep(0, qlen)), qlen) # pad the back with 0 to the length of qmm
			num.misses  <- c(num.misses, x %*% qm) #number of misses in ovrlen
			start.loc <- c(start.loc, ml[[miss]]+1) #start location of the ‘miss’ overlap section
			} #end if ml[[miss]]+ovrlap <= to
		if(ml[[miss]] >= (from+ovrlap)) #look toward the beginning from ml[[miss]]
			{
			x <- tail(c(rep(0,qlen), ovrlen), ml[[miss]]-from) # pad the front with 0
			x <- head(c(x, rep(0, qlen)), qlen) # pad the back with 0 to the length of qmm
			num.misses  <- c(num.misses, x %*% qm) #number of misses in ovrlen
			start.loc <- c(start.loc, ml[[miss]]-ovrlap) #start location of the ‘miss’ overlap section
			} #end if ml[[miss]]-ovrlap >= from
		} # end for miss in ml
	## Smallest value in num.misses (minmis) to calculate highest per.id for a given overlap length
	start <- which(num.misses %in% min(num.misses))[[1]]  # first or only index of minmis
	if (strand == "+") {
		hit.start <- hs + (start.loc[[start]] - from) } else { 
		hit.start <- hs - (start.loc[[start]] - from) }
	out <- c(100*(ovrlap - min(num.misses))/ovrlap, hit.start, ovrlap)
	} # end else ml <> -1 
  } #end ovrlap <= qlen 
out
} # end overlap
