getTimeBias <- function(alnsF, alnsIdx, start.time, num.bin = 100)
{
	#library(seqPatch)
	#load('~/dump/getTimeBiasTest.Rdata')
	#start.time <- start.time.tmp
	#num.bin <- 100
	
	time.min <- min(sapply(start.time, min, na.rm=TRUE))
	time.max <- max(sapply(start.time, max, na.rm=TRUE))
	time.min <- 0
	step.len <- (time.max - time.min)/num.bin
	movie.names <- unique(alnsIdx$movieName)
	
	IPDbyBin <- list()
	for (i in 1:length(movie.names)){
		idx.m <- which(alnsIdx$movieName==movie.names[i])
		alnsF.m <- alnsF[idx.m]
		start.time.m <- start.time[idx.m]	
		cat('movie : ',movie.names[i], '\n')
		IPDbyBin[[i]] <- .Call('getTimeBias', alnsF.m, start.time.m, step.len, as.integer(num.bin))
	}
	names(IPDbyBin) <- movie.names
	IPDbyBin
}






