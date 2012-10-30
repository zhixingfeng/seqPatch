correctMovieEffect <- function(alnsF, alnsIdx, movie.effect)
{
	baseline.idx <- sapply(alnsIdx$movieName, function(x, movie.effect.name) which(x==movie.effect.name), movie.effect.name=movie.effect$movie.name)
	for (i in 1:length(baseline.idx)){
		alnsF[[i]]$IPD <- alnsF[[i]]$IPD - movie.effect$movie.IPD.mean[baseline.idx[i]]
		if(i==1000*floor(i/1000))
			cat('processed ',i,' reads','\r')
	}
	cat('processed ',length(baseline.idx),' reads','\n')

	return(alnsF)
}
checkMovieEffect <- function(alnsF, alnsIdx)
{
	movies <- as.character(unique(alnsIdx$movieName))
	movies <- sort(movies)
	movies.read.idx <- list()
	for (i in 1:length(movies)){
        	movies.read.idx[[i]] <- which(alnsIdx$movieName==movies[i])
	}

	movie.IPD.mean <- numeric(length(movies.read.idx))
	movie.IPD.sd <- numeric(length(movies.read.idx))
	movie.IPD.len <- numeric(length(movies.read.idx))

	for (i in 1:length(movies.read.idx)){
		
        	tmp.IPD <- sapply(alnsF[movies.read.idx[[i]]], function(x) sum( x$IPD[!is.na(x$IPD)]) )
        	tmp.IPD.2 <- sapply(alnsF[movies.read.idx[[i]]], function(x) sum( x$IPD[!is.na(x$IPD)]^2 ))
        	tmp.len <- sapply(alnsF[movies.read.idx[[i]]], function(x) sum(!is.na(x$IPD)))
        	movie.IPD.mean[i] <- sum(tmp.IPD)/sum(tmp.len)
        	movie.IPD.sd[i] <- sqrt( sum(tmp.IPD.2)/sum(tmp.len) - (movie.IPD.mean[i])^2 )
        	movie.IPD.len[i] <- sum(tmp.len)

        	cat(movies[i],'\n')
	}
	result <- list()
	result$movie.IPD.mean <- movie.IPD.mean
        result$movie.IPD.sd <- movie.IPD.sd
        result$movie.IPD.len <- movie.IPD.len
	result$movie.name <- movies
	result	
}





