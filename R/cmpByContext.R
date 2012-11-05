cmpByContext <- function(data.1, data.2)
{
	motif <- intersect(names(data.1), names(data.2))	
	rl <- list()	
	rl$x <- sapply(motif, function(x, cur.data) cur.data[[x]],cur.data=data.1)
	rl$y <- sapply(motif, function(x, cur.data) cur.data[[x]],cur.data=data.2)
	rl	
}

#cmpByContext <- function(data.1, data.2)
#{
#	motif.len <- nchar(names(data.1)[1])
#	motif <- enum_seq(motif.len)
#	x <- list()
#	y <- list()
#	idx <- 1
#	for (i in 1:length(motif)){
#		if (i == 1000*floor(i/1000))
#			cat('processed ',i,' motifs\r')
#		cur.motif <- motif[i]
#		tmp.1 <- data.1[[cur.motif]]
#		tmp.2 <- data.2[[cur.motif]]
#		if (is.null(tmp.1) | is.null(tmp.2))
#			next
#		if (is.na(tmp.1) | is.na(tmp.2))
#			next
#		x[idx] <- tmp.1
#		y[idx] <- tmp.2
#		idx <- idx + 1
#	}
#	cat('processed ', length(motif), ' motifs\n')
#	result <- list()
#	result$x <- as.numeric(x)
#	result$y <- as.numeric(y)
#	result
#}


