# get context of significant sites (if cutoff is not NULL, top is not used)
getDetectedContext <- function(detection, genomeSeq, stat='LR_log', top=1000, cutoff=NULL, left.len = 15, right.len = 15)
{
	if ( any(sort(names(detection$pos)) != sort(names(genomeSeq$pos))) )
		stop('unmatched detection and genomeSeq.')
	
	# merage signals
	x <- numeric(0)
	for (i in 1:length(detection$pos)){
		x <- c(x, detection$pos[[i]][[stat]])
	}
	for (i in 1:length(detection$neg)){
                x <- c(x, detection$neg[[i]][[stat]])
        }	
	
	# get cutoff
	if (is.null(cutoff)){
		x.sorted <- sort(x,decreasing=TRUE)
		if (top>length(x.sorted)) stop('top is too large.')
		cutoff <- x.sorted[top]
	}

	# get significant contexts
	detected.context <- character(0)
	for (i in 1:length(detection$pos)){
		idx <-which(detection$pos[[i]][[stat]]>=cutoff)
		cur.ref <- names(detection$pos[i])
		cur.loc <- idx + detection$genome.start.pos[[cur.ref]] - 1
		cur.context <- sapply(cur.loc, function(loc, cur.genomeSeq, left, right) substr(cur.genomeSeq, loc - left, loc + right) ,
			cur.genomeSeq = genomeSeq$pos[[cur.ref]], left=left.len, right=right.len)
		detected.context <- c(detected.context, cur.context)
	} 
	
	for (i in 1:length(detection$neg)){
                idx <-which(detection$neg[[i]][[stat]]>=cutoff)
                cur.ref <- names(detection$neg[i])
                cur.loc <- idx + detection$genome.start.neg[[cur.ref]] - 1
                cur.context <- sapply(cur.loc, function(loc, cur.genomeSeq, left, right) reverseSeq(substr(cur.genomeSeq, loc - right, loc + left)) ,
                	cur.genomeSeq=genomeSeq$neg[[cur.ref]], left=left.len, right=right.len)
                detected.context <- c(detected.context, cur.context)
        }

	detected.context	
}

write.context <- function(context, file)
{
	cat('>seq', 1, '\n', context[1], '\n', sep='', file=file, append=FALSE )
	for (i in 2:length(context)){	
		cat('>seq', i, '\n', context[i], '\n', sep='', file=file, append=TRUE)
	}
}

# (old functions)
bootstrap.motif.freq.bg <- function(motif, genomeSeq, idx.sel, strand, tol, B, n.context)
{
	motif.len <- nchar(motif[1])
	.Call('bootstrap_bg_motif', motif, genomeSeq, as.integer(idx.sel - 1),
	        as.integer(strand), as.integer(tol), as.integer(motif.len), as.integer(B), as.integer(n.context) )

}
getSigfreq <- function(motif, context.sig)
{
        motif.len <- nchar(motif[1])
        sig.count <- list()
        sig.freq <- list()
        for (i in 1:length(motif)){
                cur.idx <- gregexpr(motif[i], context.sig)
                is.find <- sapply(cur.idx, function(x) all(x!=-1))
                sig.count[[ motif[i] ]] <- sum(is.find)
                sig.freq[[ motif[i] ]] <- sig.count[[ motif[i] ]] / length(context.sig)
                cat(i,'\r')
        }
        cat(length(motif),'\n')
        list(sig.count=sig.count, sig.freq=sig.freq, total = length(context.sig))
}




