# report IPD distribution and number of sites of motifs
report.motif.stat <- function(detection, genomeSeq, motif, shift, ref.name, prefix, stat='LR_log', mar=c(10,4,4,4))	
{
	if (length(motif)!=length(shift)) stop('incompatible motif and shift')

	n <- 1
	stat.motif.list <- list()
	for (i in 1:length(motif)){
		idx.motif <-  getIdxByMotif(genomeSeq, motif[i], shift[i])
		stat.motif <- getStatByIdx(detection, idx.motif, stat) 
		stat.motif.list[[n]] <- stat.motif$pos[[ref.name]]
		stat.motif.list[[n+1]] <- stat.motif$neg[[ref.name]] 
		names(stat.motif.list)[n] <- paste(motif[i],'_pos',sep='')
		names(stat.motif.list)[n+1] <- paste(motif[i],'_neg',sep='')
		n <- n + 2
	}
	len <- length(stat.motif.list)
	stat.motif.list [[len+1]] <- detection$pos[[ref.name]][[stat]]
	names(stat.motif.list)[len + 1] <- 'All'

	file.path <- paste(prefix, '_', ref.name, '.pdf', sep='')
	pdf(file=file.path)
	par(mar=mar)
	boxplot(stat.motif.list, outline=FALSE, las=2, ylab=stat)
	dev.off()
}
getDetectedContextByIndex <- function(detected.sites, genomeSeq, left.len=15, right.len=15, is.merge=TRUE) 
{
	detected.context <- list()
	detected.context$pos <- list()
	detected.context$neg <- list()
	# forward strand 
	for (i in 1:length(detected.sites$pos)){
		cur.ref <- names(detected.sites$pos[i])
		cur.genomeseq <- genomeSeq$pos[[i]]
		detected.context$pos[[cur.ref]] <- character(0)
		cur.idx <- detected.sites$pos[[i]]
		for (j in 1:length(cur.idx)){
			detected.context$pos[[cur.ref]] <- c(detected.context$pos[[cur.ref]], 
				substr(cur.genomeseq, cur.idx[j]-left.len, cur.idx[j]+right.len))
		}
	}

	# backward strand
        for (i in 1:length(detected.sites$neg)){
                cur.ref <- names(detected.sites$neg[i])
                cur.genomeseq <- genomeSeq$neg[[i]]
                detected.context$neg[[cur.ref]] <- character(0)
                cur.idx <- detected.sites$neg[[i]]
                for (j in 1:length(cur.idx)){
                        detected.context$neg[[cur.ref]] <- c(detected.context$neg[[cur.ref]],
                                reverseSeq(substr(cur.genomeseq, cur.idx[j]-right.len, cur.idx[j]+left.len)))
                }
        }
	if (is.merge==TRUE){
		detected.context.pool <- character(0)
		for (i in 1:length(detected.context$pos)){
			for (j in 1:length(detected.context$pos[[i]]))
				detected.context.pool <- c(detected.context.pool, detected.context$pos[[i]][j])
		}		
		for (i in 1:length(detected.context$neg)){
                        for (j in 1:length(detected.context$neg[[i]]))
                                detected.context.pool <- c(detected.context.pool, detected.context$neg[[i]][j])
                }
		return(detected.context.pool)
	}else{
		return(detected.context)
	}
	
}
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

# write context into a fasta file
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




