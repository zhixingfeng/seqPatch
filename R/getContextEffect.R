summaryContextEffect <- function(context.effect.full)
{
	context.effect <- list()
	context.effect <- lapply(context.effect.full, function(x) {rl <- list(); rl$mean <- as.numeric(sapply(x,mean));
								   rl$var <- as.numeric(sapply(x,var)); rl$len <- as.numeric(sapply(x,length));rl } )		
	context.effect

}

#getContextEffect <- function(alnsF, left.len, right.len)
#{
#	if (left.len<1 | right.len < 1)
#		stop('left.len and rigth.len should be >= 1 ')
#	.Call('getContextEffect',alnsF, as.integer(left.len), as.integer(right.len))
#}

getContextEffect.mean <- function(context.effect)
{
	lapply(context.effect, function(x) {  sum(sapply(x,sum, na.rm=TRUE) )/sum(sapply(x,function(t) sum(!is.na(t))))  }  )
}
getContextEffectForRead <- function(alnsF, context.effect, left.len, right.len)
{
	if (left.len<1 | right.len < 1)
                stop('left.len and rigth.len should be >= 1 ')
	if (left.len + right.len + 1 != nchar(names(context.effect)[1]))
		stop('left.len + rigth.len + 1 should be equal to context length\n')
	.Call('getContextEffectForRead',alnsF, as.numeric(context.effect), names(context.effect), as.integer(left.len), as.integer(right.len) )
	
}
getContextEffectByPos <- function(genomeF, genomeSeq, left.len=6, right.len=1, min.cvg=15, min.pos=10, is.filter=TRUE, is.summary=TRUE, is.extend=FALSE)
{
	if (any(names(genomeSeq$pos) != names(genomeSeq$neg)) | any(names(genomeF$features$ipd_pos) != names(genomeF$features$ipd_neg)) ){
		cat('Error : refSeq of forward strand and negative strand should be the same\n')
		return (NULL)
	}

		
	# get context effect for each refseq
	context.effect.byrefseq <- list()
	context.ex.byrefseq <- list()
	for (i in 1:length(genomeF$features$ipd_pos)){
		cur_refSeq <- names(genomeF$features$ipd_pos)[i]
		cat('get context effect of ', cur_refSeq, '\n')
		ipd_pos <- genomeF$features$ipd_pos[[i]]	
		ipd_neg <- genomeF$features$ipd_neg[[i]]
		if (is.null(genomeSeq$pos[[cur_refSeq]]) | is.null(genomeSeq$neg[[cur_refSeq]]))
			stop('could not find ',cur_refSeq,' in genomeSeq\n')
		genomeSeq_pos <- genomeSeq$pos[[cur_refSeq]]
		genomeSeq_neg <- genomeSeq$neg[[cur_refSeq]]
		if (is.null(genomeSeq_pos) | is.null(genomeSeq_neg)){
			cat('Error: can not find \'',cur_refSeq, '\' in the genome.\n', sep='')
			return(NULL)	
		}
		context.effect.byrefseq[[cur_refSeq]] <- .Call('getContextEffectByPos', ipd_pos, ipd_neg, genomeSeq_pos, genomeSeq_neg, as.integer(genomeF$genome.start.pos[[i]]), as.integer(genomeF$genome.start.neg[[i]]), as.integer(left.len), as.integer(right.len) )
		if (is.extend==TRUE)
			context.ex.byrefseq[[cur_refSeq]] <- .Call('getContextExByPos', ipd_pos, ipd_neg, genomeSeq_pos, genomeSeq_neg, as.integer(genomeF$genome.start.pos[[i]]), as.integer(genomeF$genome.start.neg[[i]]), as.integer(left.len), as.integer(right.len) )
	}

	# merge context 
	context.effect <- context.effect.byrefseq[[1]]
	if(is.extend==TRUE) context.ex <- context.ex.byrefseq[[1]]

	if (length(context.effect.byrefseq)>=2){
		for (i in 2:length(context.effect.byrefseq)){
			context.effect <- .Call('getContextEffectByPosMerge', context.effect
			, context.effect.byrefseq[[i]], names(context.effect), names(context.effect.byrefseq[[i]]) )
		}
		if (is.extend==TRUE){
			for (i in 2:length(context.ex.byrefseq)){
		                context.ex <- .Call('getContextExByPosMerge', context.ex
                	        , context.ex.byrefseq[[i]], names(context.ex), names(context.ex.byrefseq[[i]]) )
        		}
	
		}
	}
	
	# remove low coverage positions and contexts with too few positions
	if (is.filter == TRUE){
		cat('remove low coverage positions and contexts with too few positions.\n')
		context.effect <- trim.context.effect(context.effect, min.cvg)
		len <- sapply(context.effect, length)
		context.effect <- context.effect[len>=min.pos]
	}

	if (is.extend == FALSE){
		if (is.summary == TRUE){
			return (summaryContextEffect(context.effect))
		}else{
			return (context.effect)
		}
	}else{
		return (list(context.effect=context.effect, context.ex = context.ex))	
	}
}


#mergeContextEffect <- function(context.effect.list)
#{
#	context.all <- character()
#	# get context union
#	for (i in 1:length(context.effect.list)){
#		context.all <- union(context.all, names(context.effect.list[[i]]))	
#	}
#
#	# get context effect union
#	idx <- which (!grepl('N',context.all))
#	context.effect.all <- list()
#	for (i in 1:length(idx)){
#		cur.context <- context.all[idx[i]]
#		context.effect.all[[cur.context]]$mean <- numeric(0)
#		context.effect.all[[cur.context]]$var <- numeric(0)
#		context.effect.all[[cur.context]]$len <- numeric(0)
#		for (j in 1:length(context.effect.list)){
#			context.effect.all[[cur.context]]$mean <- c(context.effect.all[[cur.context]]$mean,
#									context.effect.list[[j]][[cur.context]]$mean)
#			context.effect.all[[cur.context]]$var <- c(context.effect.all[[cur.context]]$var,
#                                                                        context.effect.list[[j]][[cur.context]]$var)
#			context.effect.all[[cur.context]]$len <- c(context.effect.all[[cur.context]]$len,
#                                                                        context.effect.list[[j]][[cur.context]]$len)
#		}
#		if (i == 100*floor(i/100)){
#			cat(i,'\r')
#		}	
#	}
#	cat(i,'\n')
#	context.effect.all
#}

mergeContextEffect <- function(context.effect.list)
{
	if (length(context.effect.list)<=1) return(context.effect.list[[1]])
	context.effect.all <- context.effect.list[[1]]
	for (i in 2:length(context.effect.list)){
		context.effect.all <- .Call('mergeContextEffectSummary', context.effect.all, context.effect.list[[i]], 
					names(context.effect.all), names(context.effect.list[[i]] ))
	}
	idx <- which (!grepl('N', names(context.effect.all)))
	context.effect.all[idx]
		
}


getContextEffectSet <- function(genomeF, genomeSeq, left.len=c(7,6,6), right.len=c(2,2,1), min.cvg=15, min.pos=10, is.filter=TRUE, is.summary=TRUE, is.extend=FALSE)
{
	if (length(left.len)!=length(right.len)) stop('length of left.len and right.len should be the same.')

	context.effect.set <- list()	
	for (i in 1:length(left.len)){
		cat('context_effect_',left.len[i],'_',right.len[i],'\n', sep='')
		context.effect.set[[i]] <- getContextEffectByPos(genomeF, genomeSeq, left.len[i], right.len[i],
					 min.cvg, min.pos, is.filter, is.summary, is.extend)		
	}

	names(context.effect.set) <- paste('context_effect_', left.len, '_', right.len, sep='')
	context.effect.set[[length(left.len)+1]] <- left.len
	context.effect.set[[length(left.len)+2]] <- right.len
	names(context.effect.set)[length(left.len) + 1:2] <- c('left.len', 'right.len')
	context.effect.set
}

mergeContextEffectSet <- function(context.effect.set.list)
{
	if (length(context.effect.set.list)<=1) return(context.effect.set.list[[1]])
	# check consistency of left.len and right.len	
	for (i in 2:length(context.effect.set.list)){
		if (any(context.effect.set.list[[i]]$left.len != context.effect.set.list[[i-1]]$left.len) |
		    any(context.effect.set.list[[i]]$right.len != context.effect.set.list[[i-1]]$right.len)|
		    length(context.effect.set.list[[i]])-2 != length(context.effect.set.list[[i]]$left.len)|
		    length(context.effect.set.list[[i-1]])-2 != length(context.effect.set.list[[i-1]]$left.len) ){	
			stop('inconsistent left.len or right.len.')
		}
	}

	context.effect.set.pool <- list()	
	for (i in 1:length(context.effect.set.list[[1]]$left.len)){
		cat('merge ', names(context.effect.set.list[[1]])[i], '\n', sep='')
		cur.pool <- list()
		for (j in 1:length(context.effect.set.list)){
			cur.pool[[j]] <- context.effect.set.list[[j]][[i]]
		}
		context.effect.set.pool[[i]] <- mergeContextEffect(cur.pool)
	
	}		
	names(context.effect.set.pool) <- names(context.effect.set.list[[1]])[1:length(context.effect.set.list[[1]]$left.len)]
	context.effect.set.pool$left.len <- context.effect.set.list[[1]]$left.len
	context.effect.set.pool$right.len <- context.effect.set.list[[1]]$right.len

	context.effect.set.pool
}

filterContextEffect <- function(context.effect, min.len = 10)
{
	# remove context with too few sites
	len <- sapply(context.effect, function(x) length(x$mean))		
	context.effect[len>=min.len]
}

filterContextEffectSet <- function(context.effect.set, min.len = 10)
{
        # remove context with too few sites
	for (i in 1:(length(context.effect.set)-2)){
		context.effect.set[[i]] <- filterContextEffect(context.effect.set[[i]], min.len)
	}
	context.effect.set
}


# get hyperprameters of each context
getContextHyperparameters <- function(context.effect, max.iter=10)
{
	lapply(context.effect, function(x, max.iter) {rl <- hieModel_core(x$mean, x$var, x$len, max.iter); 
				hyperPara <- as.list(rl$hyperPara); names(hyperPara)=c('theta','kappa','upsilon','tau2');
				hyperPara$num.sites <- length(x$mean); hyperPara$median.cvg <- median(x$len);
				hyperPara
				}, max.iter=max.iter)		
}






