summaryContextEffect <- function(context.effect.full)
{
	context.effect <- list()
	context.effect <- lapply(context.effect.full, function(x) {rl <- list(); rl$mean <- as.numeric(sapply(x,mean));
								   rl$var <- as.numeric(sapply(x,var)); rl$len <- as.numeric(sapply(x,length));rl } )		
	context.effect

}

getContextEffect <- function(alnsF, left.len, right.len)
{
	if (left.len<1 | right.len < 1)
		stop('left.len and rigth.len should be >= 1 ')
	.Call('getContextEffect',alnsF, as.integer(left.len), as.integer(right.len))
}

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
	if (any(names(genomeSeq$pos) != names(genomeSeq$neg)) | any(names(genomeF$features$ipd_pos) != names(genomeF$features$ipd_neg))
		| any(names(genomeSeq$pos) != names(genomeF$features$ipd_pos)) ){
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



