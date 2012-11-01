detectModification.NC <- function(genomeF.native, genomeSeq, context.effect, left.len=6, right.len=1, method='hieModel', is.cvg.all=1)
{
	### check method 
	if (method!='hieModel' & method != 'pooling')
                stop('method should be \'hieModel\' or  \'pooling\'')
	### check left.len and right.len
	if (left.len + right.len + 1 != nchar(names(context.effect)[1])) 
		stop('left.len + right.len + 1 should be equal to length of context\n')

	
	detection <- list()
	for (i in 1:length(genomeF.native$features$ipd_pos) ){
		cat('detection modifications in forward strand : ', names(genomeF.native$features$ipd_pos)[i],'\n' )
		detection$pos[[i]] <- .Call('hieModel_NC', genomeF.native$features$ipd_pos[[i]],as.integer(genomeF.native$genome.start.pos[[i]]),
 			genomeSeq$pos[[names(genomeF.native$features$ipd_pos)[i]]],context.effect, names(context.effect), 
			as.integer(0), as.integer(left.len), as.integer(right.len), method,as.integer(is.cvg.all) )
		names(detection$pos[[i]])[names(detection$pos[[i]])=='t_stat'] <- 't.stat'
		names(detection$pos[[i]])[names(detection$pos[[i]])=='ipd_ratio'] <- 'ipd.ratio'
	}
	names(detection$pos) <- names(genomeF.native$features$ipd_pos)
	for (i in 1:length(genomeF.native$features$ipd_neg) ){
                cat('detection modifications in backward strand : ', names(genomeF.native$features$ipd_neg)[i],'\n' )
		detection$neg[[i]] <- .Call('hieModel_NC', genomeF.native$features$ipd_neg[[i]],as.integer(genomeF.native$genome.start.neg[[i]]),
                	genomeSeq$neg[[names(genomeF.native$features$ipd_neg)[i]]], context.effect, names(context.effect),
	 		as.integer(1), as.integer(left.len), as.integer(right.len), method ,as.integer(is.cvg.all))
		names(detection$neg[[i]])[names(detection$neg[[i]])=='t_stat'] <- 't.stat'
                names(detection$neg[[i]])[names(detection$neg[[i]])=='ipd_ratio'] <- 'ipd.ratio'
	}
	names(detection$neg) <- names(genomeF.native$features$ipd_neg)	
	
	detection$genome.start.pos <- genomeF.native$genome.start.pos
        detection$genome.start.neg <- genomeF.native$genome.start.neg
	
	detection
}

