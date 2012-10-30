hieModelEB <- function(genomeF.native, genomeF.ctrl=NULL, genomeSeq, context.effect=NULL, left.len=6, right.len=1,
		 is.filter.sample=FALSE, is.filter.context.effect=FALSE, 
		is.trim.cvg = FALSE, is.trim.pos = FALSE, min.cvg=5, min.pos=4,  method='hieModel', is.perm=FALSE)
{
	if (method!='hieModel.nc' & method!='hieModel' & method != 'CC' &  method !='pooling')
		stop('method should be \'hieModel.nc\', \'hieModel\', \'CC\' or \'pooling\'')
	if (method!='hieModel.nc' & is.null(genomeF.ctrl))
		stop('genomeF.ctrl should be specified')
	if (method=='hieModel.nc' & !is.null(genomeF.ctrl))
		warning('genomeF.ctrl would not be used when method=\'hieModel.nc\'')
	### remove positions with low coverage
	if (is.trim.cvg == TRUE & !is.null(context.effect)){
		cat('remove low coverage position from context effect\n')
        	context.effect <- trim.context.effect(context.effect, min.cvg)
	}
	### remove context with too few positions 
	if (is.trim.pos == TRUE & !is.null(context.effect)){
		cat('remove context with too few positions\n')
		len <- sapply(context.effect, length)
		context.effect <- context.effect[len>=min.pos]
	}
	### filter data
	##cat('filter data\n')
	if (is.filter.sample==TRUE){
		cat("remove outliers \n")
		for (i in 1:length(genomeF.native$features$ipd_pos))
			genomeF.native$features$ipd_pos[[i]] <- filter.outlier.byGenomeF(genomeF.native$features$ipd_pos[[i]])
        	for (i in 1:length(genomeF.native$features$ipd_neg))
			genomeF.native$features$ipd_neg[[i]] <- filter.outlier.byGenomeF(genomeF.native$features$ipd_neg[[i]])
		for (i in 1:length(genomeF.ctrl$features$ipd_pos))
        		genomeF.ctrl$features$ipd_pos[[i]] <- filter.outlier.byGenomeF(genomeF.ctrl$features$ipd_pos[[i]])
		for (i in 1:length(genomeF.ctrl$features$ipd_neg))
        		genomeF.ctrl$features$ipd_neg[[i]] <- filter.outlier.byGenomeF(genomeF.ctrl$features$ipd_neg[[i]])
	}

	## get result
	detection <- list()
	cat('detect modifiations\n')
	if (is.perm==TRUE) print('detect permuted data')
	for (i in 1:length(genomeF.native$features$ipd_pos) ){
		cat('detection modifications in forward strand : ', names(genomeF.native$features$ipd_pos)[i],'\n' )
		if (is.perm==TRUE){
			detection$pos[[i]] <- .Call('hieModel_CC_perm', genomeF.native$features$ipd_pos[[i]], genomeF.ctrl$features$ipd_pos[[i]],
                                as.integer(genomeF.native$genome.start.pos[[i]]), as.integer(genomeF.ctrl$genome.start.pos[[i]]) ,
                                genomeSeq$pos[[names(genomeF.native$features$ipd_pos)[i]]], context.effect, names(context.effect), 
				as.integer(0), as.integer(left.len), as.integer(right.len), method )
		}else{
			detection$pos[[i]] <- .Call('hieModel_CC', genomeF.native$features$ipd_pos[[i]], genomeF.ctrl$features$ipd_pos[[i]],
                        	as.integer(genomeF.native$genome.start.pos[[i]]), as.integer(genomeF.ctrl$genome.start.pos[[i]]) ,
                        	genomeSeq$pos[[names(genomeF.native$features$ipd_pos)[i]]], context.effect, names(context.effect), 
				as.integer(0), as.integer(left.len), as.integer(right.len), method )
		}
		names(detection$pos[[i]])[names(detection$pos[[i]])=='t_stat'] <- 't.stat'
		names(detection$pos[[i]])[names(detection$pos[[i]])=='ipd_ratio'] <- 'ipd.ratio'
	}
	names(detection$pos) <- names(genomeF.native$features$ipd_pos)

	for (i in 1:length(genomeF.native$features$ipd_neg) ){
                cat('detection modifications in backward strand : ', names(genomeF.native$features$ipd_neg)[i],'\n' )
		if (is.perm==TRUE){
			detection$neg[[i]] <- .Call('hieModel_CC_perm', genomeF.native$features$ipd_neg[[i]], genomeF.ctrl$features$ipd_neg[[i]],
                        	as.integer(genomeF.native$genome.start.neg[[i]]), as.integer(genomeF.ctrl$genome.start.neg[[i]]) ,
                         	genomeSeq$neg[[names(genomeF.native$features$ipd_neg)[i]]], context.effect, names(context.effect), as.integer(1), as.integer(left.len), as.integer(right.len), method )

		}else{
			detection$neg[[i]] <- .Call('hieModel_CC', genomeF.native$features$ipd_neg[[i]], genomeF.ctrl$features$ipd_neg[[i]],
                        	as.integer(genomeF.native$genome.start.neg[[i]]), as.integer(genomeF.ctrl$genome.start.neg[[i]]) ,
                         	genomeSeq$neg[[names(genomeF.native$features$ipd_neg)[i]]], context.effect, names(context.effect),
	 			as.integer(1), as.integer(left.len), as.integer(right.len), method )
		}
		names(detection$neg[[i]])[names(detection$neg[[i]])=='t_stat'] <- 't.stat'
                names(detection$neg[[i]])[names(detection$neg[[i]])=='ipd_ratio'] <- 'ipd.ratio'
	}	
	names(detection$neg) <- names(genomeF.native$features$ipd_neg)
	detection$genome.start.pos <- genomeF.native$genome.start.pos
	detection$genome.start.neg <- genomeF.native$genome.start.neg
	
		
	detection
	
}



