hieModelEB <- function(genomeF.native, genomeF.ctrl=NULL, genomeSeq, context.effect=NULL, left.len=6, right.len=1, method='hieModel', is.perm=FALSE)
{
	if (method!='hieModel.nc' & method!='hieModel' & method != 'CC' &  method !='pooling')
		stop('method should be \'hieModel.nc\', \'hieModel\', \'CC\' or \'pooling\'')
	if (method!='hieModel.nc' & is.null(genomeF.ctrl))
		stop('genomeF.ctrl should be specified')
	if (method=='hieModel.nc' & !is.null(genomeF.ctrl))
		warning('genomeF.ctrl would not be used when method=\'hieModel.nc\'')
	
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
                         	genomeSeq$neg[[names(genomeF.native$features$ipd_neg)[i]]], context.effect, names(context.effect),
				as.integer(1), as.integer(left.len), as.integer(right.len), method )

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


hieModel_core <- function(IPD_mean, IPD_var, IPD_n, max.iter=10)
{
	# call hieModel_core defined in seqPatch2.cpp
	.Call('R_API_hieModel_core', IPD_mean, IPD_var, IPD_n, as.integer(max.iter) )
}



