hieModelEB.NH <- function(genomeF.native, genomeF.ctrl, genomeSeq, left.len=6, right.len=1,
		 is.filter.sample=FALSE,  min.cvg=5, min.pos=4)
{
	### filter data
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
	for (i in 1:length(genomeF.native$features$ipd_pos) ){
		cat('detection modifications in forward strand : ', names(genomeF.native$features$ipd_pos)[i],'\n' )
		detection$pos[[i]] <- .Call('hieModel_NH', genomeF.native$features$ipd_pos[[i]], genomeF.ctrl$features$ipd_pos[[i]],
                        	as.integer(genomeF.native$genome.start.pos[[i]]), as.integer(genomeF.ctrl$genome.start.pos[[i]]) ,
                        	genomeSeq$pos[[names(genomeF.native$features$ipd_pos)[i]]], as.integer(0), as.integer(left.len), as.integer(right.len),
				as.integer(min.pos), as.integer(min.cvg))
		names(detection$pos[[i]])[names(detection$pos[[i]])=='t_stat'] <- 't.stat'
		names(detection$pos[[i]])[names(detection$pos[[i]])=='ipd_ratio'] <- 'ipd.ratio'
	}
	names(detection$pos) <- names(genomeF.native$features$ipd_pos)
	for (i in 1:length(genomeF.native$features$ipd_neg) ){
                cat('detection modifications in backward strand : ', names(genomeF.native$features$ipd_neg)[i],'\n' )
		detection$neg[[i]] <- .Call('hieModel_NH', genomeF.native$features$ipd_neg[[i]], genomeF.ctrl$features$ipd_neg[[i]],
                        	as.integer(genomeF.native$genome.start.neg[[i]]), as.integer(genomeF.ctrl$genome.start.neg[[i]]) ,
                         	genomeSeq$neg[[names(genomeF.native$features$ipd_neg)[i]]], as.integer(1), as.integer(left.len), as.integer(right.len), 
				as.integer(min.pos), as.integer(min.cvg))
		names(detection$neg[[i]])[names(detection$neg[[i]])=='t_stat'] <- 't.stat'
                names(detection$neg[[i]])[names(detection$neg[[i]])=='ipd_ratio'] <- 'ipd.ratio'
	}	
	names(detection$neg) <- names(genomeF.native$features$ipd_neg)
	
	detection$genome.start.pos <- genomeF.native$genome.start.pos
	detection$genome.start.neg <- genomeF.native$genome.start.neg
		
	detection
	
}



