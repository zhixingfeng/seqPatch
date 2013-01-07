detectModProp.NC <- function(genomeF.native, genomeSeq, context.hyperPara, left.len = 6, right.len = 1, max.iter=50)
{
	# check consistence of context.hyperPara and left.len, right.len
	if (nchar(names(context.hyperPara)[1]) != left.len + right.len + 1)
		stop ('left.len + right.len + 1 should be equal to length of context')	
	
	# detect forward strand 
	detection <- list()
	
	for (i in 1:length(genomeF.native$features$ipd_pos) ){
		cat('detection modifications proportion in forward strand : ', names(genomeF.native$features$ipd_pos)[i],'\n' )
		detection$pos[[i]] <- .Call('R_API_DetectModProp_NC', genomeF.native$features$ipd_pos[[i]], genomeF.native$features$moleculeID_pos[[i]],
				as.integer(genomeF.native$genome.start.pos[[i]]), genomeSeq$pos[[names(genomeF.native$features$ipd_pos)[i]]], 
				context.hyperPara, names(context.hyperPara), as.integer(0), as.integer(left.len), as.integer(right.len), as.integer(max.iter) )
	}
	names(detection$pos) <- names(genomeF.native$features$ipd_pos)	

	detection

}



