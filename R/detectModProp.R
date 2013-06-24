detectModProp.CC <- function(genomeF.native, genomeF.ctrl, max.iter=5000)
{
        # detect forward strand
        detection <- list()

        for (i in 1:length(genomeF.native$features$ipd_pos) ){
                cat('detection modifications proportion in forward strand : ', names(genomeF.native$features$ipd_pos)[i],'\n' )
                
		detection$pos[[i]] <- .Call('R_API_DetectModProp_CC', genomeF.native$features$ipd_pos[[i]], genomeF.native$features$moleculeID_pos[[i]],
				genomeF.ctrl$features$ipd_pos[[i]], as.integer(genomeF.native$genome.start.pos[[i]]), 
				as.integer(genomeF.ctrl$genome.start.pos[[i]]), as.integer(max.iter) )
        }
        names(detection$pos) <- names(genomeF.native$features$ipd_pos)

        # detect backward strand
        for (i in 1:length(genomeF.native$features$ipd_neg) ){
                cat('detection modifications proportion in backward strand : ', names(genomeF.native$features$ipd_neg)[i],'\n' )
		
		detection$neg[[i]] <- .Call('R_API_DetectModProp_CC', genomeF.native$features$ipd_neg[[i]], genomeF.native$features$moleculeID_neg[[i]],
                                genomeF.ctrl$features$ipd_neg[[i]], as.integer(genomeF.native$genome.start.neg[[i]]),
                                as.integer(genomeF.ctrl$genome.start.neg[[i]]), as.integer(max.iter) )

        }
        names(detection$neg) <- names(genomeF.native$features$ipd_neg)

        detection$genome.start.pos <- genomeF.native$genome.start.pos
        detection$genome.start.neg <- genomeF.native$genome.start.neg

        detection

}


detectModProp.NC <- function(genomeF.native, genomeSeq, context.hyperPara, left.len = 6, right.len = 1, max.iter=5000)
{
	# check consistence of context.hyperPara and left.len, right.len
	if (nchar(names(context.hyperPara)[1]) != left.len + right.len + 1)
		stop ('left.len + right.len + 1 should be equal to length of context')	
	
	# detect forward strand 
	detection <- list()
	
	for (i in 1:length(genomeF.native$features$ipd_pos) ){
		cat('detection modifications proportion in forward strand : ', names(genomeF.native$features$ipd_pos)[i],'\n' )
		if (is.null(genomeSeq$pos[[names(genomeF.native$features$ipd_pos)[i]]] )){
                        cat('Warning: can NOT find names',(genomeF.native$features$ipd_pos)[i],'in genomeSeq.\n')
                        continue;
                }

		detection$pos[[i]] <- .Call('R_API_DetectModProp_NC', genomeF.native$features$ipd_pos[[i]], genomeF.native$features$moleculeID_pos[[i]],
				as.integer(genomeF.native$genome.start.pos[[i]]), genomeSeq$pos[[names(genomeF.native$features$ipd_pos)[i]]], 
				context.hyperPara, names(context.hyperPara), as.integer(0), as.integer(left.len), as.integer(right.len), as.integer(max.iter) )
	}
	names(detection$pos) <- names(genomeF.native$features$ipd_pos)	

	# detect backward strand
        for (i in 1:length(genomeF.native$features$ipd_neg) ){
                cat('detection modifications proportion in backward strand : ', names(genomeF.native$features$ipd_neg)[i],'\n' )
		if (is.null(genomeSeq$neg[[names(genomeF.native$features$ipd_neg)[i]]] )){
			cat('Warning: can NOT find names',(genomeF.native$features$ipd_neg)[i],'in genomeSeq.\n')
			continue; 
		}
                detection$neg[[i]] <- .Call('R_API_DetectModProp_NC', genomeF.native$features$ipd_neg[[i]], genomeF.native$features$moleculeID_neg[[i]],
                                as.integer(genomeF.native$genome.start.neg[[i]]), genomeSeq$neg[[names(genomeF.native$features$ipd_neg)[i]]],
                                context.hyperPara, names(context.hyperPara), as.integer(1), as.integer(left.len), as.integer(right.len), as.integer(max.iter) )
        }
        names(detection$neg) <- names(genomeF.native$features$ipd_neg)

	detection$genome.start.pos <- genomeF.native$genome.start.pos
        detection$genome.start.neg <- genomeF.native$genome.start.neg

	detection

}



