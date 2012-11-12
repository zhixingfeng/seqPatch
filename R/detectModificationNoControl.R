detectModification.NC.adaptive <- function(genomeF.native, genomeSeq, context.effect.set, is.merge=TRUE)
{
	method <- 'hieModel'
        # check consistency of context.effect.set
        context <- names(context.effect.set)[1:(length(context.effect.set)-2)]
        context.split <- strsplit(context, '_')
        context.left.len <- sapply(context.split, function(x) as.numeric(x[3]) )
        context.right.len <- sapply(context.split, function(x) as.numeric(x[4]) )
        if (any(context.left.len!=context.effect.set$left.len) | any(context.right.len!=context.effect.set$right.len))
                stop('inconsistent context and left.len rigth.len.')
        if (any(context.effect.set$left.len!= sort(context.effect.set$left.len,decreasing=TRUE)) |
                any(context.effect.set$right.len!= sort(context.effect.set$right.len,decreasing=TRUE)))
                stop('left.len and rigth.len should be decreasing.')
        if (length(context.effect.set)-2 <= 1) stop('length of context.effect.set should be > 3')

	# start to detect
        detection <- list()
        cur.context.effect <- context.effect.set[[1]]
        detection[[1]] <- detectModification.NC(genomeF.native, genomeSeq, cur.context.effect,
                        left.len = context.effect.set$left.len[1], right.len = context.effect.set$right.len[1], method )

	for (i in 2:(length(context.effect.set)-2) ){
                # mask genome locus that have been detected
                for (j in 1:length(genomeF.native$features$ipd_pos)){
                        cur.chr <- names(genomeF.native$features$ipd_pos[j])
                        if ( is.null(detection[[i-1]]$pos[[cur.chr]]) ) stop('inconsistent genomeF and detection.')
                        idx.mask <- which(detection[[i-1]]$pos[[cur.chr]]$is_findContext == 1)
                        for (k in idx.mask){
                                genomeF.native$features$ipd_pos[[j]][[k]] <- numeric(0)
                        }
                }

                for (j in 1:length(genomeF.native$features$ipd_neg)){
                        cur.chr <- names(genomeF.native$features$ipd_neg[j])
                        if ( is.null(detection[[i-1]]$neg[[cur.chr]]) ) stop('inconsistent genomeF and detection.')
                        idx.mask <- which(detection[[i-1]]$neg[[cur.chr]]$is_findContext == 1)
                        for (k in idx.mask){
                                genomeF.native$features$ipd_neg[[j]][[k]] <- numeric(0)
                        }
                }

                cur.context.effect <- context.effect.set[[i]]
                detection[[i]] <- detectModification.NC(genomeF.native, genomeSeq, cur.context.effect,
                                left.len = context.effect.set$left.len[i], right.len = context.effect.set$right.len[i], method )

        }

	# merge detection results
        if (is.merge == TRUE){
                detection.merged <- detection[[1]]
                for (i in 2:length(detection) ){
			# positive strand
                        for (j in 1:length(detection[[i]]$pos)){
                                cur.idx <- which(detection[[i]]$pos[[j]]$is_findContext == 1)
                                for (k in cur.idx){
                                        detection.merged$pos[[j]]$mu_native[k] <- detection[[i]]$pos[[j]]$mu_native[k];
                                        detection.merged$pos[[j]]$mu_ctrl[k] <- detection[[i]]$pos[[j]]$mu_ctrl[k];
                                        detection.merged$pos[[j]]$sigma2_native[k] <- detection[[i]]$pos[[j]]$sigma2_native[k];
                                        detection.merged$pos[[j]]$sigma2_ctrl[k] <- detection[[i]]$pos[[j]]$sigma2_ctrl[k];
                                        detection.merged$pos[[j]]$sampleSize_native[k] <- detection[[i]]$pos[[j]]$sampleSize_native[k];
                                        detection.merged$pos[[j]]$sampleSize_ctrl[k] <- detection[[i]]$pos[[j]]$sampleSize_ctrl[k];
                                        detection.merged$pos[[j]]$t.stat[k] <- detection[[i]]$pos[[j]]$t.stat[k];
                                        detection.merged$pos[[j]]$L_log_alt[k] <- detection[[i]]$pos[[j]]$L_log_alt[k];
                                        detection.merged$pos[[j]]$L_log_null[k] <- detection[[i]]$pos[[j]]$L_log_null[k];
                                        detection.merged$pos[[j]]$LR_log[k] <- detection[[i]]$pos[[j]]$LR_log[k];
                                        detection.merged$pos[[j]]$is_findContext[k] <- 1
                                        #print(k)
                                }
                        }
			# negative strand 
			for (j in 1:length(detection[[i]]$neg)){
                                cur.idx <- which(detection[[i]]$neg[[j]]$is_findContext == 1)
                                for (k in cur.idx){
                                        detection.merged$neg[[j]]$mu_native[k] <- detection[[i]]$neg[[j]]$mu_native[k];
                                        detection.merged$neg[[j]]$mu_ctrl[k] <- detection[[i]]$neg[[j]]$mu_ctrl[k];
                                        detection.merged$neg[[j]]$sigma2_native[k] <- detection[[i]]$neg[[j]]$sigma2_native[k];
                                        detection.merged$neg[[j]]$sigma2_ctrl[k] <- detection[[i]]$neg[[j]]$sigma2_ctrl[k];
                                        detection.merged$neg[[j]]$sampleSize_native[k] <- detection[[i]]$neg[[j]]$sampleSize_native[k];
                                        detection.merged$neg[[j]]$sampleSize_ctrl[k] <- detection[[i]]$neg[[j]]$sampleSize_ctrl[k];
                                        detection.merged$neg[[j]]$t.stat[k] <- detection[[i]]$neg[[j]]$t.stat[k];
                                        detection.merged$neg[[j]]$L_log_alt[k] <- detection[[i]]$neg[[j]]$L_log_alt[k];
                                        detection.merged$neg[[j]]$L_log_null[k] <- detection[[i]]$neg[[j]]$L_log_null[k];
                                        detection.merged$neg[[j]]$LR_log[k] <- detection[[i]]$neg[[j]]$LR_log[k];
                                        detection.merged$neg[[j]]$is_findContext[k] <- 1
                                        #print(k)
                                }
                        }

                }
                return(detection.merged)

        }else{
                return(detection)
        }
	

}


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

