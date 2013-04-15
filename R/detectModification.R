detectModification.adaptive <- function(genomeF.native, genomeF.ctrl, genomeSeq, context.effect.set, is.merge=TRUE)
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
	detection[[1]] <- detectModification(genomeF.native, genomeF.ctrl, genomeSeq, cur.context.effect,
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
		detection[[i]] <- detectModification(genomeF.native, genomeF.ctrl, genomeSeq, cur.context.effect,
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


detectModification <- function(genomeF.native, genomeF.ctrl, genomeSeq=NULL, context.effect=NULL, left.len=6, right.len=1, method='hieModel')
{
	### check method 
	if (method!='hieModel' & method != 'CC')
                stop('method should be \'hieModel\' or  \'CC\'')
	

	### check reference genome
	if (!is.null(genomeSeq)){
		if (names(genomeF.native$features$ipd_pos) != names(genomeSeq$pos) | names(genomeF.ctrl$features$ipd_pos) != names(genomeSeq$pos) | 
			names(genomeF.native$features$ipd_neg) != names(genomeSeq$neg) | names(genomeF.ctrl$features$ipd_neg) != names(genomeSeq$neg))
			stop('Wrong reference genome, check tag of the fasta file, i.e. string after \'>\'\n')
	}
	
	# check where genomeSeq provide for hierarchical model 
	if (method=='hieModel' & is.null(genomeSeq)){
		stop('genomeSeq could not be NULL when method=\'hieModel\' ')
	}
	
	
	detection <- list()
	if (method == 'CC'){
		left.len <- 0
		right.len <- 0	
		cat('detect DNA modification by case-control method\n')
		# generate a fake genomeSeq and context.effect to avoid change the C code
		fake.genomeSeq <- list()
		fake.genomeSeq$pos <- list()
		for (i in 1:length(genomeF.native$genome.start.pos)){
			len <- max(length(genomeF.native$feature$ipd_pos[[i]]) + genomeF.native$genome.start.pos[[i]]-1, 
					length(genomeF.native$feature$ipd_neg[[i]]) + genomeF.native$genome.start.neg[[i]]-1)
			fake.genomeSeq$pos[[i]] <- paste(rep('N',len),collapse='')
		}
		names(fake.genomeSeq$pos) <- names(genomeF.native$genome.start.pos)
		fake.genomeSeq$neg <- fake.genomeSeq$pos

		# OK, let me also make a fake context.effect 
		fake.context.effect<-list()
		fake.context <- paste(rep('N',left.len + right.len + 1), collapse = '')
		fake.context.effect[[fake.context]] <- list();
		fake.context.effect[[fake.context]]$mean <- 0
		fake.context.effect[[fake.context]]$var <- 0
		fake.context.effect[[fake.context]]$len <- 0
		detection <- hieModelEB(genomeF.native, genomeF.ctrl, fake.genomeSeq, fake.context.effect, left.len, right.len, method, is.perm=FALSE )
	}	
	if (method == 'hieModel' & !is.null(context.effect)){
		cat('detect DNA modification by hierarchical model (using historical data)\n')
		detection <- hieModelEB(genomeF.native, genomeF.ctrl, genomeSeq, context.effect, left.len, right.len, method, is.perm=FALSE )
	}
	
	if (method == 'hieModel' & is.null(context.effect)){
		cat('detect DNA modification by hierarchical model (no historical data)\n')
		detection <- hieModelEB.NH(genomeF.native, genomeF.ctrl, genomeSeq, left.len, right.len,
					FALSE, min.cvg, min.pos)
	}
	
	if (method=='CC') detection <- getZvalue.from.t.stat(detection, TRUE)

	detection
}

