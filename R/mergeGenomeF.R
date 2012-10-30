mergeGenomeF <- function(genomeF.list)
{
	### check whether all element in genomeF.list have moleculeID or all not have 
	is.moleculeID <- !is.null(genomeF.list[[1]]$feature[['moleculeID_pos']])
	for (i in 1:length(genomeF.list)){
		if (is.moleculeID != !is.null(genomeF.list[[i]]$feature[['moleculeID_pos']]) )
			stop('Elements in genomeF.list must all have moleculeID or all have not')
	}

	genomeF <- genomeF.list[[1]]
	
	### get union of reference names
	ref.names <- character()	
	for (i in 1:length(genomeF.list) ){
		ref.names <- union(ref.names, names(genomeF.list[[i]]$features$ipd_pos))
	}	

	### merge genome.start and genome.end
	genome.start.pos.byRef <- list(); genome.start.neg.byRef <- list();	
	genome.end.pos.byRef <- list(); genome.end.neg.byRef <- list();
	for (i in 1:length(genomeF.list)){
		# forward strand
		for (j in 1:length(genomeF.list[[i]]$genome.start.pos)){
			cur.ref <- names(genomeF.list[[i]]$genome.start.pos[j])
			genome.start.pos.byRef[[cur.ref]] <- c(genome.start.pos.byRef[[cur.ref]], genomeF.list[[i]]$genome.start.pos[[j]])
			genome.end.pos.byRef[[cur.ref]] <- c(genome.end.pos.byRef[[cur.ref]], genomeF.list[[i]]$genome.start.pos[[j]] + 
								length(genomeF.list[[i]]$features$ipd_pos[[j]])-1)
		}
		
		# reverse strand
		for (j in 1:length(genomeF.list[[i]]$genome.start.neg)){
                        cur.ref <- names(genomeF.list[[i]]$genome.start.neg[j])
                        genome.start.neg.byRef[[cur.ref]] <- c(genome.start.neg.byRef[[cur.ref]], genomeF.list[[i]]$genome.start.neg[[j]])
                	genome.end.neg.byRef[[cur.ref]] <- c(genome.end.neg.byRef[[cur.ref]], genomeF.list[[i]]$genome.start.neg[[j]] +
								 length(genomeF.list[[i]]$features$ipd_neg[[j]])-1)
		}
	}
	
	for (i in 1:length(genome.start.pos.byRef)){
		genomeF$genome.start.pos[[names(genome.start.pos.byRef[i])]] <- as.integer(min(genome.start.pos.byRef[[i]]))
	}

	for (i in 1:length(genome.start.neg.byRef)){
                genomeF$genome.start.neg[[names(genome.start.neg.byRef[i])]] <- as.integer(min(genome.start.neg.byRef[[i]]))
        }

	genome.end.pos <- list()
	genome.end.neg <- list()
	for (i in 1:length(genome.end.pos.byRef)){
		genome.end.pos[[names(genome.end.pos.byRef[i])]] <- as.integer(max(genome.end.pos.byRef[[i]]))
	}

	for (i in 1:length(genome.end.neg.byRef)){
                genome.end.neg[[names(genome.end.neg.byRef[i])]] <- as.integer(max(genome.end.neg.byRef[[i]]))
        }


	### merge features  
	for (i in 1:length(ref.names)){
		 cat('merge',ref.names[i],'\n')
		# merge ipd_pos
		data.list <- list()
		genome.start.list <- integer(0)
		idx <- 1
		for (j in 1:length(genomeF.list)){
			if (!is.null(genomeF.list[[j]]$features$ipd_pos[[ref.names[i]]])){
				data.list[[idx]] <- genomeF.list[[j]]$features$ipd_pos[[ref.names[i]]]
				genome.start.list[idx] <- genomeF.list[[j]]$genome.start.pos[[ref.names[i]]]
				idx <- idx + 1
			}
		}
		# cat('merge',ref.names[i],': genome.start.list is ',genome.start.list,'\n')
		
		cur.genome.start <- genomeF$genome.start.pos[[ref.names[i]]]
		cur.genome.end <- genome.end.pos[[ref.names[i]]]
		if (length(data.list)==1){
			genomeF$features$ipd_pos[[ref.names[i]]] <- data.list[[1]]
		}else{
			genomeF$features$ipd_pos[[ref.names[i]]] <- mergeGenomeF.core(data.list, genome.start.list, cur.genome.start, cur.genome.end) 
		}
		# merge ipd_neg
                data.list <- list()
                genome.start.list <- integer(0)
                idx <- 1
                for (j in 1:length(genomeF.list)){
                        if (!is.null(genomeF.list[[j]]$features$ipd_neg[[ref.names[i]]])){
                                data.list[[idx]] <- genomeF.list[[j]]$features$ipd_neg[[ref.names[i]]]
                                genome.start.list[idx] <- genomeF.list[[j]]$genome.start.neg[[ref.names[i]]]
                                idx <- idx + 1
                        }
                }
                cur.genome.start <- genomeF$genome.start.neg[[ref.names[i]]]
                cur.genome.end <- genome.end.neg[[ref.names[i]]]
		if(length(data.list)==1){
			genomeF$features$ipd_neg[[ref.names[i]]] <- data.list[[1]]
		}else{
                	genomeF$features$ipd_neg[[ref.names[i]]] <- mergeGenomeF.core(data.list, genome.start.list, cur.genome.start, cur.genome.end)
		}
		if (is.moleculeID==TRUE){
			# merge ipd_pos
	                data.list <- list()
        	        genome.start.list <- integer(0)
                	idx <- 1
                	for (j in 1:length(genomeF.list)){
                        	if (!is.null(genomeF.list[[j]]$features$moleculeID_pos[[ref.names[i]]])){
                                	data.list[[idx]] <- genomeF.list[[j]]$features$moleculeID_pos[[ref.names[i]]]
                                	genome.start.list[idx] <- genomeF.list[[j]]$genome.start.pos[[ref.names[i]]]
                                	idx <- idx + 1
                        	}
                	}
                	cur.genome.start <- genomeF$genome.start.pos[[ref.names[i]]]
                	cur.genome.end <- genome.end.pos[[ref.names[i]]]
                	if (length(data.list)==1){
				genomeF$features$moleculeID_pos[[ref.names[i]]] <- data.list[[1]]
			}else{
				genomeF$features$moleculeID_pos[[ref.names[i]]] <- mergeGenomeF.core(data.list, genome.start.list, cur.genome.start, cur.genome.end)
			}
			# merge ipd_neg
                        data.list <- list()
                        genome.start.list <- integer(0)
                        idx <- 1
                        for (j in 1:length(genomeF.list)){
                                if (!is.null(genomeF.list[[j]]$features$moleculeID_neg[[ref.names[i]]])){
                                        data.list[[idx]] <- genomeF.list[[j]]$features$moleculeID_neg[[ref.names[i]]]
                                        genome.start.list[idx] <- genomeF.list[[j]]$genome.start.neg[[ref.names[i]]]
                                        idx <- idx + 1
                                }
                        }
                        cur.genome.start <- genomeF$genome.start.neg[[ref.names[i]]]
                        cur.genome.end <- genome.end.neg[[ref.names[i]]]
                        if (length(data.list)==1){
				genomeF$features$moleculeID_neg[[ref.names[i]]] <- data.list[[1]]
			}else{
				genomeF$features$moleculeID_neg[[ref.names[i]]] <- mergeGenomeF.core(data.list, genome.start.list, cur.genome.start, cur.genome.end)	
			}
		}
	}
	genomeF	
	#data.list <- list()
	#genome.start.list <- integer(0)
	#data.list[[1]] <- genomeF.1$features$ipd_pos$ctg7180000000111
	#data.list[[2]] <- genomeF.2$features$ipd_pos$ctg7180000000111		
	#genome.start.list[1] <- genomeF.1$genome.start.pos$ctg7180000000111
	#genome.start.list[2] <- genomeF.2$genome.start.pos$ctg7180000000111
	#cur.genome.start <- genomeF$genome.start.pos$ctg7180000000111
	#cur.genome.end <- genome.end.pos$ctg7180000000111
}	

mergeGenomeF.core <- function(data.list, genome.start.list, cur.genome.start, cur.genome.end)
{
	if (min(genome.start.list) < cur.genome.start)
		stop('genome.start should be smaller than any elements in genome.start.list.\n')	
	if (length(data.list)!=length(genome.start.list))
		stop('length of data.list should be the same with genome.start.list.\n');
	.Call('mergeGenomeF', data.list, genome.start.list, cur.genome.start, cur.genome.end)
		
}



