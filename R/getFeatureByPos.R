getFeaturesAlongGenome <- function(alnsF, alnsIdx, is.use.CCS = FALSE, is.filter.outlier=FALSE, is.bc=FALSE)
{
	if (length(alnsF) != nrow(alnsIdx)) stop('length of alnsF and alnsIdx should be the same!\n')
	refseq.names <- levels(as.factor(alnsIdx$fullRefName))
	genome.start.pos <- list()
        genome.start.neg <- list()

	for (i in 1:length(refseq.names)){
		genome.start.pos[i] <- min(alnsIdx$tStart[alnsIdx$alignedStrand==0&alnsIdx$fullRefName==refseq.names[i]])
	        genome.start.neg[i] <- min(alnsIdx$tStart[alnsIdx$alignedStrand==1&alnsIdx$fullRefName==refseq.names[i]])
	}
	names(genome.start.pos) <- refseq.names
	names(genome.start.neg) <- refseq.names

	result.feature <- .Call('getFeaturesAlongGenome',alnsF, alnsIdx,names(alnsF[[1]]), names(alnsIdx), as.integer(is.use.CCS))

	result <- list()
	result$genome.start.pos <- genome.start.pos
	result$genome.start.neg <- genome.start.neg	
	result$features <- result.feature
	
	if (any(names(result$genome.start.pos) != names(result$features$ipd_pos)) | any(names(result$genome.start.neg) != names(result$features$ipd_neg)))
		stop('unmatched genome.start and features names.\n')
	
	if (is.use.CCS==TRUE){
		if (any(names(result$genome.start.pos) != names(result$features$moleculeID_pos)) | any(names(result$genome.start.neg) != names(result$features$moleculeID_neg)))
                	stop('unmatched genome.start and moleculeID names.\n')
	}
	if (is.bc==TRUE){
		cat('logarithm transform IPDs\n')
		for (i in 1:length(result$features$ipd_pos))
        		result$features$ipd_pos[[i]] <- sapply(result$features$ipd_pos[[i]],function(x) log(x+0.01))
		for (i in 1:length(result$features$ipd_neg))
        		result$features$ipd_neg[[i]] <- sapply(result$features$ipd_neg[[i]],function(x) log(x+0.01))
	}

	if (is.filter.outlier==TRUE){
		cat('filter outliers\n')
		for (i in 1:length(result$features$ipd_pos))
        		result$features$ipd_pos[[i]] <- filter.outlier.byGenomeF(result$features$ipd_pos[[i]])
		for (i in 1:length(result$features$ipd_neg))
       			result$features$ipd_neg[[i]] <- filter.outlier.byGenomeF(result$features$ipd_neg[[i]])
	}
	result
}






