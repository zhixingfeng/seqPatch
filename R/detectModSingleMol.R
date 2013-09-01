#mergeZscoreSingleMolecule <- function(z.score)
#{
#        z.score.pool <- numeric(0)
#	idx.pool <- numeric(0)
        # forward strand
#        for (i in 1:length(z.score$pos)){
#                cur.z.score.pool <- .Call('R_API_unlist', z.score$pos[[i]])	
#		z.score.pool <- c(z.score.pool, cur.z.score.pool$z)
#		idx.pool <- c(idx.pool, cur.z.score.pool$idx)
#        }
        # backward strand
#        for (i in 1:length(z.score$neg)){
#                cur.z.score.pool <- .Call('R_API_unlist', z.score$neg[[i]])
#                z.score.pool <- c(z.score.pool, cur.z.score.pool$z)
#        	idx.pool <- c(idx.pool, cur.z.score.pool$idx)
#	}
#        z.score.pool
#}

mergeZscoreSingleMolecule <- function(z.score)
{
	z.score.pool <- numeric(0)
	# forward strand 
	for (i in 1:length(z.score$pos)){
		cur.z.score.pool <- unlist(z.score$pos[[i]])
		z.score.pool <- c(z.score.pool, cur.z.score.pool)
	}
	# backward strand 
	for (i in 1:length(z.score$neg)){
                cur.z.score.pool <- unlist(z.score$neg[[i]])
                z.score.pool <- c(z.score.pool, cur.z.score.pool)
        }
	z.score.pool	
}
getLocfdrSingleMolecule <- function(z.score, x, fdr)
{
	fdr.all <- list()
	fdr.all$pos <- list()
	fdr.all$neg <- list()
	
	# forward strand 
	cat('forward strand\n')
	for (i in 1:length(z.score$pos)){
		cur.ref <- names(z.score$pos[i])
		cat(cur.ref,'\n')
		fdr.all$pos[[cur.ref]] <- .Call('R_API_getLocfdrSingleMolecule', z.score$pos[[i]], x, fdr)
	}				

	# backward strand
	cat('backward strand\n')
        for (i in 1:length(z.score$neg)){
                cur.ref <- names(z.score$neg[i])
		cat(cur.ref,'\n')
                fdr.all$neg[[cur.ref]] <- .Call('R_API_getLocfdrSingleMolecule', z.score$neg[[i]], x, fdr)
        }
	fdr.all
}
getZscoreSingleMolecule.CC <- function(genomeF.native, genomeF.ctrl, is.merge=FALSE, is.z=TRUE, is.delta=FALSE)
{
	z.score <- list()
	z.score$pos <- list()
	z.score$neg <- list()
	if (any(names(genomeF.native$features$ipd_pos) != names(genomeF.native$features$ipd_neg) ) ){
		stop('incompatible genomeF.native')
	}
		
	# forward strand 
	cat('forward strand\n')
	for (i in 1:length(genomeF.native$features$ipd_pos)){
		cur.ref <- names(genomeF.native$features$ipd_pos[i])

		start.native <- genomeF.native$genome.start.pos[[cur.ref]]
		start.ctrl <- genomeF.ctrl$genome.start.pos[[cur.ref]]

		ipd.native <- genomeF.native$features$ipd_pos[[cur.ref]]
		ipd.ctrl <- genomeF.ctrl$features$ipd_pos[[cur.ref]]

		mol.id.native <- genomeF.native$features$moleculeID_pos[[cur.ref]]	
		mol.id.ctrl <- genomeF.ctrl$features$moleculeID_pos[[cur.ref]]
			
		cat(cur.ref, '\n')
		z.score$pos[[cur.ref]] <- .Call('R_API_getZscoreSingleMolecule_CC', ipd.native, ipd.ctrl,
			 mol.id.native, mol.id.ctrl, as.integer(start.native), as.integer(start.ctrl), as.integer(is.z),
			as.integer(is.delta))
	}		
	
	# backward strand
	cat('backward strand\n')
	for (i in 1:length(genomeF.native$features$ipd_neg)){
                cur.ref <- names(genomeF.native$features$ipd_neg[i])

                start.native <- genomeF.native$genome.start.neg[[cur.ref]]
                start.ctrl <- genomeF.ctrl$genome.start.neg[[cur.ref]]

                ipd.native <- genomeF.native$features$ipd_neg[[cur.ref]]
                ipd.ctrl <- genomeF.ctrl$features$ipd_neg[[cur.ref]]

                mol.id.native <- genomeF.native$features$moleculeID_neg[[cur.ref]]
                mol.id.ctrl <- genomeF.ctrl$features$moleculeID_neg[[cur.ref]]
		
		cat(cur.ref, '\n')
                z.score$neg[[cur.ref]] <- .Call('R_API_getZscoreSingleMolecule_CC', ipd.native, ipd.ctrl,
                         mol.id.native, mol.id.ctrl, as.integer(start.native), as.integer(start.ctrl), as.integer(is.z),
			as.integer(is.delta))
        }
	z.score$genome.start.pos <- genomeF.native$genome.start.pos
	z.score$genome.start.neg <- genomeF.native$genome.start.neg
	if (is.merge==TRUE)
		return(mergeZscoreSingleMolecule(z.score))
	else  
		return(z.score)
}

getZscoreSingleMolecule.CC.one.site <- function(x, x.idx, y, is.z=TRUE)
{
	x.group <- split(x,x.idx)
	if (is.z==TRUE){
		sapply(x.group, function(x,y) {t.stat <- t.test(x,y,var.equal=TRUE)$stat; t.df=t.test(x,y,var.equal=TRUE)$para;qnorm(pt(t.stat,t.df))}, y=y)
	}else{
		sapply(x.group, function(x,y) t.test(x,y,var.equal=TRUE)$stat, y=y)
	}
}




