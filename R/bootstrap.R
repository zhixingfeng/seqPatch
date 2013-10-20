get.bootstrap.FDR <- function(detected.prop, prop.null)
{
	prop.cutoff <- seq(0.99,0.01,-0.01)
	FDR <- rep(NaN,length(prop.cutoff))
	detected.prop.pool <- collapse.rl(detected.prop)
	
	N <- length(detected.prop.pool)
	for (i in 1:length(prop.cutoff)){
		num.detected.null <- sum(prop.null >= prop.cutoff[i]) * (N/length(prop.null))
		num.detected <- sum(detected.prop.pool >= prop.cutoff[i], na.rm=TRUE)
		FDR[i] <- num.detected.null / num.detected
	}
	list(FDR=FDR, prop.cutoff=prop.cutoff)
}


get.n.ccs.dist <- function(genomeF.native)
{
	mol.id.group <- list()
        mol.id.group$pos <- list()
        mol.id.group$neg <- list()

        for (i in 1:length(genomeF.native$features$moleculeID_pos)){
                cur.ref <- names(genomeF.native$features$moleculeID_pos[i])
                mol.id.group$pos[[cur.ref]] <- .Call('R_API_getNumCCS',genomeF.native$features$ipd_pos[[i]],
                                                genomeF.native$features$moleculeID_pos[[i]], as.integer(1))

        }

        for (i in 1:length(genomeF.native$features$moleculeID_neg)){
                cur.ref <- names(genomeF.native$features$moleculeID_neg[i])
                mol.id.group$neg[[cur.ref]] <- .Call('R_API_getNumCCS',genomeF.native$features$ipd_neg[[i]],
                                                genomeF.native$features$moleculeID_neg[[i]], as.integer(1))

        }
        unlist(mol.id.group)

}


get.n.mol.dist <- function(genomeF.native)
{
        mol.id.group <- list()
        mol.id.group$pos <- list()
        mol.id.group$neg <- list()

        for (i in 1:length(genomeF.native$features$moleculeID_pos)){
                cur.ref <- names(genomeF.native$features$moleculeID_pos[i])
                mol.id.group$pos[[cur.ref]] <- .Call('R_API_getNumMol',genomeF.native$features$ipd_pos[[i]],
                                                genomeF.native$features$moleculeID_pos[[i]], as.integer(1))

        }

        for (i in 1:length(genomeF.native$features$moleculeID_neg)){
                cur.ref <- names(genomeF.native$features$moleculeID_neg[i])
                mol.id.group$neg[[cur.ref]] <- .Call('R_API_getNumMol',genomeF.native$features$ipd_neg[[i]],
                                                genomeF.native$features$moleculeID_neg[[i]], as.integer(1))

        }
        unlist(mol.id.group)

}

bootstrap.genomeF <- function(genomeF.native, genomeF.wga, B = 1000000, n.CCS=NULL, n.mol = NULL)
{
	# get distribution of number of CCS
	if (is.null(n.CCS)){
		cat('get distribution of number of CCS\n')
		n.CCS <- get.n.ccs.dist(genomeF.native)
	}
	# get distribution of number of molecules	
	if (is.null(n.mol)){
		cat('get distribution of number of molecules\n')
		n.mol <- get.n.mol.dist(genomeF.native)
	}
	# get mean and standard variance 
	cat('get mu.0\n')
	mu.0 <- sapply(genomeF.wga$features$ipd_pos[[1]], mean)
	cat('get sigma.0\n')
	sd.0 <- sapply(genomeF.wga$features$ipd_pos[[1]], sd)

	# bootstrap
	genomeF.native.B <- list()
	genomeF.native.B$features$ipd_pos$chr <- list()
	genomeF.native.B$features$moleculeID_pos$chr <- list()
	genomeF.native.B$genome.start.pos$chr <- 1
	genomeF.native.B$genome.start.neg <- genomeF.native.B$genome.start.pos

	genomeF.wga.B <- list()
        genomeF.wga.B$features$ipd_pos$chr <- list()
        genomeF.wga.B$features$moleculeID_pos$chr <- list()
        genomeF.wga.B$genome.start.pos$chr <- 1
        genomeF.wga.B$genome.start.neg <- genomeF.wga.B$genome.start.pos

	n.mol.B <- sample(n.mol[n.mol>2], B, replace=TRUE)
	mu.0.B <- sample(mu.0[!is.na(mu.0)], B, replace=TRUE)
	sd.0.B <- sample(sd.0[!is.na(sd.0)], B, replace=TRUE)

	cat('boostrap\n')	
	for (i in 1:B){
		cur.n.ccs <- sample(n.CCS, n.mol.B[i], replace=TRUE)
		cur.idx <- mapply(function(idx,n.ccs) rep(idx,n.ccs) , idx=1:n.mol.B[i], n.ccs=cur.n.ccs)	
		cur.idx.pool <- unlist(cur.idx)
		
		genomeF.native.B$features$ipd_pos$chr[[i]] <- rnorm(length(cur.idx.pool), mu.0.B[i], sd.0.B[i])
		genomeF.native.B$features$moleculeID_pos$chr[[i]] <- as.numeric(cur.idx.pool)
		
		genomeF.wga.B$features$ipd_pos$chr[[i]] <- rnorm(length(cur.idx.pool), mu.0.B[i], sd.0.B[i])
		genomeF.wga.B$features$moleculeID_pos$chr[[i]] <- as.numeric(cur.idx.pool)
		
		if (i == 1000*floor(i/1000)) cat(i,'\r')
	}
	cat(B,'\n')
	genomeF.native.B$features$ipd_neg <- genomeF.native.B$features$ipd_pos
	genomeF.native.B$features$moleculeID_neg <- genomeF.native.B$features$moleculeID_pos
	genomeF.wga.B$features$ipd_neg <- genomeF.wga.B$features$ipd_pos
	genomeF.wga.B$features$moleculeID_neg <- genomeF.wga.B$features$moleculeID_pos	
	
	rl <- list()
	rl$native <- genomeF.native.B
	rl$wga <- genomeF.wga.B
	rl	
}

np.bootstrap.genomeF <- function(genomeF.native, genomeF.wga, min.cvg=50, n.CCS=NULL, n.mol = NULL)
{
	# get distribution of number of CCS
        if (is.null(n.CCS)){
                cat('get distribution of number of CCS\n')
                n.CCS <- get.n.ccs.dist(genomeF.native)
        }
        # get distribution of number of molecules
        if (is.null(n.mol)){
                cat('get distribution of number of molecules\n')
                n.mol <- get.n.mol.dist(genomeF.native)
        }
		
	cvg.wga <- sapply(genomeF.wga$features$ipd_pos[[1]], length)
	idx.cvg <- which(cvg.wga>=min.cvg)
	ipd.wga.cvg <- genomeF.wga$features$ipd_pos[[1]][idx.cvg]

	n.mol.B <- sample(n.mol[n.mol>2], length(ipd.wga.cvg), replace=TRUE)

	genomeF.native.B <- list()
        genomeF.native.B$features$ipd_pos$chr <- list()
        genomeF.native.B$features$moleculeID_pos$chr <- list()
        genomeF.native.B$genome.start.pos$chr <- 1
        genomeF.native.B$genome.start.neg <- genomeF.native.B$genome.start.pos

        genomeF.wga.B <- list()
        genomeF.wga.B$features$ipd_pos$chr <- list()
        genomeF.wga.B$features$moleculeID_pos$chr <- list()
        genomeF.wga.B$genome.start.pos$chr <- 1
        genomeF.wga.B$genome.start.neg <- genomeF.wga.B$genome.start.pos

	
	for (i in 1:length(ipd.wga.cvg)){
		cur.n.ccs <- sample(n.CCS, n.mol.B[i], replace=TRUE)
                cur.idx <- mapply(function(idx,n.ccs) rep(idx,n.ccs) , idx=1:n.mol.B[i], n.ccs=cur.n.ccs)
                cur.idx.pool <- unlist(cur.idx)
	
		genomeF.native.B$features$ipd_pos$chr[[i]] <- sample(ipd.wga.cvg[[i]], length(cur.idx.pool) , replace=TRUE)
                genomeF.native.B$features$moleculeID_pos$chr[[i]] <- as.numeric(cur.idx.pool)

                genomeF.wga.B$features$ipd_pos$chr[[i]] <- sample(ipd.wga.cvg[[i]], replace=TRUE)
                genomeF.wga.B$features$moleculeID_pos$chr[[i]] <- 1:length(ipd.wga.cvg[[i]])

		if (i == 1000*floor(i/1000)) cat(i,'\r')
	}
	cat(length(ipd.wga.cvg),'\n')
        genomeF.native.B$features$ipd_neg <- genomeF.native.B$features$ipd_pos
        genomeF.native.B$features$moleculeID_neg <- genomeF.native.B$features$moleculeID_pos
        genomeF.wga.B$features$ipd_neg <- genomeF.wga.B$features$ipd_pos
        genomeF.wga.B$features$moleculeID_neg <- genomeF.wga.B$features$moleculeID_pos

	rl <- list()
        rl$native <- genomeF.native.B
        rl$wga <- genomeF.wga.B
        rl

}






