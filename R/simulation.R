sim.genomeF.C <- function(prop=c(rep(1,50),rep(0.5,50), rep(0,20000)), genomeF.wga, d=1.5, factor.n.ccs = 1, factor.n.mol = 1)
{
	if(length(prop)>length(genomeF.wga$features$ipd_pos[[1]]))
		stop('prop is too large')
	n.ccs.all <- .Call('R_API_getNumCCS',genomeF.wga$features$ipd_pos[[1]],
                                                genomeF.wga$features$moleculeID_pos[[1]], as.integer(1))
	n.mol.all <- .Call('R_API_getNumMol',genomeF.wga$features$ipd_pos[[1]],
                                                genomeF.wga$features$moleculeID_pos[[1]], as.integer(1))
	mu.0 <- sapply(genomeF.wga$features$ipd_pos[[1]], mean)
	sd.0 <- sapply(genomeF.wga$features$ipd_pos[[1]], sd)
	
	n.ccs.all <- round(n.ccs.all*factor.n.ccs)
	n.mol.all <- round( n.mol.all*factor.n.mol)
	n.ccs.all[n.ccs.all <= 0] <- 1
	n.mol.all[n.mol.all <= 0] <- 1
	rl <- .Call('R_API_getIPDandMolID', as.numeric(prop), as.numeric(n.ccs.all), as.numeric(n.mol.all),
			as.numeric( mu.0), as.numeric(sd.0), as.numeric(d))	

	rl
}


sim.genomeF <- function(prop=c(rep(1,50),rep(0.5,50), rep(0,20000)), n.mol.lambda = 30, n.ccs.lambda=5,
		mu.0 = 0, sigma.0 = 1, mu.1=1.5, sigma.1 = 1)
{
	genomeF.native <- list()
	genomeF.native$features$ipd_pos$chr <- list()
	genomeF.native$features$moleculeID_pos$chr <- list()
	genomeF.native$genome.start.pos$chr <- 1
	genomeF.native$genome.start.neg <- genomeF.native$genome.start.pos

	genomeF.wga <- list()
	genomeF.wga$features$ipd_pos$chr <- list()
	genomeF.wga$features$moleculeID_pos$chr <- list()
	genomeF.wga$genome.start.pos$chr <- 1
	genomeF.wga$genome.start.neg <- genomeF.native$genome.start.pos

	for (i in 1:length(prop)){
		if (i==1000*floor(i/1000)) cat('simulated', i, 'sites\r')
		# native sample
		n.mol <- rpois(1, n.mol.lambda)
		n.mol[n.mol == 0]  <- 1
		if (n.ccs.lambda == 1){ 
			n.ccs <- rep(1, n.mol)
		}else{
			n.ccs <- rpois(n.mol, n.ccs.lambda)
		}
		n.ccs[n.ccs == 0] <- 1	
		cur.idx <- mapply(function(x,idx,n.ccs) rep(idx,n.ccs) , idx=1:n.mol, n.ccs=n.ccs)
		ipd <- list()
		j.break <- round(n.mol * prop[i])
		for (j in 1:n.mol){
			if (j <= j.break){
				ipd[[j]] <- rnorm(n.ccs[j], mu.1, sigma.1)
			}else{
				ipd[[j]] <- rnorm(n.ccs[j], mu.0, sigma.0)
			}
		}
		genomeF.native$features$ipd_pos$chr[[i]] <- unlist(ipd)
		genomeF.native$features$moleculeID_pos$chr[[i]] <- as.numeric(unlist(cur.idx))
		
		# wga sample 
		n.mol <- rpois(1, n.mol.lambda)
                n.mol[n.mol == 0]  <- 1
                if (n.ccs.lambda == 1){
                        n.ccs <- rep(1, n.mol)
                }else{
                        n.ccs <- rpois(n.mol, n.ccs.lambda)
                }
                n.ccs[n.ccs == 0] <- 1
		cvg <- sum(n.ccs)
		genomeF.wga$features$ipd_pos$chr[[i]] <- rnorm(cvg, mu.0, sigma.0)
		genomeF.wga$features$moleculeID_pos$chr[[i]] <- as.numeric(1:cvg)
	}		
	cat('simulated', length(prop), 'sites\n')	
	genomeF.native$features$ipd_neg <- genomeF.native$features$ipd_pos
	genomeF.native$features$moleculeID_neg <- genomeF.native$features$moleculeID_pos
	genomeF.wga$features$ipd_neg <- genomeF.wga$features$ipd_pos
	genomeF.wga$features$moleculeID_neg <- genomeF.wga$features$moleculeID_pos
	
	rl <- list()
	rl$native <- genomeF.native
	rl$wga <- genomeF.wga
	rl		
}


