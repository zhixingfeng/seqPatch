detectModPropEB <- function(genomeF.native, genomeF.wga, f1.x, f1.y, max.iter = 100)
{
	if (is.null(genomeF.native$features$moleculeID_pos) | is.null(genomeF.native$features$moleculeID_neg))
		stop('can not find molecule ID in genomeF.native')
	
	rl <- list()
	rl$pos <- list()
	rl$neg <- list()
	
	# forward strand
	cat('detect forward strand.\n')
	for (i in 1:length(genomeF.native$features$ipd_pos)){
		cur.ref <- names(genomeF.native$features$ipd_pos[i])
		if (is.null(genomeF.wga$features$ipd_pos[[cur.ref]]))	
			stop('can not find \'',cur.ref, '\' in genomeF.wga')
		ipd.wga <- genomeF.wga$features$ipd_pos[[cur.ref]]
		rl$pos[[cur.ref]] <- .Call('R_API_detectModProp_EB',genomeF.native$features$ipd_pos[[i]],
					genomeF.native$features$moleculeID_pos[[i]], ipd.wga, 
					as.integer(genomeF.native$genome.start.pos[[i]]),
					as.integer(genomeF.wga$genome.start.pos[[cur.ref]]),
					as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))	
			
	}

	# backward strand 
	cat('detect backward strand.\n')
	for (i in 1:length(genomeF.native$features$ipd_neg)){
                cur.ref <- names(genomeF.native$features$ipd_neg[i])
                if (is.null(genomeF.wga$features$ipd_neg[[cur.ref]]))
                        stop('can not find \'',cur.ref, '\' in genomeF.wga')
                ipd.wga <- genomeF.wga$features$ipd_neg[[cur.ref]]
                rl$neg[[cur.ref]] <- .Call('R_API_detectModProp_EB',genomeF.native$features$ipd_neg[[i]],
                                        genomeF.native$features$moleculeID_neg[[i]], ipd.wga,
                                        as.integer(genomeF.native$genome.start.neg[[i]]),
                                        as.integer(genomeF.wga$genome.start.neg[[cur.ref]]),
                                        as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))

        }	

	rl
}



EB.mixture.model <- function(IPD, idx, mu.0 = 0, sigma.0 = 1, f1.x, f1.y, max.iter = 100)
{
	s <- diff(f1.x)[1]
	f1.int <- sum(f1.y*s)
	if (abs(f1.int - 1)>=1e-6) stop('f1 is not a valid density function')

	.Call('R_API_EBmixture', as.numeric(IPD), as.numeric(idx), as.numeric(mu.0), as.numeric(sigma.0),
		 as.numeric(f1.x), as.numeric(f1.y),  as.integer(max.iter))
}


detectModPropEB.pooling <- function(z.score, mu.0 = 0, sigma.0 = 1, f1.x, f1.y, max.iter = 100)
{
	rl <- list()
	rl$pos <- list()
	rl$neg <- list()
	# forward strand
	cat('forward strand:\n')
	for (i in 1:length(z.score$pos)){
		cur.ref <- names(z.score$pos[i])
		cat(cur.ref,':\n', sep='')
		rl$pos[[cur.ref]] <- .Call('R_API_detectModProp_EB_pooling', z.score$pos[[i]], as.numeric(mu.0), as.numeric(sigma.0),
			                as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))

	}

	# backward strand
	cat('backward strand:\n')
	for (i in 1:length(z.score$neg)){
                cur.ref <- names(z.score$neg[i])
		cat(cur.ref,':\n', sep='')
                rl$neg[[cur.ref]] <- .Call('R_API_detectModProp_EB_pooling', z.score$neg[[i]], as.numeric(mu.0), as.numeric(sigma.0),
                                        as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))

        }
	
	rl$genome.start.pos <- z.score$genome.start.pos
	rl$genome.start.neg <- z.score$genome.start.neg
	rl

}

EB.mixture.model.pooling <- function(z.score, mu.0 = 0, sigma.0 = 1, f1.x, f1.y, max.iter = 100)
{
	.Call('R_API_EBmixture_pooling', as.numeric(z.score), as.numeric(mu.0), as.numeric(sigma.0),
		as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))
}



estimate.f1 <- function(genomeF.native, genomeF.wga, nulltype = 1, locfdr.df = 14, locfdr.pct0 = 1/4, locfdr.bre = 120, 
			f1.df = 14, f1.bre = 120, out.dir = NULL)
{
	if (!is.null(out.dir)){
		if (!file.exists(out.dir))	
			dir.create(out.dir, recursive=TRUE)
	}
	z.score.pool <- getZscoreSingleMolecule.CC(genomeF.native, genomeF.wga, is.z=TRUE, is.merge=TRUE)
	d.pool <- getZscoreSingleMolecule.CC(genomeF.native, genomeF.wga, is.delta=TRUE, is.merge=TRUE)
		
	idx.sel <- which(!is.na(z.score.pool))
	
	# estimate local FDR
	zz <- z.score.pool[idx.sel]
	if (!is.null(out.dir)){
		file <- paste(out.dir,'/','locfdr_report.pdf', sep='')
		pdf(file)
		rl.fdr <- locfdr(zz, nulltype=nulltype, df=locfdr.df, pct0 = locfdr.pct0, bre=locfdr.bre, mlests = c(0,1), plot=2)
		dev.off()
	}else{
		rl.fdr <- locfdr(zz, nulltype=nulltype, df=locfdr.df, pct0 = locfdr.pct0, bre=locfdr.bre, mlests = c(0,1), plot=0)
	}	

	# get weighted histogram of f1
	w <- 1 - rl.fdr$fdr
	w[w<=1e-6 | zz<0] <- 0
	idx <- idx.sel[w>1e-6]
	
	breaks <- seq(min(d.pool[idx]), max(d.pool[idx]), length=f1.bre)
	if (!is.null(out.dir)){
                file <- paste(out.dir,'/','f1_d_hist.pdf', sep='')
                pdf(file)
			f1.hist <- wtd.hist(d.pool[idx], weight=w[idx], breaks=breaks)
                dev.off()
        }else{
		f1.hist <- wtd.hist(d.pool[idx], weight=w[idx], breaks=breaks, plot=FALSE)
        }
	
	# fit f1 
	data <- as.data.frame(cbind(f1.hist$counts, f1.hist$mids))
	names(data) <- c('y', 'x')
	f1 <- glm(y~ns(x,df=f1.df), quasipoisson, data=data)		

	s <- diff(breaks)[1]	
	f1$fitted.values <- f1$fitted.values / sum(f1$fitted.values*s)	

	if (!is.null(out.dir)){
                file <- paste(out.dir,'/','f1_fitted.pdf', sep='')
                pdf(file)
                        f1.hist <- wtd.hist(d.pool[idx], weight=w[idx], breaks=breaks, freq=FALSE)
			lines(f1.hist$mids, f1$fitted.values, col='blue')
                dev.off()
	}

	rl <- list()
	rl$x <- f1.hist$mids
	rl$y <- f1$fitted.values
	rl$breaks <- f1.hist$breaks	
	rl
}






