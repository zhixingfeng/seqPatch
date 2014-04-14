get.nullGenomeF <- function(genomeF.wga, genomeSeq, left=10, right=3)
{
	IPD <- genomeF.wga$features$ipd_pos[[1]]
	mol.id <- genomeF.wga$features$moleculeID_pos[[1]]
	genome_start <- genomeF.wga$genome.start.pos[[1]]

	cat('get homologous locus:\n')	
	max_index <- .Call('R_API_getNullDist_ID',genomeSeq$pos[[1]], as.integer(left), as.integer(right))
	cat('get genomeF:\n')
	rl <- .Call('R_API_getNullDist', IPD, mol.id, as.integer(genome_start), max_index)

	genomeF.wga$features$ipd_pos[[1]] <- rl$IPD
	genomeF.wga$features$ipd_neg[[1]] <- rl$IPD
	genomeF.wga$features$moleculeID_pos[[1]] <- rl$mol_id
	genomeF.wga$features$moleculeID_neg[[1]] <- rl$mol_id
	
	genomeF.wga
}



collapse.rl <- function(rl, item = 'prop')
{
	rl.list <- list()
	k <- 1
	for (i in 1:length(rl$pos)) {
		if (is.null(rl$pos[[i]][[item]])) stop('unable to find ', item)
		rl.list[[k]] <- rl$pos[[i]][[item]]
		k <- k + 1
	}
	for (i in 1:length(rl$neg)) {
		if (is.null(rl$neg[[i]][[item]])) stop('unable to find ', item)
                rl.list[[k]] <- rl$neg[[i]][[item]]
                k <- k + 1
        }
	unlist(rl.list)
}

saturation.analysis <- function(genomeF.native, genomeF.wga, rate=seq(0.1,1,0.1), out.dir=NULL, max.iter=100)
{
	f1.B <- list()
	rl.B <- list()
	for (i in 1:length(rate)){
		genomeF.sub <- sample_subreads(genomeF.native, rate[i])	
		if (!is.null(out.dir)){
			cur.out.dir <- paste(out.dir,'/rate_', rate[i], sep='')
			f1.B[[i]] <- estimate.f1(genomeF.sub, genomeF.wga, out.dir = cur.out.dir)
		}else{
			f1.B[[i]] <- estimate.f1(genomeF.sub, genomeF.wga, out.dir = NULL)
		}
		rl.B[[i]] <- detectModPropEB(genomeF.sub, genomeF.wga, f1.B[[i]]$x, f1.B[[i]]$y, max.iter)
	}
	
	names(f1.B) <- paste('rate_',rate, sep='')
	names(rl.B) <- paste('rate_',rate, sep='')

	# generate qqplot
	if (!is.null(out.dir)){
		cur.out.dir <-  paste(out.dir,'/qqplot' , sep='')		
		if (!file.exists(cur.out.dir)) 
			dir.create(cur.out.dir, recursive = TRUE)
		prop.last <- collapse.rl(rl.B[[length(rate)]], item = 'prop')
		avg.n.last <- collapse.rl(rl.B[[length(rate)]], item = 'avg_n')
		for (i in 1:(length(rate)-1)){
			prop.i <- collapse.rl(rl.B[[i]], item='prop')
			avg.n.i <- collapse.rl(rl.B[[i]], item = 'avg_n')
			file <- paste(cur.out.dir,'/rate_',rate[i],'_vs_',rate[length(rate)],'.pdf', sep='')
			pdf(file)
			qqplot(prop.i, prop.last, xlim=c(0,1), ylim=c(0,1), 
				xlab=paste('rate = ', rate[i],'\n','average number of subreads = ', round(mean(avg.n.i),2), sep=''),
				ylab=paste('rate = ', rate[length(rate)], '\n','average number of subreads = ', round(mean(avg.n.last),2), sep=''), main='')
			abline(0,1,col='blue', lwd=2)
			dev.off()			
		}		
	}
		
	list(f1.B=f1.B, rl.B=rl.B)
	
}


sample_subreads <- function(genomeF, rate)
{
	if (rate<0 | rate>1) rate<-1
	if (is.null(genomeF$features$moleculeID_pos) | is.null(genomeF$features$moleculeID_neg))
                stop('can not find molecule ID in genomeF.native')

	# forward strand
        cat('sampling forward strand, rate is', rate,'\n')
        for (i in 1:length(genomeF$features$ipd_pos)){
		cur.ref <- names(genomeF$features$ipd_pos[i])
		cat(cur.ref,'\n')
		for (j in 1:length(genomeF$features$ipd_pos[[i]])){
			IPD <- genomeF$features$ipd_pos[[i]][[j]]
			idx <- genomeF$features$moleculeID_pos[[i]][[j]]
			cur.rl <- .Call('R_API_sample_subreads', IPD, idx, rate)
			if (length(IPD)==0){
				genomeF$features$ipd_pos[[i]][[j]] <- numeric(0)
				genomeF$features$moleculeID_pos[[i]][[j]] <- numeric(0)
			}else{
				genomeF$features$ipd_pos[[i]][[j]] <- cur.rl$IPD
				genomeF$features$moleculeID_pos[[i]][[j]] <- cur.rl$idx	
			}
			if (j==1000*floor(j/1000)) cat(j,'\r')
		}
		cat(length(genomeF$features$ipd_pos[[i]]),'\n')
        }

	# backward strand
        cat('sampling backward strand, rate is', rate,'\n')
        for (i in 1:length(genomeF$features$ipd_neg)){
                cur.ref <- names(genomeF$features$ipd_neg[i])
                cat(cur.ref,'\n')
                for (j in 1:length(genomeF$features$ipd_neg[[i]])){
                        IPD <- genomeF$features$ipd_neg[[i]][[j]]
                        idx <- genomeF$features$moleculeID_neg[[i]][[j]]
                        cur.rl <- .Call('R_API_sample_subreads', IPD, idx, rate)
                        if (length(IPD)==0){
                                genomeF$features$ipd_neg[[i]][[j]] <- numeric(0)
                                genomeF$features$moleculeID_neg[[i]][[j]] <- numeric(0)
                        }else{
				genomeF$features$ipd_neg[[i]][[j]] <- cur.rl$IPD
                        	genomeF$features$moleculeID_neg[[i]][[j]] <- cur.rl$idx
			}
                        if (j==1000*floor(j/1000)) cat(j,'\r')
                }
                cat(length(genomeF$features$ipd_neg[[i]]),'\n')
        }

	genomeF	
}





detectModPropEB <- function(genomeF.native, genomeF.wga, f1.x = NULL, f1.y = NULL, is_f1_var=FALSE, max.iter = 100, is.EM=FALSE)
{
	if (is.null(genomeF.native$features$moleculeID_pos) | is.null(genomeF.native$features$moleculeID_neg))
		stop('can not find molecule ID in genomeF.native')
	if ((is.null(f1.x) | is.null(f1.y)) & is.EM == FALSE)
		stop('f1.x and f1.y should not be NULl when is.EM=FALSE')
	s <- diff(f1.x)[1]
        f1.int <- sum(f1.y*s)
        if (abs(f1.int - 1)>=1e-6) stop('f1 is not a valid density function')
        if (any(f1.x != sort(f1.x))) stop('f1.x should be sorted')

	
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
					as.numeric(f1.x), as.numeric(f1.y), as.integer(is_f1_var), as.integer(max.iter), as.integer(is.EM))	
			
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
                                        as.numeric(f1.x), as.numeric(f1.y),  as.integer(is_f1_var), as.integer(max.iter), as.integer(is.EM))

        }	

	rl
}

EB.mixture.model.EM.naive <- function(IPD, idx, mu.0 = 0, sigma.0 = 1, max.iter = 100)
{

        .Call('R_API_EBmixture_EM_naive', as.numeric(IPD), as.numeric(idx), as.numeric(mu.0), as.numeric(sigma.0), as.integer(max.iter))
}


EB.mixture.model.EM <- function(IPD, idx, mu.0 = 0, sigma.0 = 1, f1.x, f1.y, max.iter = 100)
{
        s <- diff(f1.x)[1]
        f1.int <- sum(f1.y*s)
        if (abs(f1.int - 1)>=1e-6) stop('f1 is not a valid density function')
        if (any(f1.x != sort(f1.x))) stop('f1.x should be sorted')

        .Call('R_API_EBmixture_EM', as.numeric(IPD), as.numeric(idx), as.numeric(mu.0), as.numeric(sigma.0),
                 as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))
}


EB.mixture.model <- function(IPD, idx, mu.0 = 0, sigma.0 = 1, f1.x, f1.y, max.iter = 100)
{
	s <- diff(f1.x)[1]
	f1.int <- sum(f1.y*s)
	if (abs(f1.int - 1)>=1e-6) stop('f1 is not a valid density function')
	if (any(f1.x != sort(f1.x))) stop('f1.x should be sorted')

	.Call('R_API_EBmixture', as.numeric(IPD), as.numeric(idx), as.numeric(mu.0), as.numeric(sigma.0),
		 as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))
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



estimate.f1 <- function(genomeF.native, genomeF.wga, nulltype = 1, locfdr.df = 28, locfdr.pct0 = 1/4, locfdr.bre = 120, 
			f1.df = 28, f1.bre = 120, out.dir = NULL)
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
	rm(z.score.pool);gc()
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
	w[w<=1e-3 | zz<0] <- 0
	rm(rl.fdr);gc()
	
	idx.w <- which(w>1e-3)
	idx <- idx.sel[idx.w]
	
	d.pool.idx <- d.pool[idx]
	#d.pool.idx <- round(d.pool.idx, 6)	
	w.idx <- w[idx.w]
	#w.idx <- round(w.idx,6)
	rm(d.pool);gc()
	rm(w);gc()
	
	
	breaks <- seq(min(d.pool.idx), max(d.pool.idx), length=f1.bre)
	#breaks[2:(length(breaks)-1)] <- round(breaks[2:(length(breaks)-1)], 6)
	if (!is.null(out.dir)){
                file <- paste(out.dir,'/','f1_d_hist.pdf', sep='')
                pdf(file)
			f1.hist <- m.wtd.hist(d.pool.idx, weight=w.idx, breaks=breaks)
			#f1.hist <- weighted.hist(d.pool.idx, w.idx, breaks=breaks,plot=TRUE)
                dev.off()
        }else{
		f1.hist <- m.wtd.hist(d.pool.idx, weight=w.idx, breaks=breaks, plot=FALSE)
		#f1.hist <- weighted.hist(d.pool.idx, w.idx, breaks=breaks,plot=FALSE)
        }
	
	# fit f1 
	#data <- as.data.frame(cbind(f1.hist$counts, f1.hist$mids))
	#names(data) <- c('y', 'x')
	#f1 <- glm(y~ns(x,df=f1.df), quasipoisson, data=data)		
	
	f1 <- list()
	f1$fitted.values <- f1.hist$counts
	s <- diff(breaks)[1]	
	f1$fitted.values <- f1$fitted.values / sum(f1$fitted.values*s)	

	if (!is.null(out.dir)){
                file <- paste(out.dir,'/','f1_density.pdf', sep='')
                pdf(file)
                        f1.hist <- m.wtd.hist(d.pool.idx, weight=w.idx, breaks=breaks, freq=FALSE)
			#f1.hist <- weighted.hist(d.pool.idx, w.idx, breaks=breaks, freq=FALSE,plot=TRUE)
			lines(f1.hist$mids, f1$fitted.values, col='blue')
                dev.off()
	}

	rl <- list()
	rl$x <- f1.hist$mids
	rl$y <- f1$fitted.values
	rl$breaks <- f1.hist$breaks	
	rl
}






